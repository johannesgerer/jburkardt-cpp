# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "test_opt_con.hpp";

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_OPT_CON_PRB.
//
//  Discussion:
//
//    TEST_OPT_CON_PRB calls the TEST_OPT_CON tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_OPT_CON_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_OPT_CON library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_OPT_CON_PRB\n";
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
//    18 October 2011
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
  cout << "  Problem    Title\n";
  cout << "\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "  " << setw(6) << problem
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
//    TEST02 evaluates the objective function at each starting point.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *f;
  double *fs;
  int i;
  int know;
  int m;
  int n = 100000;
  int problem;
  int problem_num;
  int seed;
  string title;
  double *x;
  double *xs;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For each problem, evaluate the function at many points.\n";
  cout << "  Number of sample points = " << n << "\n";
//
//  Get the number of problems.
//
  problem_num = p00_problem_num ( );

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    cout << "\n";
    cout << "  Problem " << problem << "\n";

    title = p00_title ( problem );

    cout << "  " << title << "\n";

    m = p00_m ( problem );

    cout << "  M =     " << m << "\n";

    a = new double[m];
    b = new double[m];
 
    p00_ab ( problem, m, a, b );

    cout << "\n";
    cout << "    I      A(i)      B(i)\n";
    cout << "\n";

    for ( i = 0; i < m; i++ )
    {
      cout << "  " << setw(4) << i
           << "  " << setw(10) << a[i]
           << "  " << setw(10) << b[i] << "\n";
    }

    seed = 123456789;
    x = r8col_uniform_new ( m, n, a, b, &seed );
    f = p00_f ( problem, m, n, x );

    cout << "\n";
    cout << "  Max(F) = " << r8vec_max ( n, f ) << "\n";
    cout << "  Min(F) = " << r8vec_min ( n, f ) << "\n";

    know = 0;
    xs = p00_sol ( problem, m, know );
    if ( know != 0 )
    {
      fs = p00_f ( problem, m, 1, xs );
      cout << "  F(X*)  = " << fs[0] << "\n";
      delete [] fs;
      delete [] xs;
    }
    else
    {
      cout << "  X* is not given.\n";
    }

    delete [] a;
    delete [] b;
    delete [] f;
    delete [] x;
  }
  return;
}
