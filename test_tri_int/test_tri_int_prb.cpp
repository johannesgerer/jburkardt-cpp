# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_tri_int.H"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_TRI_INT_PRB.
//
//  Discussion:
//
//    TEST_TRI_INT_PRB calls the TEST_TRI_INT routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_TRI_INT_PRB:\n";
  cout << "  C++ version,\n";
  cout << "  Test the TEST_TRI_INT library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_TRI_INT_PRB:\n";
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
//    TEST01 tests GET_PROB_NUM and P00_TITLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int problem;
  int prob_num;
  char *title;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  GET_PROB_NUM reports the number of problems.\n";
  cout << "  P00_TITLE returns a title for each problem.\n";

  prob_num = get_problem_num ( );

  cout << "\n";
  cout << "  The number of problems available is " << prob_num << "\n";
  cout << "\n";
  cout << "  The problem titles:\n";
  cout << "\n";

  for ( problem = 1; problem <= prob_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "  " << setw(8) << problem
         << "  " << title << "\n";
    delete [] title;
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests P00_MONTE_CARLO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double abs_error;
  double exact = 1.0;
  int n;
  int n_log;
  int n_log_max = 15;
  int problem;
  int prob_num;
  double result;
  int seed;
  char *title;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_MONTE_CARLO applies a Monte Carlo rule.\n";

  prob_num = get_problem_num ( );

  cout << "\n";
  cout << "  Problem            Exact         Seed\n";
  cout << "           Pts       Approx        Error\n";
  cout << "\n";
//
//  Pick a problem.
//
  for ( problem = 1; problem <= prob_num; problem++ )
  {
    title = p00_title ( problem );

    seed = 123456789;
//
//  Call RANDOM_INITIALIZE in case we are using the FORTRAN90
//  random number generator inside of P00_MONTE_CARLO!
//
    seed = random_initialize ( seed );

    cout << "\n";
    cout << title << "\n";
    cout << "  " << setw(6) << problem
         << "  " << setw(12) << exact
         << "  " << setw(12) << seed << "\n";
    cout << "\n";
//
//  Pick a number of points.
//
    for ( n_log = 0; n_log <= n_log_max; n_log++ )
    {
      n = i4_power ( 2, n_log );

      result = p00_monte_carlo ( problem, n, &seed );

      abs_error = r8_abs ( exact - result );

      cout << "  " << "      "
           << "  " << setw(6) << n
           << "  " << setw(12) << result
           << "  " << setw(12) << abs_error << "\n";
    }
    delete [] title;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests P00_VERTEX_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double abs_error;
  double exact = 1.0;
  int level;
  int level_max = 4;
  int n;
  int problem;
  int prob_num;
  double result;
  int singularity;
  char *title;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  P00_VERTEX_SUB applies a vertex rule with subdivision.\n";

  prob_num = get_problem_num ( );

  cout << "\n";
  cout << "  Problem            Exact\n";
  cout << "           Pts       Approx        Error\n";
  cout << "\n";
//
//  Pick a problem.
//
  for ( problem = 1; problem <= prob_num; problem++ )
  {
    title = p00_title ( problem );
    singularity = p00_singularity ( problem );

    cout << "\n";
    cout << title << "\n";
    cout << "  " << setw(6) << problem
         << "  " << setw(12) << exact << "\n";
    cout << "\n";

    if ( singularity == 1 )
    {
      cout << "  Skip this problem, it has vertex singularities.\n";
    }
    else if ( singularity == 2 )
    {
      cout << "  Skip this problem, it has edge singularities.\n";
    }
    else if ( singularity == 3 )
    {
      cout << "  Skip this problem, it has internal singularities.\n";
    }
    else
    {
//
//  Pick a number of points.
//
      n = 0;
      result = 0.0;

      for ( level = 0; level <= 4; level++ )
      {
        p00_vertex_sub ( problem, level, &n, &result );

        abs_error = r8_abs ( exact - result );

        cout << "  " << "      "
             << "  " << setw(6) << n
             << "  " << setw(12) << result
             << "  " << setw(12) << abs_error << "\n";
      }
    }
    delete [] title;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests P00_WANDZURA05_SUB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  double abs_error;
  double exact = 1.0;
  int level;
  int n;
  int problem;
  int prob_num;
  double result;
  int test;
  int test_max = 5;
  char *title;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  P00_WANDZURA05_SUB applies a Wandzura rule with subdivision.\n";

  prob_num = get_problem_num ( );

  cout << "\n";
  cout << "  Problem            Exact\n";
  cout << "           Pts       Approx        Error\n";
  cout << "\n";
//
//  Pick a problem.
//
  for ( problem = 1; problem <= prob_num; problem++ )
  {
    title = p00_title ( problem );

    cout << "\n";
    cout << title << "\n";
    cout << "  " << setw(6) << problem
         << "  " << setw(12) << exact << "\n";
    cout << "\n";
//
//  Pick a number of points.
//
    for ( test = 0; test <= test_max; test++ )
    {
      level = i4_power ( 2, test );

      result = p00_wandzura05_sub ( problem, level, &n );

      abs_error = r8_abs ( exact - result );

      cout << "  " << "    "
           << "  " << setw(6) << n
           << "  " << setw(12) << result
           << "  " << setw(12) << abs_error << "\n";

      if ( abs_error < 1.0E-06 )
      {
        cout << "                            Accuracy acceptable\n";
        break;
      }
    }
    delete [] title;
  }
  return;
}
