# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_ls.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_LS_PRB.
//
//  Discussion:
//
//    TEST_LS_PRB tests the TEST_LS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TEST_LS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_LS library.\n";
  cout << "  This test also requires the R8LIB library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_LS_PRB\n";
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
//    TEST01 summarizes the test data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double b_norm;
  int i;
  int j;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r;
  double r_norm;
  double *x;
  double x_norm;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Get each least squares test and compute the maximum residual.\n";
  cout << "  The L2 norm of the residual MUST be no greater than\n";
  cout << "  the L2 norm of the right hand side, else 0 is a better solution.\n";

  prob_num = p00_prob_num ( );

  cout << "\n";
  cout << "  Number of problems = " << prob_num << "\n";
  cout << "\n";
  cout << "  Index     M     N     ||B||         ||X||         ||R||\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    m = p00_m ( prob );
    n = p00_n ( prob );

    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x = p00_x ( prob, n );

    r = new double[m];

    for ( i = 0; i < m; i++ )
    {
      r[i] = - b[i];
      for ( j = 0; j < n; j++ )
      {
        r[i] = r[i] + a[i+j*m] * x[j];
      }
    }

    b_norm = r8vec_norm ( m, b );
    x_norm = r8vec_norm ( n, x );
    r_norm = r8vec_norm ( m, r );

    cout << "  " << setw(5) << prob
         << "  " << setw(4) << m
         << "  " << setw(4) << n
         << "  " << setw(12) << b_norm
         << "  " << setw(12) << x_norm
         << "  " << setw(12) << r_norm << "\n";

    delete [] a;
    delete [] b;
    delete [] r;
    delete [] x;
  }
  return;
}
