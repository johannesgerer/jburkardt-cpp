# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_interp_nd.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( int m, int n );
void test03 ( int m, int n );
void test04 ( int m, int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_INTERP_ND_PRB.
//
//  Discussion:
//
//    TEST_INTERP_ND_PRB calls the TEST_INTERP_ND tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;

  timestamp (  );

  cout << "\n";
  cout << "TEST_INTERP_ND_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INTERP_ND library.\n";
  cout << "  The R8LIB library is also needed.\n";

  test01 ( );

  n = 10;
  for ( m = 2; m <= 4; m++ )
  {
    test02 ( m, n );
  }

  n = 2;
  for ( m = 2; m <= 4; m++ )
  {
    test03 ( m, n );
  }

  m = 4;
  n = 10000;
  test04 ( m, n );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INTERP_ND_PRB\n";
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
//    TEST01 prints the title of each test function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;
  string title;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  P00_TITLE returns the problem title.\n";

  prob_num = p00_prob_num ( );
  cout << "\n";
  cout << "  There are a total of " << prob_num << " problems.\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );
    cout << "  " << prob << "  \"" << title << "\"\n";
  }

  return;
}
//****************************************************************************80

void test02 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 samples each function in M dimensions, at N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
{
  double *c;
  double *f;
  int i;
  int j;
  int prob;
  int prob_num;
  int seed;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_F evaluates the function.\n";
  cout << "  Here, we use spatial dimension M = " << m << "\n";
  cout << "  The number of points is N = " << n << "\n";

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem  " << prob << "\n";

    c = p00_c ( prob, m, seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, seed );
    r8vec_print ( m, w, "  W parameters:" );

    cout << "\n";
    cout << "      F(X)              X(1)      X(2) ...\n";
    cout << "\n";

    f = p00_f ( prob, m, c, w, n, x );

    for ( j = 0; j < n; j++ )
    {
      cout << "  " << setw(14) << f[j];
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << setw(10) << x[i+j*m];
      }
      cout << "\n";
    }
    delete [] c;
    delete [] f;
    delete [] w;
  }
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 samples each derivative component in M dimensions, at N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
{
  double *c;
  double *d;
  double *dmat;
  double *f;
  int i;
  int id;
  int j;
  int prob;
  int prob_num;
  int seed;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  P00_D evaluates derivative components.\n";
  cout << "  Here, we use spatial dimension M = " << m << "\n";
  cout << "  The number of points is N = " << n << "\n";

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem  " << prob << "\n";

    c = p00_c ( prob, m, seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, seed );
    r8vec_print ( m, w, "  W parameters:" );

    cout << "\n";
    cout << "                         X(1)      X(2) ...\n";
    cout << "      F(X)            dFdX(1)   dFdX(2) ...\n";

    f = p00_f ( prob, m, c, w, n, x );

    dmat = new double[m*n];

    for ( id = 0; id < m; id++ )
    {
      d = p00_d ( prob, m, id, c, w, n, x );
      for ( j = 0; j < n; j++ )
      {
        dmat[id+j*m] = d[j];
      }
      delete [] d;
    }

    for ( j = 0; j < n; j++ )
    {
      cout << "\n";
      cout << "                ";
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << setw(10) << x[i+j*m];
      }
      cout << "\n";
      cout << "  " << setw(14) << f[j];
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << setw(10) << dmat[i+j*m];
      }
      cout << "\n";
    }
    delete [] c;
    delete [] dmat;
    delete [] f;
    delete [] w;
  }

  delete [] x;

  return;
}
//****************************************************************************80

void test04 ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 estimates integrals in M dimensions, using N points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
{
  double *c;
  double *f;
  int j;
  int prob;
  int prob_num;
  double q1;
  double q2;
  int seed;
  double *w;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "P00_Q returns the integral of F over [0,1]^m.\n";
  cout << "  Here, we use spatial dimension M = " << m << "\n";
  cout << "  The number of sample points is N = " << n << "\n";

  seed = 123456789;
  x = r8mat_uniform_01_new ( m, n, seed );

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem  " << prob << "\n";

    c = p00_c ( prob, m, seed );
    r8vec_print ( m, c, "  C parameters:" );

    w = p00_w ( prob, m, seed );
    r8vec_print ( m, w, "  W parameters:" );

    cout << "\n";
    cout << "      Exact Integral     Q\n";
    cout << "\n";

    q1 = p00_q ( prob, m, c, w );

    f = p00_f ( prob, m, c, w, n, x );
    q2 = r8vec_sum ( n, f ) / ( double ) ( n );

    cout << "  " << setw(14) << q1
         << "  " << setw(14) << q2 << "\n";

    delete [] c;
    delete [] f;
    delete [] w;
  }

  delete [] x;

  return;
}
