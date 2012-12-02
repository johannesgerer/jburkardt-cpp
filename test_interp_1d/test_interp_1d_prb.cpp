# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_interp_1d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( int nd );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_INTERP_1D_TEST tests the TEST_INTERP_1D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int nd;

  timestamp ( );

  cout << "\n";
  cout << "TEST_INTERP_1D_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INTERP_1D library.\n";
  cout << "  The R8LIB library is needed.\n";

  test01 ( );

  nd = 11;
  test02 ( nd );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INTERP_1D_TEST\n";
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
//    TEST01 simply prints the title of each function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
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
  cout << "  Print the title of each function.\n";

  prob_num = p00_prob_num ( );
  
  cout << "\n";
  cout << "  There are " << prob_num << " functions available:\n";
  cout << "  Index  Title\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    title = p00_title ( prob );
    cout << "  " << setw(2) << prob << "  " << title << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( int nd )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_INTERP_1D_TEST02 evaluates each test function at ND sample points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ND, the number of sample points.
//
{
  double a;
  double b;
  double *f;
  int j;
  int prob;
  int prob_num;
  double *x;

  cout << "\n";
  cout << "TEST_INTERP_1D_TEST02\n";
  cout << "  Use P00_F to sample each function.\n";

  prob_num = p00_prob_num ( );

  a = 0.0;
  b = 1.0;
  x = r8vec_linspace_new ( nd, a, b );

  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    f = p00_f ( prob, nd, x );
    cout << "\n";
    cout << "  X, F(X) for problem " << prob << "\n";
    cout << "\n";
    for ( j = 0; j < nd; j++ )
    {
      cout << "  " << setw(2) << j
           << "  " << setw(10) << x[j]
           << "  " << setw(10) << f[j] << "\n";
    }
    delete [] f;
  }

  delete [] x;

  return;
}
