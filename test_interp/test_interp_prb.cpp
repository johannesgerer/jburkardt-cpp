# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "test_interp.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_INTERP_PRB.
//
//  Discussion:
//
//    TEST_INTERP_PRB tests the TEST_INTERP library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp (  );
  cout << "\n";
  cout << "TEST_INTERP_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INTERP library.\n";
  cout << "  This test also requires the R8LIB library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INTERP_PRB\n";
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
//    TEST01 shows how P00_STORY can be called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  P00_STORY prints the problem \"story\".\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem " << prob << "\n";

    p00_story ( prob );
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 prints the data for each problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2011
//
//  Author:
//
//    John Burkardt
//
{
  int data_num;
  int dim_num;
  double *p;
  int prob;
  int prob_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  P00_DATA_NUM returns N, the number of data points.\n";
  cout << "  P00_DIM_NUM returns M, the dimension of data.\n";
  cout << "  P00_DATA returns the actual (MxN) data.\n";

  prob_num = p00_prob_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    cout << "\n";
    cout << "  Problem  " << prob << "\n";

    data_num = p00_data_num ( prob );
    cout << "  DATA_NUM " << data_num << "\n";
    dim_num = p00_dim_num ( prob );
    cout << "  DIM_NUM  " << dim_num << "\n";

    p = p00_data ( prob, dim_num, data_num );

    r8mat_transpose_print ( dim_num, data_num, p, "  Data array:" );

    delete [] p;
  }

  return;
}
