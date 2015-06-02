# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "box_behnken.hpp"

int main ( );

void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BOX_BEHNKEN_PRB.
//
//  Discussion:
//
//    BOX_BEHNKEN_PRB tests the BOX_BEHNKEN library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 October 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BOX_BEHNKEN_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BOX_BEHNKEN library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BOX_BEHNKEN_PRB\n";
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
//    TEST01 tests BOX_BEHNKEN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3

  int dim_num = DIM_NUM;
  double range[DIM_NUM*2] = {
    0.0, 10.0,  5.0,
    1.0, 11.0, 15.0 };
  int x_num;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  BOX_BEHNKEN computes a Box-Behnken dataset.\n";

  r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

  x_num = box_behnken_size ( dim_num );

  cout << "\n";
  cout << "  For dimension DIM_NUM = " << dim_num << "\n";
  cout << "  the Box-Behnken design is of size " << x_num << "\n";

  x = box_behnken ( dim_num, x_num, range );

  r8mat_transpose_print ( dim_num, x_num, x, "  The Box-Behnken design:" );

  delete [] x;

  return;
# undef DIM_NUM
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R8MAT_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2012
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 4

  int dim_num = DIM_NUM;
  string file_out_name = "box_behnken_04_33.txt";
  double range[DIM_NUM*2] = {
    0.0, 0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0 };
  int x_num;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R8MAT_WRITE writes a Box-Behnken dataset\n";
  cout << "  to a file.\n";

  r8mat_transpose_print ( dim_num, 2, range, "  The ranges:" );

  x_num = box_behnken_size ( dim_num );

  cout << "\n";
  cout << "  For dimension DIM_NUM = " << dim_num << "\n";
  cout << "  the Box-Behnken design is of size " << x_num << "\n";

  x = box_behnken ( dim_num, x_num, range );

  r8mat_write ( file_out_name, dim_num, x_num, x );

  delete [] x;

  cout << "\n";
  cout << "  The data was written to the file \"" << file_out_name << "\".\n";

  return;
# undef DIM_NUM
}

