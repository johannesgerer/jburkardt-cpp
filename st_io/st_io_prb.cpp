# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "st_io.hpp"

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
//    MAIN is the main program for ST_IO_PRB.
//
//  Discussion:
//
//    ST_IO_PRB tests the ST_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ST_IO_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ST_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ST_IO_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void  test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests ST_WRITE.
//
//  Discussion:
//
//    The matrix is:
//
//      11  12   0   0  15
//      21  22   0   0   0
//       0   0  33   0  35
//       0   0   0  44   0
//      51   0  53   0  55
//
//    The index vectors are 1 based, and so have to be converted to
//    0-base before being written.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double ast[11] = {
    51.0, 12.0, 11.0, 33.0, 15.0, 
    53.0, 55.0, 22.0, 35.0, 44.0, 
    21.0 };
  int i_max;
  int i_min;
  int ist[11] = {
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 };
  int j_max;
  int j_min;
  int jst[11] = {
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 };
  int m = 5;
  int n = 5;
  int nst = 11;
  string output_filename = "a5by5.st";

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ST_WRITE writes an ST file.\n";

  i4vec_dec ( nst, ist );
  i4vec_dec ( nst, jst );

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  st_print ( m, n, nst, ist, jst, ast, 
    "  Sparse Triple (ST) data:" );

  st_write ( output_filename, m, n, nst, ist, jst, ast );

  cout << "\n";
  cout << "  Wrote the matrix data to '" << output_filename << "'\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ST_HEADER_READ, ST_DATA_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *ast;
  int i_max;
  int i_min;
  string input_filename = "kershaw.st";
  int *ist;
  int j_max;
  int j_min;
  int *jst;
  int m;
  int n;
  int nst;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  ST_HEADER_READ reads the header from an ST file.\n";
  cout << "  ST_DATA_READ reads the data from an ST file.\n";
  cout << "\n";
  cout << "  Read the data from '" << input_filename << "'\n";

  st_header_read ( input_filename, i_min, i_max, j_min, j_max, m, n, nst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  ast = new double[nst];
  ist = new int[nst];
  jst = new int[nst];

  st_data_read ( input_filename, m, n, nst, ist, jst, ast );

  st_print ( m, n, nst, ist, jst, ast, 
    "  Sparse Triplet (ST) data read from file:" );

  delete [] ast;
  delete [] ist;
  delete [] jst;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests ST_SORT_A.
//
//  Discussion:
//
//    The matrix is:
//
//      11  12   0   0  15
//      21  22   0   0   0
//       0   0  33   0  35
//       0   0   0  44   0
//      51   0  53   0  55
//
//    The index vectors are 1 based, and so have to be converted to
//    0-base before being written.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double ast[11] = {
    51.0, 12.0, 11.0, 33.0, 15.0, 
    53.0, 55.0, 22.0, 35.0, 44.0, 
    21.0 };
  int i_max;
  int i_min;
  int ist[11] = {
     5, 1, 1, 3, 1, 5, 5, 2, 3, 4, 2 };
  int j_max;
  int j_min;
  int jst[11] = {
     1, 2, 1, 3, 5, 3, 5, 2, 5, 4, 1 };
  int m = 5;
  int n = 5;
  int nst = 11;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  ST_SORT_A sorts an ST matrix by columns.\n";

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );

  st_print ( m, n, nst, ist, jst, ast, "  Matrix data before sorting:" );

  st_sort_a ( m, n, nst, ist, jst, ast );

  st_print ( m, n, nst, ist, jst, ast, "  Matrix data sorted by column:" );

  return;
}
