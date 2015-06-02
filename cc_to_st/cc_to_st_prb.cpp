# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cc_to_st.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CC_TO_ST_PRB.
//
//  Discussion:
//
//    CC_TO_ST_PRB tests the CC_TO_ST library.
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
  timestamp ( );
  cout << "\n";
  cout << "CC_TO_ST_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CC_TO_ST library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CC_TO_ST_PRB\n";
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
//    TEST01 tests CC_TO_ST using a 1-based matrix.
//
//  Discussion:
//
//    This test uses a trivial matrix whose full representation is:
//
//          2  3  0  0  0
//          3  0  4  0  6
//      A = 0 -1 -3  2  0
//          0  0  1  0  0
//          0  4  2  0  1
//
//    The 1-based CC representation is
//
//      #  ICC  CCC  ACC
//     --  ---  ---  ---
//      1    1    1    2
//      2    2         3
//
//      3    1    3    3
//      4    3        -1
//      5    5         4
//
//      6    2    6    4
//      7    3        -3
//      8    4         1
//      9    5         2
//
//     10    3   10    2
//
//     11    2   11    6
//     12    5         1
//
//     13    *   13
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
# define N 5
# define NCC 12

  double acc[NCC] = {
    2.0,  3.0, 
    3.0, -1.0,  4.0, 
    4.0, -3.0,  1.0, 2.0, 
    2.0, 
    6.0, 1.0 };
  double *ast;
  int ccc[N+1] = {
    1, 3, 6, 10, 11, 13 };
  int i;
  int icc[NCC] = {
    1, 2, 
    1, 3, 5, 
    2, 3, 4, 5, 
    3, 
    2, 5 };
  int *ist;
  int *jst;
  int m = 5;
  int n = N;
  int ncc = NCC;
  int nst;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Convert a 1-based CC matrix to ST format.\n";
//
//  Print the CC matrix.
//
  cc_print ( m, n, ncc, icc, ccc, acc, "  The CC matrix:" );
//
//  Convert it.
//
  ist = new int[ncc];
  jst = new int[ncc];
  ast = new double[ncc];

  cc_to_st ( m, n, ncc, icc, ccc, acc, nst, ist, jst, ast );
//
//  Print the ST matrix.
//
  st_print ( m, n, nst, ist, jst, ast, "  The ST matrix:" );
//
//  Free memory.
//
  delete [] ast;
  delete [] ist;
  delete [] jst;

  return;
# undef N
# undef NCC
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CC_TO_ST using a 0-based matrix.
//
//  Discussion:
//
//    This test uses a trivial matrix whose full representation is:
//
//          2  3  0  0  0
//          3  0  4  0  6
//      A = 0 -1 -3  2  0
//          0  0  1  0  0
//          0  4  2  0  1
//
//    The 0-based CC representation is
//
//      #  ICC  CCC  ACC
//     --  ---  ---  ---
//      0    0    0    2
//      1    1         3
//
//      2    0    2    3
//      3    2        -1
//      4    4         4
//
//      5    1    5    4
//      6    2        -3
//      7    3         1
//      8    4         2
//
//      9    2    9    2
//
//     10    1   10    6
//     11    4         1
//
//     12    *   12
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
# define N 5
# define NCC 12

  double acc[NCC] = {
    2.0,  3.0, 
    3.0, -1.0,  4.0, 
    4.0, -3.0,  1.0, 2.0, 
    2.0, 
    6.0, 1.0 };
  double *ast;
  int ccc[N+1] = {
    0, 2, 5, 9, 10, 12 };
  int i;
  int icc[NCC] = {
    0, 1, 
    0, 2, 4, 
    1, 2, 3, 4, 
    2, 
    1, 4 };
  int *ist;
  int *jst;
  int m = 5;
  int n = N;
  int ncc = NCC;
  int nst;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Convert a 0-based CC matrix to ST format.\n";
//
//  Print the CC matrix.
//
  cc_print ( m, n, ncc, icc, ccc, acc, "  The CC matrix:" );
//
//  Convert it.
//
  ist = new int[ncc];
  jst = new int[ncc];
  ast = new double[ncc];

  cc_to_st ( m, n, ncc, icc, ccc, acc, nst, ist, jst, ast );
//
//  Print the ST matrix.
//
  st_print ( m, n, nst, ist, jst, ast, "  The ST matrix:" );
//
//  Free memory.
//
  delete [] ast;
  delete [] ist;
  delete [] jst;

  return;
}
