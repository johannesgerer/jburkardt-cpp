# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "st_to_cc.hpp"

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
//    MAIN is the main program for ST_TO_CC_PRB.
//
//  Discussion:
//
//    ST_TO_CC_PRB tests the ST_TO_CC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ST_TO_CC_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the ST_TO_CC library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ST_TO_CC_PRB\n";
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
//    TEST01 tests ST_TO_CC using a tiny matrix.
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
//    A (1-based) ST representation, reading in order by rows is:
//
//      I  J   A
//     -- --  --
//      1  1   2
//      1  2   3
//
//      2  1   3
//      2  3   4
//      2  5   6
//
//      3  2  -1
//      3  3  -3
//      3  4   2
//
//      4  3   1
//
//      5  2   4
//      5  3   2
//      5  5   1
//
//    The CC representation (which goes in order by columns) is
//
//      #   I  JC   A
//     --  --  --  --
//      1   1   1   2
//      2   2       3
//
//      3   1   3   3
//      4   3      -1
//      5   5       4
//
//      6   2   6   4
//      7   3      -3
//      8   4       1
//      9   5       2
//
//     10   3  10   2
//
//     11   2  11   6
//     12   5       1
//
//     13   *  13
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
# define NST 12

  double *acc;
  double ast[NST] = {
    2.0,  3.0, 
    3.0,  4.0,  6.0, 
   -1.0, -3.0,  2.0, 
    1.0, 
    4.0,  2.0,  1.0 };
  int *ccc;
  int i_max;
  int i_min;
  int *icc;
  int ist[NST] = {
    1, 1, 
    2, 2, 2, 
    3, 3, 3, 
    4, 
    5, 5, 5 };
  int j_max;
  int j_min;
  int jst[NST] = {
    1, 2, 
    1, 3, 5, 
    2, 3, 4, 
    3, 
    2, 3, 5 };
  int m = 5;
  int n = 5;
  int ncc;
  int nst = NST;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Convert a sparse matrix from ST to CC format.\n";
  cout << "  ST: sparse triplet,    I, J,  A.\n";
  cout << "  CC: compressed column, I, CC, A.\n";

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );
//
//  Decrement the 1-based data.
//
  i4vec_dec ( nst, ist );
  i4vec_dec ( nst, jst );
//
//  Print the ST matrix.
//
  st_print ( m, n, nst, ist, jst, ast, "  The matrix in ST format:" );
//
//  Get the CC size.
//
  ncc = st_to_cc_size ( nst, ist, jst );

  cout << "\n";
  cout << "  Number of CC values = " << ncc << "\n";
//
//  Create the CC indices.
//
  icc = new int[ncc];
  ccc = new int[n+1];

  st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc );
//
//  Create the CC values.
//
  acc = st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc );
//
//  Print the CC matrix.
//
  cc_print ( m, n, ncc, icc, ccc, acc, "  CC Matrix:" );
//
//  Free memory.
//
  delete [] acc;
  delete [] ccc;
  delete [] icc;

  return;
# undef NST
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests ST_TO_CC on a matrix stored in a file.
//
//  Discussion:
//
//    We assume no prior knowledge about the matrix except the filename.
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
  double *acc;
  double *ast;
  int *ccc;
  string filename_st = "west_st.txt";
  int i_max;
  int i_min;
  int *icc;
  int *ist;
  int j_max;
  int j_min;
  int *jst;
  int m;
  int n;
  int ncc;
  int nst;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Convert a sparse matrix from ST to CC format.\n";
  cout << "  ST: sparse triplet,    I, J,  A.\n";
  cout << "  CC: compressed column, I, CC, A.\n";
  cout << "  This matrix is read from the file '" << filename_st << "'\n";
//
//  Get the size of the ST matrix.
//
  st_header_read ( filename_st, i_min, i_max, j_min, j_max, m, n, nst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );
//
//  Allocate space.
//
  ist = new int[nst];
  jst = new int[nst];
  ast = new double[nst];
//
//  Read the ST matrix.
//
  st_data_read ( filename_st, m, n, nst, ist, jst, ast );
//
//  Decrement the 1-based data.
//
  i4vec_dec ( nst, ist );
  i4vec_dec ( nst, jst );
//
//  Print the ST matrix.
//
  st_print ( m, n, nst, ist, jst, ast, "  The matrix in ST format:" );
//
//  Get the CC size.
//
  ncc = st_to_cc_size ( nst, ist, jst );

  cout << "\n";
  cout << "  Number of CC values = " << ncc << "\n";
//
//  Create the CC indices.
//
  icc = new int[ncc];
  ccc = new int[n+1];

  st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc );
//
//  Create the CC values.
//
  acc = st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc );
//
//  Print the CC matrix.
//
  cc_print ( m, n, ncc, icc, ccc, acc, "  CC Matrix:" );
//
//  Free memory.
//
  delete [] acc;
  delete [] ast;
  delete [] ccc;
  delete [] icc;
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
//    TEST03 creates a CC sparse matrix file from an ST file.
//
//  Discussion:
//
//    We assume no prior knowledge about the matrix except the filename.
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
  double *acc;
  double *ast;
  int *ccc;
  string filename_acc = "west_acc.txt";
  string filename_ccc = "west_ccc.txt";
  string filename_icc = "west_icc.txt";
  string filename_st = "west_st.txt";
  int i_max;
  int i_min;
  int *icc;
  int *ist;
  int j_max;
  int j_min;
  int *jst;
  int m;
  int n;
  int ncc;
  int nst;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Convert a sparse matrix from ST to CC format.\n";
  cout << "  ST: sparse triplet,    I, J,  A.\n";
  cout << "  CC: compressed column, I, CC, A.\n";
  cout << "  The ST matrix is read from the file '" << filename_st << "'\n";
  cout << "  and the CC matrix is written to the files:\n";
  cout << "    '" << filename_icc << "',\n";
  cout << "    '" << filename_ccc << "', and\n";
  cout << "    '" << filename_acc << "'.\n";
//
//  Get the size of the ST matrix.
//
  st_header_read ( filename_st, i_min, i_max, j_min, j_max, m, n, nst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );
//
//  Allocate space.
//
  ist = new int[nst];
  jst = new int[nst];
  ast = new double[nst];
//
//  Read the ST matrix.
//
  st_data_read ( filename_st, m, n, nst, ist, jst, ast );
//
//  Decrement the 1-based data.
//
  i4vec_dec ( nst, ist );
  i4vec_dec ( nst, jst );
//
//  Get the CC size.
//
  ncc = st_to_cc_size ( nst, ist, jst );

  cout << "\n";
  cout << "  Number of CC values = " << ncc << "\n";
//
//  Create the CC indices.
//
  icc = new int[ncc];
  ccc = new int[n+1];

  st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc );
//
//  Create the CC values.
//
  acc = st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc );
//
//  Write the CC matrix.
//
  i4vec_write ( filename_icc, ncc, icc );
  i4vec_write ( filename_ccc, n + 1, ccc );
  r8vec_write ( filename_acc, ncc, acc );
//
//  Free memory.
//
  delete [] acc;
  delete [] ast;
  delete [] ccc;
  delete [] icc;
  delete [] ist;
  delete [] jst;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 works with a CC sparse matrix with many repeated index pairs.
//
//  Discussion:
//
//    To complete this test, I want to compare AST * X and ACC * X.
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
  double *acc;
  double *ast;
  double *b1;
  double *b2;
  int *ccc;
  int i_max;
  int i_min;
  int *icc;
  int *ist;
  int j_max;
  int j_min;
  int *jst;
  int m;
  int n;
  int ncc;
  int nst;
  int nx;
  int ny;
  double r;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Convert a sparse matrix from ST to CC format.\n";
  cout << "  ST: sparse triplet,    I, J,  A.\n";
  cout << "  CC: compressed column, I, CC, A.\n";
  cout << "  The ST matrix is the Wathen finite element matrix.\n";
  cout << "  It has many repeated index pairs.\n";
  cout << "  To check, compare ACC*X - AST*X for a random X.\n";
//
//  Get the size of the ST matrix.
//
  nx = 3;
  ny = 3;
  nst = wathen_st_size ( nx, ny );

  cout << "\n";
  cout << "  Number of ST values = " << nst << "\n";
//
//  Set the formal matrix size
//
  m = 3 * nx * ny + 2 * nx + 2 * ny + 1;
  n = m;
//
//  Set a random vector.
//
  seed = 123456789;
  x = r8vec_uniform_01_new ( n, seed );
//
//  Allocate space.
//
  ist = new int[nst];
  jst = new int[nst];
//
//  Create the ST matrix.
//
  seed = 123456789;
  ast = wathen_st ( nx, ny, nst, seed, ist, jst );

  i_min = i4vec_min ( nst, ist );
  i_max = i4vec_max ( nst, ist );
  j_min = i4vec_min ( nst, jst );
  j_max = i4vec_max ( nst, jst );

  st_header_print ( i_min, i_max, j_min, j_max, m, n, nst );
//
//  Compute B1 = AST * X
//
  b1 = st_mv ( m, n, nst, ist, jst, ast, x );
//
//  Get the CC size.
//
  ncc = st_to_cc_size ( nst, ist, jst );

  cout << "  Number of CC values = " << ncc << "\n";
//
//  Create the CC indices.
//
  icc = new int[ncc];
  ccc = new int[n+1];
  st_to_cc_index ( nst, ist, jst, ncc, n, icc, ccc );
//
//  Create the CC values.
//
  acc = st_to_cc_values ( nst, ist, jst, ast, ncc, n, icc, ccc );
//
//  Compute B2 = ACC * X.
//
  b2 = cc_mv ( m, n, ncc, icc, ccc, acc, x );
//
//  Compare B1 and B2.
//
  r = r8vec_diff_norm ( m, b1, b2 );
  cout << "  || ACC*X - AST*X|| = " << r << "\n";
//
//  Free memory.
//
  delete [] acc;
  delete [] ast;
  delete [] b1;
  delete [] b2;
  delete [] ccc;
  delete [] icc;
  delete [] ist;
  delete [] jst;
  delete [] x;

  return;
}
