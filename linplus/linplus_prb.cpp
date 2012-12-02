# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "linplus.hpp"

int main ( );
void test019 ( );
void test0192 ( );
void test0193 ( );
void test0194 ( );
void test0195 ( );
void test0196 ( );
void test02 ( );
void test03 ( );
void test035 ( );
void test037 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test105 ( );
void test11 ( );
void test115 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test151 ( );
void test152 ( );
void test153 ( );
void test154 ( );
void test155 ( );
void test1565 ( );
void test1566 ( );
void test157 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );
void test193 ( );
void test195 ( );
void test197 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test235 ( );
void test24 ( );
void test25 ( );
void test26 ( );
void test265 ( );
void test2655 ( );
void test27 ( );
void test275 ( );
void test28 ( );
void test285 ( );
void test29 ( );
void test295 ( );
void test30 ( );
void test31 ( );
void test315 ( );
void test317 ( );
void test32 ( );
void test33 ( );
void test34 ( );
void test345 ( );
void test35 ( );
void test36 ( );
void test37 ( );
void test38 ( );
void test385 ( );
void test39 ( );
void test40 ( );
void test41 ( );
void test42 ( );
void test422 ( );
void test423 ( );
void test425 ( );
void test426 ( );
void test428 ( );
void test43 ( );
void test44 ( );
void test443  ( );
void test445 ( );
void test45 ( );
void test46 ( );
void test47 ( );
void test48 ( );
void test485 ( );
void test49 ( );
void test50 ( );
void test505 ( );
void test51 ( );
void test515 ( );
void test517 ( );
void test52 ( );
void test525 ( );
void test527 ( );
void test53 ( );
void test534 ( );
void test535 ( );
void test54 ( );
void test55 ( );
void test555 ( );
void test56 ( );
void test57 ( );
void test5705 ( );
void test571 ( );
void test572 ( );
void test5722 ( );
void test5724 ( );
void test5725 ( );
void test573 ( );
void test574 ( );
void test5745 ( );
void test575 ( );
void test577 ( );
void test58 ( );
void test581 ( );
void test583 ( );
void test585 ( );
void test587 ( );
void test589 ( );
void test59 ( );
void test60 ( );
void test605 ( );
void test61 ( );
void test62 ( );
void test63 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    LINPLUS_PRB tests routines from the LINPLUS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "LINPLUS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LINPLUS library.\n";

  test019 ( );
  test0192 ( );
  test0193 ( );
  test0194 ( );
  test0195 ( );
  test0196 ( );
  test02 ( );
  test03 ( );
  test035 ( );
  test037 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test105 ( );
  test11 ( );
  test115 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test154 ( );
  test155 ( );
  test1565 ( );
  test1566 ( );
  test157 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );
  test193 ( );
  test195 ( );
  test197 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test235 ( );
  test24 ( );
  test25 ( );
  test26 ( );
  test265 ( );
  test2655 ( );
  test27 ( );
  test275 ( );
  test28 ( );
  test285 ( );
  test29 ( );
  test295 ( );

  test30 ( );
  test31 ( );
  test315 ( );
  test317 ( );
  test32 ( );
  test33 ( );
  test34 ( );
  test345 ( );
  test35 ( );
  test36 ( );
  test37 ( );
  test38 ( );
  test39 ( );
  test385 ( );
  test39 ( );

  test40 ( );
  test41 ( );
  test42 ( );
  test422 ( );
  test423 ( );
  test425 ( );
  test426 ( );
  test428 ( );
  test43 ( );
  test44 ( );
  test443( );
  test445 ( );
  test45 ( );

  test46 ( );
  test47 ( );
  test48 ( );
  test485 ( );
  test49 ( );

  test50 ( );
  test505 ( );
  test51 ( );
  test515 ( );
  test517 ( );
  test52 ( );
  test525 ( );
  test527 ( );
  test53 ( );
  test534 ( );
  test535 ( );
  test54 ( );
  test55 ( );
  test555 ( );
  test56 ( );
  test57 ( );
  test5705 ( );
  test571 ( );
  test572 ( );
  test5722 ( );
  test5724 ( );
  test5725 ( );
  test573 ( );
  test574 ( );
  test5745 ( );
  test575 ( );
  test577 ( );
  test58 ( );
  test581 ( );
  test583 ( );
  test585 ( );
  test587 ( );
  test589 ( );
  test59 ( );

  test60 ( );
  test605 ( );
  test61 ( );
  test62 ( );
  test63 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LINPLUS_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test019 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST019 tests C8VEC_UNITY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N_MAX 5

  int n;
  double *x;

  cout << "\n";
  cout << "TEST019\n";
  cout << "  C8VEC_UNITY returns the N complex roots of unity.\n";
  cout << "\n";

  for ( n = 1; n <= N_MAX; n++ )
  {
    x = c8vec_unity ( n );

    c8vec_print ( n, x, "  Roots of unity:" );
    delete [] x;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void test0192 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0192 tests R8CC_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  string a_file = "r8cc_a.txt";
  int col[N+1] = { 1, 4, 6, 8, 10, 13 };
  string col_file = "r8cc_col.txt";
  int row[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  string row_file = "r8cc_row.txt";

  cout << "\n";
  cout << "TEST0192\n";
  cout << "  For a matrix in the R8CC format,\n";
  cout << "  (double precision compressed column sparse)\n";
  cout << "  R8CC_WRITE writes the matrix to 3 files.\n";
  cout << "\n";
  cout << "  Matrix rows M     = " << M << "\n";
  cout << "  Matrix columns N  = " << N << "\n";
  cout << "  Nonzeros NZ_NUM   = " << NZ_NUM << "\n";

  i4vec_print ( N+1, col, "  The COL vector:" );

  i4vec_print ( NZ_NUM, row, "  The ROW vector:" );

  a = r8cc_indicator ( M, N, NZ_NUM, col, row );

  r8cc_print ( M, N, NZ_NUM, col, row, a, "  The R8CC matrix:" );

  r8cc_write ( col_file, row_file, a_file, M, N, NZ_NUM, col, row, a );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test0193 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0193 tests R8CC_READ_SIZE and R8CC_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  string a_file = "r8cc_a.txt";
  int base;
  int *col;
  string col_file = "r8cc_col.txt";
  int m;
  int n;
  int nz_num;
  int *row;
  string row_file = "r8cc_row.txt";

  cout << "\n";
  cout << "TEST0193\n";
  cout << "  For a matrix in the R8CC format,\n";
  cout << "  (double precision compressed column sparse)\n";
  cout << "  R8CC_READ_SIZE reads the sizes of the data;\n";
  cout << "  R8CC_READ reads the data.\n";

  r8cc_read_size ( col_file, row_file, &m, &n, &nz_num, &base );

  cout << "\n";
  cout << "  Matrix rows M     = " << m      << "\n";
  cout << "  Matrix columns N  = " << n      << "\n";
  cout << "  Nonzeros NZ_NUM   = " << nz_num << "\n";
  cout << "  Index base (0/1)  = " << base   << "\n";

  a = new double[nz_num];
  col = new int[n+1];
  row = new int[nz_num];

  r8cc_read ( col_file, row_file, a_file, m, n, nz_num, col, row, a );

  i4vec_print ( n+1, col, "  The COL vector:" );

  i4vec_print ( nz_num, row, "  The ROW vector:" );

  r8cc_print ( m, n, nz_num, col, row, a, "  The R8CC matrix:" );

  delete [] a;
  delete [] col;
  delete [] row;
//
//  Delete the files.  
//  If you want to look at the files, comment out these lines
//  and rerun the problem!
//
  file_delete ( a_file );
  file_delete ( col_file );
  file_delete ( row_file );

  return;
}
//****************************************************************************80

void test0194 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0194 tests R8VEC_TO_R8CB, R8CB_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int m = 5;
  int n = 8;
  int ml = 2;
  int mu = 1;

  double *a;
  int i;
  int j;
  int k;
  double *x;

  cout << "\n";
  cout << "TEST0194\n";
  cout << "  For a compressed banded matrix,\n";
  cout << "  R8VEC_TO_R8CB converts a real vector to an R8CB matrix.\n";
  cout << "  R8CB_TO_R8VEC converts an R8CB matrix to a real vector.\n";
  cout << "\n";
  cout << "  Matrix rows M =       " << m  << "\n";
  cout << "  Matrix columns N =    " << n  << "\n";
  cout << "  Lower bandwidth ML  = " << ml << "\n";
  cout << "  Upper bandwidth MU  = " << mu << "\n";

  a = r8cb_indicator ( m, n, ml, mu );

  r8cb_print ( m, n, ml, mu, a, "  The R8CB indicator matrix:" );

  x = r8cb_to_r8vec ( m, n, ml, mu, a );

  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= ml+mu+1; i++ )
    {
      k = k + 1;
      cout << setw(6)  << i      << "  "
           << setw(6)  << j      << "  "
           << setw(6)  << k      << "  "
           << setw(12) << x[k-1] << "\n";
    }
  }

  delete [] a;
  a = r8vec_to_r8cb ( m, n, ml, mu, x );

  r8cb_print ( m, n, ml, mu, a, "  The recovered R8CB indicator matrix:" );

  delete [] a;
  delete [] x;

  return;
}
//****************************************************************************80

void test0195 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0195 tests R8VEC_TO_R8GB, R8GB_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int m = 5;
  int n = 8;
  int ml = 2;
  int mu = 1;

  double *a;
  int i;
  int j;
  int k;
  double *x;

  cout << "\n";
  cout << "TEST0195\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8VEC_TO_R8GB converts a real vector to an R8GB matrix.\n";
  cout << "  R8GB_TO_R8VEC converts an R8GB matrix to a real vector.\n";
  cout << "\n";
  cout << "  Matrix rows M =       " << m  << "\n";
  cout << "  Matrix columns N =    " << n  << "\n";
  cout << "  Lower bandwidth ML  = " << ml << "\n";
  cout << "  Upper bandwidth MU  = " << mu << "\n";

  a = r8gb_indicator ( m, n, ml, mu );

  r8gb_print ( m, n, ml, mu, a, "  The R8GB indicator matrix:" );

  x = r8gb_to_r8vec ( m, n, ml, mu, a );

  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= 2*ml+mu+1; i++ )
    {
      k = k + 1;
      cout << setw(6)  << i      << "  "
           << setw(6)  << j      << "  "
           << setw(6)  << k      << "  "
           << setw(12) << x[k-1] << "\n";
    }
  }

  delete [] a;
  a = r8vec_to_r8gb ( m, n, ml, mu, x );

  r8gb_print ( m, n, ml, mu, a, "  The recovered R8GB indicator matrix:" );

  delete [] a;
  delete [] x;

  return;
}
//****************************************************************************80

void test0196 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0196 tests R8VEC_TO_R8GE, R8GE_TO_R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int m = 4;
  int n = 6;

  double *a;
  int i;
  int j;
  int k;
  double *x;

  cout << "\n";
  cout << "TEST0196\n";
  cout << "  For a general matrix,\n";
  cout << "  R8VEC_TO_R8GE converts a real vector to an R8GE matrix.\n";
  cout << "  R8GE_TO_R8VEC converts an R8GE matrix to a real vector.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << m << "\n";
  cout << "  Matrix columns N = " << n << "\n";

  a = r8ge_indicator ( m, n );

  r8ge_print ( m, n, a, " The R8GE indicator matrix:" );

  x = r8ge_to_r8vec ( m, n, a );

  k = 0;
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      k = k + 1;
      cout << setw(6) << i      << "  "
           << setw(6) << j      << "  "
           << setw(6) << k      << "  "
           << setw(6) << x[k-1] << "\n";
    }
  }

  delete [] a;
  a = r8vec_to_r8ge ( m, n, x );

  r8ge_print ( m, n, a, "  The recovered R8GE indicator matrix:" );

  delete [] a;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R83_CR_FA, R83_CR_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[3*N];
  double *a_cr;
  double b[N];
  int i;
  int j;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R83_CR_FA factors a real tridiagonal matrix;\n";
  cout << "  R83_CR_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Demonstrate multiple system solution method.\n";
//
//  Set the matrix values.
//
  a[0+0*3] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < N-1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(N-1)*3] = 0.0;
//
//  Factor the matrix once.
//
  a_cr = r83_cr_fa ( N, a );

  for ( j = 1; j <= 2; j++ )
  {
    cout << "\n";
    cout << "  Solve linear system number #" << j << ".\n";

    if ( j == 1 )
    {
      for ( i = 0; i < N-1; i++ )
      {
        b[i] = 0.0;
      }
      b[N-1] = ( double ) ( N + 1 );
    }
    else
    {
      b[0] = 1.0;
      for ( i = 1; i < N-1; i++ )
      {
        b[i] = 0.0;
      }
      b[N-1] = 1.0;
    }
//
//  Solve the linear system.
//
    x = r83_cr_sl ( N, a_cr, b );

    r8vec_print_some ( N, x, 10, "  Solution:" );

    delete [] x;
  }

  delete [] a_cr;

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R83_CR_FA, R83_CR_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[3*N];
  double *a_cr;
  double *b;
  int i;
  int j;
  double *x;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  For a real tridiagonal matrix,\n";
  cout << "  using CYCLIC REDUCTION,\n";
  cout << "  R83_CR_FA factors;\n";
  cout << "  R83_CR_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  The matrix is NOT symmetric.\n";
//
//  Set the matrix values.
//
  a[0+0*3] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*3] = ( double ) ( j + 1 );
  }

  for ( j = 0; j < N; j++ )
  {
    a[1+j*3] = 4.0 * ( double ) ( j + 1 );
  }

  for ( j = 0; j < N-1; j++ )
  {
    a[2+j*3] = ( double ) ( j + 1 );
  }
  a[2+(N-1)*3] = 0.0;

  r83_print ( N, a, "  The matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r83_mxv ( N, a, x );
//
//  Factor the matrix.
//
  a_cr = r83_cr_fa ( N, a );
//
//  Solve the linear system.
//
  x = r83_cr_sl ( N, a_cr, b );

  r8vec_print_some ( N, x, 10, "  Solution:" );

  delete [] a_cr;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 tests R83_GS_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  double a[3*N];
  double *b;
  int i;
  int j;
  int job;
  int maxit = 1000;
  double *x;

  cout << "\n";
  cout << "TEST035\n";
  cout << "  For a real tridiagonal system,\n";
  cout << "  R83_GS_SL solves a linear system using Gauss-Seidel iteration.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Iterations per call = " << maxit << "\n";
  cout << "\n";
//
//  Set the matrix values.
//
  a[0+0*3] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < N-1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(N-1)*3] = 0.0;

  for ( job = 0; job <= 1; job++ )
  {
    if ( job == 0 )
    {
      cout << "\n";
      cout << "  Solving A * x = b.\n";
      cout << "\n";
    }
    else
    {
      cout << "\n";
      cout << "  Solving A' * x = b.\n";
      cout << "\n";
    }
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r83_mxv ( N, a, x );
    }
    else
    {
      b = r83_vxm ( N, a, x );
    }
//
//  Set the starting solution.
//
    for ( i = 0; i < N; i++ )
    {
      x[i] = 0.0;
    }
//
//  Solve the linear system.
//
    for ( i = 1; i <= 3; i++ )
    {
      r83_gs_sl ( N, a, b, x, maxit, job );

      r8vec_print_some ( N, x, 10, "  Current estimated solution:" );
    }

    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests R83_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;

  cout << "\n";
  cout << "TEST037\n";
  cout << "  R83_INDICATOR sets up an S3 indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r83_indicator ( N );

  r83_print ( N, a, "  The S3 indicator matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests R83_JAC_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 100

  double a[3*N];
  double *b;
  int i;
  int j;
  int job;
  int maxit = 1000;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  For a real tridiagonal system,\n";
  cout << "  R83_JAC_SL solves a linear system using Jacobi iteration.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Iterations per call = " << maxit << "\n";
  cout << "\n";
//
//  Set the matrix values.
//
  a[0+0*3] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*3] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*3] = 2.0;
  }
  for ( j = 0; j < N-1; j++ )
  {
    a[2+j*3] = -1.0;
  }
  a[2+(N-1)*3] = 0.0;

  for ( job = 0; job <= 1; job++ )
  {
    if ( job == 0 )
    {
      cout << "\n";
      cout << "  Solving A * x = b.\n";
      cout << "\n";
    }
    else
    {
      cout << "\n";
      cout << "  Solving A' * x = b.\n";
      cout << "\n";
    }
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r83_mxv ( N, a, x );
    }
    else
    {
      b = r83_vxm ( N, a, x );
    }

    r8vec_print_some ( N, b, 10, "  The right hand side:" );
//
//  Set the starting solution.
//
    for ( i = 0; i < N; i++ )
    {
      x[i] = 0.0;
    }
//
//  Solve the linear system.
//
    for ( i = 1; i <= 3; i++ )
    {
      r83_jac_sl ( N, a, b, x, maxit, job );

      r8vec_print_some ( N, x, 10, "  Current estimated solution:" );
    }

    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests R83_NP_DET, R83_NP_FA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  double det;
  int info;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a tridiagonal matrix that can be factored\n";
  cout << "    with no pivoting,\n";
  cout << "  R83_NP_FA factors,\n";
  cout << "  R83_NP_DET computes the determinant.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r83_random ( N, &seed );

  r83_print ( N, a, "  The S3 matrix:" );
//
//  Copy the matrix into general storage.
//
  b = r83_to_r8ge ( N, a );
//
//  Factor the matrix.
//
  info = r83_np_fa ( N, a );

  r83_print ( N, a, "  The factored S3 matrix:" );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST05 - Warning!\n";
    cout << "  R83_NP_FA returns INFO = " << info << "\n";
  }
//
//  Compute the determinant.
//
  det = r83_np_det ( N, a );

  cout << "\n";
  cout << "  R83_NP_DET computes determinant =  " << det << "\n";
//
//  Factor the matrix in R8GE storage.
//
  info = r8ge_np_fa ( N, b );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST05 - Warning!\n";
    cout << "  R8GE_NP_FA returns INFO = " << info << "\n";
  }
//
//  Compute the determinant of the matrix in R8GE storage.
//
  det = r8ge_np_det ( N, b );

  cout << "  R8GE_NP_DET computes determinant = " << det << "\n";

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
// TEST06 tests R83_NP_FA, R83_NP_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  For a tridiagonal matrix that can be factored\n";
  cout << "    with no pivoting,\n";
  cout << "  R83_NP_FA factors;\n";
  cout << "  R83_NP_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r83_random ( N, &seed );

  r83_print ( N, a, "  The tridiagonal matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r83_mxv ( N, a, x );
  delete [] x;
//
//  Factor the matrix.
//
  info = r83_np_fa ( N, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST06 - Fatal error!\n";
    cout << "  The test matrix is singular.\n";
    return;
  }
//
//  Solve the linear system.
//
  job = 0;
  x = r83_np_sl ( N, a, b, job );
 
  r8vec_print ( N, x, "  Solution to A*x=b:" );
//
//  Set the desired solution
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side, using the factored matrix.
//
  job = 1;
  delete [] b;
  b = r83_np_ml ( N, a, x, job );
//
//  Solve the linear system.
//
  job = 1;
  delete [] x;
  x = r83_np_sl ( N, a, b, job );
 
  r8vec_print ( N, x, "  Solution to A'*x=b:" );
 
  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests R83_NP_FS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  R83_NP_FS factors and solves a tridiagonal\n";
  cout << "    linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix elements.
//
  a = r83_random ( N, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute b = A * x.
//
  b = r83_mxv ( N, a, x );
//
//  Wipe out the solution.
//
  delete [] x;
//
//  Solve the system.
//
  x = r83_np_fs ( N, a, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests R83_NP_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  R83_NP_ML computes A*x or A'*x\n";
  cout << "    where A has been factored by R83_FA.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r83_random ( N, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r83_mxv ( N, a, x );
    }
    else
    {
      b = r83_vxm ( N, a, x );
    }
//
//  Factor the matrix.
//
    info = r83_np_fa ( N, a );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST08 - Fatal error!\n";
      cout << "  R83_NP_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      return;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r83_np_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x:" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] a;
    delete [] b;
    delete [] b2;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests R83P_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double *a;
  double *b;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  R83P_DET, determinant of a tridiagonal\n";
  cout << "    periodic matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r83p_random ( N, &seed );

  r83p_print ( N, a, "  The periodic tridiagonal matrix:" );
//
//  Copy the matrix into a general array.
//
  b = r83p_to_r8ge ( N, a );
//
//  Factor the matrix.
//
  info = r83p_fa ( N, a, work2, work3, &work4 );
//
//  Compute the determinant.
//
  det = r83p_det ( N, a, work4 );

  cout << "\n";
  cout << "  R83P_DET computes the determinant = " << det << "\n";
//
//  Factor the general matrix.
//
  info = r8ge_fa ( N, b, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, b, pivot );

  cout << "  R8GE_DET computes the determinant = " << det << "\n";

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests R83P_FA, R83P_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;
  double *x;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  R83P_FA factors a tridiagonal periodic system.\n";
  cout << "  R83P_SL solves a tridiagonal periodic system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r83p_random ( N, &seed );
//
//  Factor the matrix.
//
  info = r83p_fa ( N, a, work2, work3, &work4 );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  TEST10 - Fatal error!\n";
    cout << "  R83P_FA returns INFO = " << info << "\n";
    return;
  }

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    b = r83p_ml ( N, a, x, job );
//
//  Solve the linear system.
//
    delete [] x;

    x = r83p_sl ( N, a, b, job, work2, work3, work4 );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'*x=b:" );
    }

    delete [] x;
    delete [] b;
  }
 
  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests R83P_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  R83P_INDICATOR sets up an S3P indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r83p_indicator ( N );

  r83p_print ( N, a, "  The S3P indicator matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests R83P_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double work2[N-1];
  double work3[N-1];
  double work4;
  double *x;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  R83P_ML computes A*x or A'*X\n";
  cout << "    where A has been factored by R83P_FA.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r83p_random ( N, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r83p_mxv ( N, a, x );
    }
    else
    {
      b = r83p_vxm ( N, a, x );
    }
//
//  Factor the matrix.
//
    info = r83p_fa ( N, a, work2, work3, &work4 );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST11 - Fatal error!\n";
      cout << "  R83P_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      return;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r83p_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] a;
    delete [] x;
    delete [] b;
    delete [] b2;
  }

  return;
# undef N
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests R85_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;

  cout << "\n";
  cout << "TEST115\n";
  cout << "  R85_INDICATOR sets up an R85 indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r85_indicator ( N );

  r85_print ( N, a, "  The R85 indicator matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests R85_NP_FS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  R85_NP_FS factors and solves a pentadiagonal\n";
  cout << "    linear system, with no pivoting.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix to a random value.
//
  a = r85_random ( N, &seed );

  r85_print ( N, a, "  The pentadiagonal matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute b = A * x.
//
  b = r85_mxv ( N, a, x );

  r8vec_print ( N, b, "  Right hand side:" );
//
//  Wipe out the solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = 0.0;
  }
//
//  Solve the system.
//
  delete [] x;
  x = r85_np_fs ( N, a, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests R8BB_FA, R8BB_PRINT, R8BB_RANDOM, R8BB_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 8
# define N2 2
# define ML 1
# define MU 1
# define N N1+N2

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  For a border banded matrix:\n";
  cout << "  R8BB_FA factors;\n";
  cout << "  R8BB_PRINT prints;\n";
  cout << "  R8BB_RANDOM randomizes;\n";
  cout << "  R8BB_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N  << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8bb_random ( N1, N2, ML, MU, &seed );

  r8bb_print ( N1, N2, ML, MU, a, "  The border-banded matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );

  r8vec_print ( N, b, "  The right hand side vector:" );
//
//  Factor the matrix.
//
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );

  r8bb_print ( N1, N2, ML, MU, a, "  The FACTORED border-banded matrix:" );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST13 - Fatal error!\n";
    cout << "  R8BB_FA claims the matrix is singular.\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  delete [] x;
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests R8BB_FA, R8BB_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 0
# define N2 10
# define N N1+N2
# define ML 0
# define MU 0

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  For a border banded matrix:\n";
  cout << "  R8BB_FA factors;\n";
  cout << "  R8BB_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N  << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8bb_random ( N1, N2, ML, MU, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );
//
//  Factor the matrix.
//
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST14 - Fatal error!\n";
    cout << "  R8BB_FA claims the matrix is singular.\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  delete [] x;
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );
 
  delete [] a;
  delete [] b; 
  delete [] x;

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests R8BB_FA, R8BB_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 10
# define N2 0
# define N N1+N2
# define ML 1
# define MU 1

  double *a;
  double *b;
  int i;
  int info;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For a border banded matrix:\n";
  cout << "  R8BB_FA factors;\n";
  cout << "  R8BB_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N  << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8bb_random ( N1, N2, ML, MU, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8bb_mxv ( N1, N2, ML, MU, a, x );
//
//  Factor the matrix
//
  info = r8bb_fa ( N1, N2, ML, MU, a, pivot );
 
  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST15 - Fatal error!\n";
    cout << "  R8BB_FA claims the matrix is singular.\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  delete [] x;
  x = r8bb_sl ( N1, N2, ML, MU, a, pivot, b );

  r8vec_print ( N, x, "  Solution to A*x=b:" );
 
  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef ML
# undef MU
# undef N
# undef N1
# undef N2
}
//****************************************************************************80

void test151 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST151 tests R8BB_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 6
# define N2 2
# define ML 1
# define MU 1

  double *a;

  cout << "\n";
  cout << "TEST151\n";
  cout << "  R8BB_INDICATOR sets up an R8BB indicator matrix.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N1 + N2 << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8bb_indicator ( N1, N2, ML, MU );

  r8bb_print ( N1, N2, ML, MU, a, "  The R8BB indicator matrix:" );

  delete [] a;

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
//****************************************************************************80

void test152 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST152 tests R8BLT_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int n = 6;
  int ml = 2;

  cout << "\n";
  cout << "TEST152\n";
  cout << "  R8BLT_INDICATOR sets up an R8BLT indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  Lower bandwidth ML = " << ml << "\n";

  a = r8blt_indicator ( n, ml );

  r8blt_print ( n, ml, a, "  The R8BLT indicator matrix:" );

  delete [] a;

  return;
}
//****************************************************************************80

void test153 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST153 tests R8BLT_MXV, R8BLT_PRINT, R8BLT_RANDOM, R8BLT_SL, R8BLT_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ML 3
# define N 10

  double *a;
  double *b;
  int i;
  int j;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST153\n";
  cout << "  For a band matrix in lower triangular storage,\n";
  cout << "  R8BLT_RANDOM sets a random value;\n";
  cout << "  R8BLT_SL solves systems;\n";
  cout << "  R8BLT_MXV computes A*x;\n";
  cout << "  R8BLT_VXM computes A'*x;\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";

  a = r8blt_random ( N, ML, &seed );

  r8blt_print ( N, ML, a, "  The R8BLT matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8blt_mxv ( N, ML, a, x );
    }
    else
    {
      b = r8blt_vxm ( N, ML, a, x );
    }

    r8vec_print ( N, b, "  The right hand side:" );
//
//  Solve the linear system.
//
    delete [] x;
    x = r8blt_sl ( N, ML, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'*x=b:" );
    }

    delete [] b;
    delete [] x;
  }

  delete [] a;

  return;
# undef ML
# undef N
}
//****************************************************************************80

void test154 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST154 tests R8BTO_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define L 3
# define M 2

  double *a;

  cout << "\n";
  cout << "TEST154\n";
  cout << "  For a real block Toeplitz matrix,\n";
  cout << "  R8BTO_INDICATOR sets up an indicator matrix.\n";
  cout << "\n";
  cout << "  Block order M =  " << M << "\n";
  cout << "  Block number L = " << L << "\n";
  cout << "  Matrix order N = " << M * L << "\n";

  a = r8bto_indicator ( M, L );

  r8bto_print ( M, L, a, "  The block Toeplitz matrix:" );

  delete [] a;

  return;
# undef L
# undef M
}
//****************************************************************************80

void test155 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST155 tests R8BTO_MXV, R8BTO_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2004
//
//  Author:
//
//    John Burkardt
//
{
# define L 3
# define M 2
# define N ( M * L )

  double a[M*M*(2*L-1)] = {
    1.0, 5.0, 2.0, 5.0, 
    3.0, 6.0, 4.0, 6.0, 
    5.0, 7.0, 6.0, 7.0, 
    7.0, 8.0, 8.0, 8.0, 
    9.0, 9.0, 0.0, 9.0 };
  double *b;
  double *x;

  cout << "\n";
  cout << "TEST155\n";
  cout << "  For a real block Toeplitz matrix,\n";
  cout << "  R8BTO_MXV computes A * x.\n";
  cout << "  R8BTO_VXM computes A'* x.\n";
  cout << "\n";
  cout << "  Block order M =  " << M << "\n";
  cout << "  Block number L = " << L << "\n";
  cout << "  Matrix order N = " << N << "\n";

  r8bto_print ( M, L, a, "  The block Toeplitz matrix:" );

  x = r8ge_indicator ( M, L );

  r8ge_print ( M, L, x, "  The 'vector' x:" );

  b = r8bto_mxv ( M, L, a, x );

  r8ge_print ( M, L, b, "  The product A*x:" );

  delete [] b;
  b = r8bto_vxm ( M, L, a, x );

  r8ge_print ( M, L, b, "  The product A'*x:" );

  delete [] b;
  delete [] x;

  return;
# undef L
# undef M
# undef N
}
//****************************************************************************80

void test1565 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1565 tests R8BUT_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6
# define MU 2

  double *a;

  cout << "\n";
  cout << "TEST1565\n";
  cout << "  R8BUT_INDICATOR sets up an SBUT indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N  << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8but_indicator ( N, MU );

  r8but_print ( N, MU, a, "  The R8BUT indicator matrix:" );

  delete [] a;

  return;
# undef N
# undef MU
}
//****************************************************************************80

void test1566 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1566 tests R8BUT_MXV, R8BUT_PRINT, R8BUT_RANDOM, R8BUT_SL, R8BUT_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MU 3
# define N 10

  double *a;
  double *b;
  int i;
  int j;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST1566\n";
  cout << "  For a band matrix in upper triangular storage,\n";
  cout << "  R8BUT_RANDOM sets a random value;\n";
  cout << "  R8BUT_SL solves systems;\n";
  cout << "  R8BUT_MXV computes matrix-vector products;\n";
  cout << "  R8BUT_VXM computes vector-matrix products;\n";
  cout << "\n";
  cout << "  Matrix order N =     " << N  << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8but_random ( N, MU, &seed );

  r8but_print ( N, MU, a, "  The R8BUT matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8but_mxv ( N, MU, a, x );
    }
    else
    {
      b = r8but_vxm ( N, MU, a, x );
    }

    r8vec_print ( N, b, "  The right hand side:" );
//
//  Solve the linear system.
//
    delete [] x;
    x = r8but_sl ( N, MU, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to the untransposed system:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to the transposed system:" );
    }
    delete [] b;
    delete [] x;

  }

  delete [] a;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test157 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST157 tests R8CB_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define N 10
# define ML 2
# define MU 3

  double *a;

  cout << "\n";
  cout << "TEST157\n";
  cout << "  R8CB_INDICATOR sets up an R8CB indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8cb_indicator ( M, N, ML, MU );

  r8cb_print ( M, N, ML, MU, a, "  The R8CB indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef ML
# undef MU
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests R8CB_NP_FA, R8CB_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define ML 2
# define MU 3

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  For a compact band matrix, no pivoting:\n";
  cout << "  R8CB_NP_FA factors;\n";
  cout << "  R8CB_DET computes the determinant;\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8cb_random ( N, ML, MU, &seed );

  r8cb_print ( N, N, ML, MU, a, "  The compact band matrix:" );
//
//  Copy the matrix into a general array.
//
  a2 = r8cb_to_r8ge ( N, ML, MU, a );
//
//  Factor the matrix.
//
  info = r8cb_np_fa ( N, ML, MU, a );
//
//  Compute the determinant.
//
  det = r8cb_det ( N, ML, MU, a );

  cout << "\n";
  cout << "  R8CB_DET computes the determinant = " << det << "\n";
//
//  Factor the general matrix.
//
  info = r8ge_fa ( N, a2, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, a2, pivot );

  cout << "  R8GE_DET computes the determinant = " << det << "\n";

  delete [] a;
  delete [] a2;

  return;
# undef N
# undef ML
# undef MU
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests R8CB_NP_FA, R8CB_NP_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define ML 1
# define MU 2

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  For a compact band matrix, no pivoting:\n";
  cout << "  R8CB_NP_FA factors;\n";
  cout << "  R8CB_NP_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r8cb_random ( N, ML, MU, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the right hand side.
//
    if ( job == 0 )
    {
      b = r8cb_mxv ( N, ML, MU, a, x );
    }
    else
    {
      b = r8cb_vxm ( N, ML, MU, a, x );
    }
//
//  Factor the matrix.
//
    info = r8cb_np_fa ( N, ML, MU, a );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST17 - Fatal error!\n";
      cout << "  R8CB_NP_FA claims the matrix is singular.\n";
      cout << "  The value of info is " << info << "\n";
      return;
    }
//
//  Solve the system.
//
    delete [] x;

    x = r8cb_np_sl ( N, ML, MU, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'x=b:" );
    }

    delete [] a;
    delete [] b;
    delete [] x;

  }

  return;
# undef N
# undef ML
# undef MU
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests R8CB_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define ML 1
# define MU 2

  double *a;
  double *b;
  double *b2;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  For a compact band matrix:\n";
  cout << "  R8CB_ML computes A*x or A'*X\n";
  cout << "    where A has been factored by R8CB_FA.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r8cb_random ( N, ML, MU, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8cb_mxv ( N, ML, MU, a, x );
    }
    else
    {
      b = r8cb_vxm ( N, ML, MU, a, x );
    }
//
//  Factor the matrix.
//
    info = r8cb_np_fa ( N, ML, MU, a );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST18 - Fatal error!\n";
      cout << "  R8CB_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      return;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r8cb_ml ( N, ML, MU, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] a;
    delete [] b;
    delete [] b2;
    delete [] x;

  }

  return;
# undef N
# undef ML
# undef MU
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests R8CBB_FA, R8CBB_PRINT, R8CBB_RANDOM, R8CBB_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 8
# define N2 2
# define ML 1
# define MU 1

  double *a;
  double *b;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  For a compressed border banded matrix:\n";
  cout << "  R8CBB_RANDOM randomly generates;\n";
  cout << "  R8CBB_PRINT prints;\n";
  cout << "  R8CBB_FA factors (no pivoting);\n";
  cout << "  R8CBB_SL solves.\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N1 + N2  << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8cbb_random ( N1, N2, ML, MU, &seed );

  r8cbb_print ( N1, N2, ML, MU, a, "  The R8CBB matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N1 + N2 );
//
//  Compute the corresponding right hand side.
//
  b = r8cbb_mxv ( N1, N2, ML, MU, a, x );
//
//  Factor the matrix
//
  info = r8cbb_fa ( N1, N2, ML, MU, a );

  r8cbb_print ( N1, N2, ML, MU, a, "  The factored R8CBB matrix:" );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8CBB_FA claims the matrix is singular.\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//
  r8vec_print ( N1 + N2, b, "  The right hand side vector b:" );

  delete [] x;
  x = r8cbb_sl ( N1, N2, ML, MU, a, b );

  r8vec_print ( N1 + N2, x, "  Solution to A*x=b:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
//****************************************************************************80

void test193 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST193 tests R8CBB_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N1 6
# define N2 2
# define ML 1
# define MU 1

  double *a;

  cout << "\n";
  cout << "TEST193\n";
  cout << "  R8CBB_INDICATOR sets up an R8CBB indicator matrix.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Matrix order N     = " << N1 + N2 << "\n";
  cout << "  Matrix suborder N1 = " << N1 << "\n";
  cout << "  Matrix suborder N2 = " << N2 << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8cbb_indicator ( N1, N2, ML, MU );

  r8cbb_print ( N1, N2, ML, MU, a, "  The R8CBB indicator matrix:" );

  delete [] a;

  return;
# undef ML
# undef MU
# undef N1
# undef N2
}
//****************************************************************************80

void test195 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST195 tests R8CI_EVAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *lambda;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST195\n";
  cout << "  R8CI_EVAL finds the eigenvalues of \n";
  cout << "  a real circulant system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ci_random ( N, &seed );

  r8ci_print ( N, a, "  The circulant matrix:" );

  lambda = r8ci_eval ( N, a );

  c8vec_print ( N, lambda, "  The eigenvalues:" );

  delete [] a;
  delete [] lambda;

  return;
# undef N
}
//****************************************************************************80

void test197 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST197 tests R8CI_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST197\n";
  cout << "  R8CI_INDICATOR sets up an R8CI indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8ci_indicator ( N );

  r8ci_print ( N, a, "  The circulant matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests R8CI_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  R8CI_SL solves a circulant system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ci_random ( N, &seed );

  r8ci_print ( N, a, "  The circulant matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8ci_mxv ( N, a, x );
    }
    else
    {
      b = r8ci_vxm ( N, a, x );
    }
//
//  Solve the linear system.
//
    delete [] x;
    x = r8ci_sl ( N, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution to A*x=b:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to A'*x=b:" );
    }

    delete [] b;
    delete [] x;

  }
 
  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests R8GB_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define ML 3
# define MU 2
# define N 10

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_DET computes the determinant.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML  = " << ML << "\n";
  cout << "  Upper bandwidth MU  = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Copy the matrix into a general array.
//
  a2 = r8gb_to_r8ge ( M, N, ML, MU, a );
//
//  Print the matrix just to show what it looks like.
//
  r8ge_print ( M, N, a, "  The banded matrix:" );
//
//  Factor the matrix.
//
  info = r8gb_fa ( N, ML, MU, a, pivot );
//
//  Compute the determinant.
//
  det = r8gb_det ( N, ML, MU, a, pivot );

  cout << "\n";
  cout << "  R8GB_DET computes the determinant = " << det << "\n";
//
//  Factor the general matrix.
//
  info = r8ge_fa ( N, a2, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, a2, pivot );

  cout << "  R8GE_DET computes the determinant = " << det << "\n";

  delete [] a;
  delete [] a2;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests R8GB_FA, R8GB_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define ML 1
# define MU 2
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_FA computes the PLU factors.\n";
  cout << "  R8GB_SL solves a factored linear system.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8gb_random ( M, N, ML, MU, &seed );

  r8gb_print ( M, N, ML, MU, a, "  The banded matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8gb_mxv ( M, N, ML, MU, a, x );
//
//  Factor the matrix.
//
  info = r8gb_fa ( N, ML, MU, a, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8GB_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  job = 0;
  delete [] x;
  x = r8gb_sl ( N, ML, MU, a, pivot, b, job );
 
  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 1;
  delete [] b;
  b = r8gb_ml ( N, ML, MU, a, pivot, x, job );
  r8vec_print ( N, b, "  Right hand side of transposed system:" );
//
//  Solve the linear system.
//
  job = 1;
  delete [] x;
  x = r8gb_sl ( N, ML, MU, a, pivot, b, job );
 
  r8vec_print ( N, x, "  Solution to transposed system:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests R8GB_FA, R8GB_TRF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define ML 1
# define MU 1
# define N 5

  double *a;
  int info;
  int pivot[N];
  int seed;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_FA factors, using LINPACK conventions;\n";
  cout << "  R8GB_TRF factors, using LAPACK conventions;\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  seed = 123456789;
  a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Factor the matrix.
//
  info = r8gb_fa ( N, ML, MU, a, pivot );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB_FA factors:" );
//
//  Set the matrix.
//
  delete [] a;
  seed = 123456789;
  a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Factor the matrix.
//
  info = r8gb_trf ( M, N, ML, MU, a, pivot );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB_TRF factors:" );

  delete [] a;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test235 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST235 tests R8GB_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define ML 3
# define MU 2
# define N 8

  double *a;

  cout << "\n";
  cout << "TEST235\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_INDICATOR sets up an indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8gb_indicator ( M, N, ML, MU );

  r8ge_print ( 2*ML+MU+1, N, a, "  The banded matrix in R8GE format:" );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests R8GB_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define ML 1
# define MU 2
# define N 10

  double *a;
  double *b;
  double *b2;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_ML computes A*x or A'*X\n";
  cout << "    where A has been factored by R8GB_FA.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8gb_mxv ( M, N, ML, MU, a, x );
    }
    else
    {
      b = r8gb_vxm ( M, N, ML, MU, a, x );
    }
//
//  Factor the matrix.
//
    info = r8gb_fa ( N, ML, MU, a, pivot );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST24 - Fatal error!\n";
      cout << "  R8GB_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      return;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r8gb_ml ( N, ML, MU, a, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] a;
    delete [] b;
    delete [] b2;
    delete [] x;
  }

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests R8GB_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define ML 1
# define MU 3
# define N 10

  double *a;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_PRINT prints the matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The banded matrix:" );

  delete [] a;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests R8GB_NZ_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 10
# define ML 1
# define MU 2

  double *a;
  int diag;
  int i;
  int j;
  int nz_num;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_NZ_NUM counts the nonzero entries.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Make some zero entries.
//
  for ( j = 0; j < N; j++ )
  {
    for ( diag = 0; diag < 2*ML+MU+1; diag++ )
    {
      if ( a[diag+j*(2*ML+MU+1)] < 0.3 )
      {
        a[diag+j*(2*ML+MU+1)] = 0.0;
      }
    }
  }

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  nz_num = r8gb_nz_num ( M, N, ML, MU, a );

  cout << "\n";
  cout << "  Nonzero entries = " << nz_num << "\n";

  delete [] a;

  return;
# undef M
# undef N
# undef ML
# undef MU
}
//****************************************************************************80

void test265 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST265 tests R8GB_TO_R8GE and R8GE_TO_R8GB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define ML 2
# define MU 1
# define N 8

  double *a;
  double *b;
  double *c;

  cout << "\n";
  cout << "TEST265\n";
  cout << "  R8GB_TO_R8GE copies a R8GB matrix to a R8GE matrix.\n";
  cout << "  R8GE_TO_R8GB copies a R8GE matrix to a R8GB matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  b = r8gb_to_r8ge ( M, N, ML, MU, a );

  r8ge_print ( M, N, b, "  The R8GE matrix:" );

  c = r8ge_to_r8gb ( M, N, ML, MU, b );

  r8gb_print ( M, N, ML, MU, c, "  The recovered R8GB matrix:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test2655 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST2655 tests R8GB_TO_R8S3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define ML 2
# define MU 1
# define N 8

  double *a;
  double *b;
  int *col;
  int isym;
  int nz_num;
  int *row;

  cout << "\n";
  cout << "TEST2655\n";
  cout << "  R8GB_TO_R8S3 copies a R8GB matrix to a R8S3 matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  a = r8gb_indicator ( M, N, ML, MU );

  r8gb_print ( M, N, ML, MU, a, "  The R8GB matrix:" );

  nz_num = r8gb_nz_num ( M, N, ML, MU, a );

  cout << "  Nonzeros NZ_NUM =    " << nz_num << "\n";

  row = new int[nz_num];
  col = new int[nz_num];
  b = new double[nz_num];

  r8gb_to_r8s3 ( M, N, ML, MU, a, nz_num, &isym, row, col, b );

  r8s3_print ( M, N, nz_num, isym, row, col, b, "  The R8S3 matrix:" );

  delete [] row;
  delete [] col;
  delete [] b;
  delete [] a;

  return;
# undef M
# undef ML
# undef MU
# undef N
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests R8GB_TRF, R8GB_TRS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 10
# define ML 1
# define MU 2
# define NRHS 1

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST27\n";
  cout << "  For a general banded matrix,\n";
  cout << "  R8GB_TRF computes the PLU factors.\n";
  cout << "  R8GB_TRS solves a factored linear system.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Lower bandwidth ML = " << ML << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8gb_random ( M, N, ML, MU, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8gb_mxv ( M, N, ML, MU, a, x );
//
//  Factor the matrix.
//
  info = r8gb_trf ( M, N, ML, MU, a, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST27 - Fatal error!\n";
    cout << "  R8GB_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  delete [] x;
  x = r8gb_trs ( N, ML, MU, NRHS, 'N', a, pivot, b );

  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 1;
  delete [] b;
  b = r8gb_mu ( N, ML, MU, a, pivot, x, job );
//
//  Solve the linear system.
//
  delete [] x;
  x = r8gb_trs ( N, ML, MU, NRHS, 'T', a, pivot, b );

  r8vec_print ( N, x, "  Solution to transposed system:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef M
# undef N
# undef ML
# undef MU
# undef NRHS
}
//****************************************************************************80

void test275 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST275 tests R8GD_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NDIAG 4

  double *a;
  int offset[NDIAG] = { -2, 0, 1, N-1 };

  cout << "\n";
  cout << "TEST275\n";
  cout << "  For a general diagonal matrix:\n";
  cout << "  R8GD_INDICATOR sets up an indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix order N            = " << N << "\n";
  cout << "  Number of diagonals NDIAG = " << NDIAG << "\n";

  i4vec_print ( NDIAG, offset, "  The offset vector:" );

  a = r8gd_indicator ( N, NDIAG, offset );

  r8gd_print ( N, NDIAG, offset, a, "  The R8GD indicator matrix:" );

  delete [] a;

  return;
# undef N
# undef NDIAG
}
//****************************************************************************80

void test28 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST28 tests R8GD_MXV. R8GD_PRINT, R8GD_RANDOM, R8GD_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NDIAG 4

  double *a;
  double *b;
  int i;
  int j;
  int offset[NDIAG] = { -2, 0, 1, N-1 };
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST28\n";
  cout << "  For a general diagonal matrix:\n";
  cout << "  R8GD_MXV computes A * x;\n";
  cout << "  R8GD_PRINT prints it;\n";
  cout << "  R8GD_RANDOM randomly generates one;\n";
  cout << "  R8GD_VXM computes A'*x;\n";
  cout << "\n";
  cout << "  Matrix order N            = " << N << "\n";
  cout << "  Number of diagonals NDIAG = " << NDIAG << "\n";

  i4vec_print ( NDIAG, offset, "  The offset vector:" );
//
//  Set the matrix.
//
  a = r8gd_random ( N, NDIAG, offset, &seed );

  r8ge_print ( N, NDIAG, a, "  The raw matrix: " );

  r8gd_print ( N, NDIAG, offset, a, "  The general diagonal matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8gd_mxv ( N, NDIAG, offset, a, x );

  r8vec_print ( N, b, "  A * x:" );
//
//  Compute the corresponding right hand side.
//
  delete [] b;
  b = r8gd_vxm ( N, NDIAG, offset, a, x );

  r8vec_print ( N, b, "  A' * x:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
# undef NDIAG
}
//****************************************************************************80

void test285 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST285 tests R8GE_CO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  double *a_inverse;
  double a_inverse_norm_l1;
  double a_lu[N*N];
  double a_norm_l1;
  double cond_l1;
  int i;
  int info;
  int j;
  int pivot[N];
  double rcond;
  double row_sum;
  double x = 2.0;
  double y = 3.0;

  cout << "\n";
  cout << "TEST285\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_CO estimates the condition number.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }

  a_norm_l1 = 0.0;
  for ( j = 0; j < N; j++ )
  {
    row_sum = 0.0;
    for ( i = 0; i < N; i++ )
    {
      row_sum = row_sum + r8_abs ( a[i+j*N] );
    }
    a_norm_l1 = r8_max ( a_norm_l1, row_sum );
  }

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_fa ( N, a_lu, pivot );

  a_inverse = r8ge_inverse ( N, a_lu, pivot );

  a_inverse_norm_l1 = 0.0;
  for ( j = 0; j < N; j++ )
  {
    row_sum = 0.0;
    for ( i = 0; i < N; i++ )
    {
      row_sum = row_sum + r8_abs ( a_inverse[i+j*N] );
    }
    a_inverse_norm_l1 = r8_max ( a_inverse_norm_l1, row_sum );
  }

  cond_l1 = a_norm_l1 * a_inverse_norm_l1;

  cout << "\n";
  cout << "  The L1 condition number is " << cond_l1 << "\n";
//
//  Factor the matrix.
//
  rcond = r8ge_co ( N, a, pivot );

  cout << "\n";
  cout << "  The R8GE_CO estimate is     " << 1.0 / rcond << "\n";

  delete [] a_inverse;

  return;
# undef N
}
//****************************************************************************80

void test29 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST29 tests R8GE_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  double det;
  double exact;
  int i;
  int info;
  int j;
  int pivot[N];
  double x = 2.0;
  double y = 3.0;

  cout << "\n";
  cout << "TEST29\n";
  cout << "  R8GE_DET, determinant of a general matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }
//
//  Factor the matrix.
//
  info = r8ge_fa ( N, a, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, a, pivot );
  exact = pow ( x, N-1 ) * ( x + ( ( double ) N ) * y );

  cout << "\n";
  cout << "  R8GE_DET computes the determinant = " << det   << "\n";
  cout << "  Correct determinant =              " << exact << "\n";

  return;
# undef N
}
//****************************************************************************80

void test295 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST295 tests R8GE_DILU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 3
# define NROW 3
# define N NROW * NCOL
# define M N

  double a[M*N];
  double *d;
  int i;
  int j;
  int k;

  cout << "\n";
  cout << "TEST295\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_DILU returns the DILU factors.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  for ( i = 1; i <= NROW * NCOL; i++ )
  {
    for ( j = 1; j <= NROW * NCOL; j++ )
    {
      if ( i == j )
      {
        a[i-1+(j-1)*M] = 4.0;
      }
      else if ( 
        i == j + 1 ||
        i == j - 1 ||
        i == j + NROW ||
        i == j - NROW )
      {
        a[i-1+(j-1)*M] = -1.0;
      }
      else
      {
        a[i-1+(j-1)*M] = 0.0;
      }
    }
  }

  r8ge_print ( M, N, a, "  Matrix A:" );
//
//  Compute the incomplete LU factorization.
//
  d = r8ge_dilu ( M, N, a );

  r8vec_print ( M, d, "  DILU factor of A:" );

  delete [] d;

  return;
# undef M
# undef N
# undef NCOL
# undef NROW
}
//********************************************************************

void test30 ( )

//********************************************************************
//
//  Purpose:
//
//    TEST30 tests R8GE_FA, R8GE_SL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;
//
  cout << "\n";
  cout << "TEST30\n";
  cout << "  R8GE_FA factors a general linear system,\n";
  cout << "  R8GE_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  Random matrix A:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mxv ( N, N, a, x );
//
//  Factor the matrix.
//
  info = r8ge_fa ( N, a, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8GE_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  job = 0;
  delete [] x;
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
//
//  Compute the corresponding right hand side.
//
  job = 0;
  delete [] b;
  b = r8ge_ml ( N, a, pivot, x, job );
//
//  Solve the system
//
  job = 0;
  delete [] x;
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 1;
  delete [] b;
  b = r8ge_ml ( N, a, pivot, x, job );
//
//  Solve the system
//
  job = 1;
  delete [] x;
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution of transposed system:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test31 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST31 tests R8GE_FA, R8GE_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST31\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_FA computes the LU factors,\n";
  cout << "  R8GE_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  The matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mxv ( N, N, a, x );
//
//  Factor the matrix.
//
  info = r8ge_fa ( N, a, pivot );
 
  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST31 - Fatal error!\n";
    cout << "  R8GE_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Display the gory details.
//
  r8mat_print ( N, N, a, "  The compressed LU factors:" );

  i4vec_print ( N, pivot, "  The pivot vector P:" );
//
//  Solve the linear system.
//
  job = 0;
  delete [] x;
  x = r8ge_sl_new ( N, a, pivot, b, job );

  r8vec_print ( N, x, "  Solution:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test315 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST315 tests R8GE_ILU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NCOL 3
# define NROW 3
# define N NROW * NCOL
# define M N

  double a[M*N];
  int i;
  int j;
  int k;
  double l[M*M];
  double lu[M*N];
  double u[M*N];

  cout << "\n";
  cout << "TEST315\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_ILU returns the ILU factors.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  for ( i = 0; i < NROW * NCOL; i++ )
  {
    for ( j = 0; j < NROW * NCOL; j++ )
    {
      if ( i == j )
      {
        a[i+j*M] = 4.0;
      }
      else if ( 
        i == j + 1 | 
        i == j - 1 | 
        i == j + NROW | 
        i == j - NROW 
      )
      {
        a[i+j*M] = -1.0;
      }
      else
      {
        a[i+j*M] = 0.0;
      }
    }
  }

  r8ge_print ( M, N, a, "  Matrix A:" );
//
//  Compute the incomplete LU factorization.
//
  r8ge_ilu ( M, N, a, l, u );

  r8ge_print ( M, M, l, "  Factor L:" );

  r8ge_print ( M, N, u, "  Factor U:" );

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < M; i++ )
    {
      lu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        lu[i+j*M] = lu[i+j*M] + l[i+k*M] * u[k+j*M];
      }
    }
  }

  r8ge_print ( M, N, lu, "  Product L*U:" );

  return;
# undef M
# undef N
# undef NCOL
# undef NROW
}
//****************************************************************************80

void test317 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST317 tests R8GE_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST317\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_INDICATOR sets up an indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  a = r8ge_indicator ( M, N );

  r8ge_print ( M, N, a, "  The R8GE indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test32 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST32 tests R8GE_NP_FA, R8GE_NP_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST32\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_NP_FA computes the LU factors without pivoting,\n";
  cout << "  R8GE_NP_SL solves factored systems.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( N, N, &seed );
//
//  Set the desired solution.
//
  x = new double[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mxv ( N, N, a, x );
//
//  Factor the matrix.
//
  info = r8ge_np_fa ( N, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST32 - Fatal error!\n";
    cout << "  R8GE_NP_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  job = 0;
  delete [] x;
  x = r8ge_np_sl ( N, a, b, job );
 
  r8vec_print_some ( N, x, 10, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 0;
  delete [] b;
  b = r8ge_np_ml ( N, a, x, job );
//
//  Solve the system
//
  job = 0;
  delete [] x;
  x = r8ge_np_sl ( N, a, b, job );

  r8vec_print_some ( N, x, 10, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 1;
  delete [] b;
  b = r8ge_np_ml ( N, a, x, job );
//
//  Solve the system
//
  job = 1;
  delete [] x;
  x = r8ge_np_sl ( N, a, b, job );

  r8vec_print ( N, x, "  Solution of transposed system:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test33 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST33 tests R8GE_NP_FA, R8GE_NP_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *a_lu;
  double *b;
  double *c;
  int i;
  int info;
  int seed = 123456789;
  int j;

  cout << "\n";
  cout << "TEST33\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_NP_FA computes LU factors without pivoting,\n";
  cout << "  R8GE_NP_INVERSE computes the inverse.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( N, N, &seed );

  r8ge_print ( N, N, a, "  The random matrix:" );
//
//  Factor and invert the matrix.
//
  a_lu = new double[N*N];

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_np_fa ( N, a_lu );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST33 - Fatal error!\n";
    cout << "  R8GE_NP_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  b = r8ge_np_inverse ( N, a_lu );

  r8ge_print ( N, N, b, "  The inverse matrix:" );
//
//  Compute A * B = C.
//
  c = r8ge_mxm ( N, a, b );

  r8ge_print ( N, N, c, "  The product:" );

  delete [] a;
  delete [] a_lu;
  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test34 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST34 tests R8GE_FS_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST34\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_FS_NEW factors and solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( N, N, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mxv ( N, N, a, x );
//
//  Factor and solve the system.
//
  delete [] x;
  x = r8ge_fs_new ( N, a, b );
  
  if ( x == NULL )
  {
    cout << "\n";
    cout << "TEST34 - Fatal error!\n";
    cout << "  R8GE_FS_NEW reports the matrix is singular.\n";
    return;
  }
  r8vec_print ( N, x, "  Solution:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test345 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST345 tests R8GE_FSS_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 209
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NB 3

  double *a;
  double *b;
  int i;
  int info;
  int j;
  int k;
  int n = N;
  int nb = NB;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST345\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_FSS_NEW factors and solves multiple linear systems.\n";
  cout << "\n";
  cout << "  Matrix order N = " << n << "\n";
  cout << "  Number of systems NB = " << nb << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( n, n, &seed );
//
//  Set the desired solutions.
//
  b = new double[n * nb];

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = 1.0;
  }
  k = 0;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 );
  }
  k = 1;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( i % 3 ) + 1;
  }
  k = 2;
  for ( i = 0; i < n; i++ )
  {
    b[i+k*n] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i+k*n] = b[i+k*n] + a[i+j*n] * x[j];
    }
  }
//
//  Factor and solve the system.
//
  delete [] x;

  x = r8ge_fss_new ( n, a, nb, b );
  
  r8mat_print ( n, nb, x, "  Solutions:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef N
# undef NB
}
//****************************************************************************80

void test35 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST35 tests R8GE_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  double a_lu[N*N];
  double *a_inverse;
  double *c;
  int i;
  int info;
  int j;
  int pivot[N];
  double x = 2.0;
  double y = 3.0;

  cout << "\n";
  cout << "TEST35\n";
  cout << "  R8GE_INVERSE inverts a general matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      if ( i == j )
      {
        a[i+j*N] = x + y;
      }
      else
      {
        a[i+j*N] = y;
      }
    }
  }

  r8ge_print ( N, N, a, "  Matrix A:" );
//
//  Factor and invert the matrix.
//
  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < N; i++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }

  info = r8ge_fa ( N, a_lu, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST35 - Fatal error!\n";
    cout << "  R8GE_FA reports the matrix is singular.\n";
    return;
  }
  a_inverse = r8ge_inverse ( N, a_lu, pivot );

  r8ge_print ( N, N, a_inverse, "  Inverse matrix B:" );
//
//  Check.
//
  c = r8ge_mxm ( N, a, a_inverse );

  r8ge_print ( N, N, c, "  Product matrix:" );

  delete [] a_inverse;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test36 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST36 tests R8GE_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int pivot[N];
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST36\n";
  cout << "  R8GE_ML computes A*x or A'*X\n";
  cout << "    where A has been factored by R8GE_FA.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r8ge_random ( N, N, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8ge_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ge_vxm ( N, N, a, x );
    }
//
//  Factor the matrix.
//
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "  Fatal error!\n";
      cout << "  R8GE_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      continue;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r8ge_ml ( N, a, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x\n" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x\n" );
    }

    delete [] a;
    delete [] b;
    delete [] b2;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test37 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST37 tests R8GE_NP_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST37\n";
  cout << "  For a matrix in general storage,\n"; 
  cout << "  R8GE_NP_ML computes A*x or A'*X\n";
  cout << "    where A has been factored by R8GE_NP_FA.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    a = r8ge_random ( N, N, &seed );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8ge_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ge_vxm ( N, N, a, x );
    }
//
//  Factor the matrix.
//
    info = r8ge_np_fa ( N, a );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "TEST37 - Fatal error!\n";
      cout << "  R8GE_NP_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      continue;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r8ge_np_ml ( N, a, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] a;
    delete [] b;
    delete [] b2;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test38 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST38 tests R8GE_MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 3

  double *amn;
  double *anm;
  double *bm;
  double *bn;
  double *cm;
  double *cn;
  int info;
  int i;
  int job;
  int pivot[M+N];
  int seed = 123456789;
  char trans;
  double *xm;
  double *xn;

  cout << "\n";
  cout << "TEST38\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_MU computes A*x or A'*X\n";
  cout << "    where A has been factored by R8GE_TRF.\n";
//
//  First test.
//  A is 5 x 3, and we compute A * x.
//
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  job = 0;
  trans = 'N';

  amn = r8ge_random ( M, N, &seed );
  xn = r8vec_indicator_new ( N );
  cm = r8ge_mxv ( M, N, amn, xn );

  info = r8ge_trf ( M, N, amn, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST38 - Fatal error!\n";
    cout << "  R8GE_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  bm = r8ge_mu ( M, N, amn, trans, pivot, xn );

  r8vec2_print_some ( M, cm, bm, 10, "  A*x and PLU*x" );

  delete [] bm;
  delete [] cm;
  delete [] xn;
  delete [] amn;
//
//  Second test.
//  A is 5 x 3, and we compute A' * x.
//
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  job = 1;
  trans = 'T';

  amn = r8ge_random ( M, N, &seed );

  xm = r8vec_indicator_new ( M );

  cn = r8ge_vxm ( M, N, amn, xm );

  info = r8ge_trf ( M, N, amn, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST38 - Fatal error!\n";
    cout << "  R8GE_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  bn = r8ge_mu ( M, N, amn, trans, pivot, xm );

  r8vec2_print_some ( N, cn, bn, 10, "  A'*x and (PLU)'*x" );

  delete [] bn;
  delete [] cn;
  delete [] xm;
  delete [] amn;
//
//  Third test.
//  A is 3 x 5, and we compute A * x.
//
  cout << "\n";
  cout << "  Matrix rows M =    " << N << "\n";
  cout << "  Matrix columns N = " << M << "\n";

  job = 0;
  trans = 'N';

  anm = r8ge_random ( N, M, &seed );

  xm = r8vec_indicator_new ( M );
  cn = r8ge_mxv ( N, M, anm, xm );
 
  info = r8ge_trf ( N, M, anm, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST38 - Fatal error!\n";
    cout << "  R8GE_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  bn = r8ge_mu ( N, M, anm, trans, pivot, xm );

  r8vec2_print_some ( N, cn, bn, 10, "  A*x and PLU*x" );

  delete [] bn;
  delete [] cn;
  delete [] xm;
  delete [] anm;
//
//  Fourth test.
//  A is 3 x 5, and we compute A' * x.
//
  cout << "\n";
  cout << "  Matrix rows M =    " << N << "\n";
  cout << "  Matrix columns N = " << M << "\n";

  job = 1;
  trans = 'T';

  anm = r8ge_random ( N, M, &seed );

  xn = r8vec_indicator_new ( N );
  cm = r8ge_vxm ( N, M, anm, xn );

  info = r8ge_trf ( N, M, anm, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST38 - Fatal error!\n";
    cout << "  R8GE_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  bm = r8ge_mu ( N, M, anm, trans, pivot, xn );

  r8vec2_print_some ( M, cm, bm, 10, "  A'*x and (PLU)'*x" );

  delete [] bm;
  delete [] cm;
  delete [] xn;
  delete [] anm;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test385 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST385 tests R8GE_PLU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 4

  double *a;
  int i;
  int j;
  int k;
  double l[M*M];
  double lu[M*N];
  double p[M*M];
  double plu[M*N];
  int seed;
  double u[M*N];

  cout << "\n";
  cout << "TEST385\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_PLU returns the PLU factors of a matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  seed = 123456789;

  a = r8ge_random ( M, N, &seed );

  r8ge_print ( M, N, a, "  Matrix A:" );
//
//  Compute the PLU factors.
//
  r8ge_plu ( M, N, a, p, l, u );

  r8ge_print ( M, M, p, "  Factor P:" );

  r8ge_print ( M, M, l, "  Factor L:" );

  r8ge_print ( M, N, u, "  Factor U:" );

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      lu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        lu[i+j*M] = lu[i+j*M] + l[i+k*M] * u[k+j*M];
      }
    }
  }

  for ( i = 0; i < M; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      plu[i+j*M] = 0.0;
      for ( k = 0; k < M; k++ )
      {
        plu[i+j*M] = plu[i+j*M] + p[i+k*M] * lu[k+j*M];
      }
    }
  }
  r8ge_print ( M, N, plu, "  Product P*L*U:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test39 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST39 tests R8GE_POLY.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 12

  double a[N*N];
  int i;
  int j;
  double *p;
  double p_true[N+1] = {   
         1.0,    -23.0,    231.0,  -1330.0,    4845.0, 
    -11628.0,  18564.0, -19448.0,  12870.0,   -5005.0, 
      1001.0,    -78.0,      1.0 };

  cout << "\n";
  cout << "TEST39\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_POLY computes the characteristic polynomial.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) i4_min ( i + 1, j + 1 );
    }
  }
//
//  Get the characteristic polynomial.
//
  p = r8ge_poly ( N, a );
//
//  Compare.
//
  r8vec2_print_some ( N+1, p, p_true, 10, "I, P(I), True P(I)" );

  delete [] p;

  return;
# undef N
}
//****************************************************************************80

void test40 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST40 tests R8GE_SL_IT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 6

  double *a;
  double a_lu[N*N];
  double b[N];
  int i;
  int info;
  int j;
  int job;
  int pivot[N];
  double *r;
  double *x;
  double *x_new;

  cout << "\n";
  cout << "TEST40\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_SL_IT applies one step of iterative \n";
  cout << "    refinement to an R8GE_SL solution.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the coefficient matrix.
//
  a = hilbert_inverse ( N );
//
//  Set the right hand side b.
//
  for ( i = 0; i < N-1; i++ )
  {
    b[i] = 0.0;
  }
  b[N-1] = 1.0;
//
//  It is necessary to keep both an unfactored and factored copy
//  of the coefficient matrix.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a_lu[i+j*N] = a[i+j*N];
    }
  }
//
//  Compute the factored coefficient matrix.
//
  info = r8ge_fa ( N, a_lu, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST40 - Fatal error!\n";
    cout << "  R8GE_FA declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the system.
//  (Careful!  R8GE_SL overwrites the right hand side with the solution!)
//
  job = 0;
  x = r8ge_sl_new ( N, a_lu, pivot, b, job );
//
//  Compute and print the residual.
//
  r = r8ge_res ( N, N, a, x, b );

  r8vec2_print_some ( N, x, r, 10, "  i, x, b-A*x" );
//
//  Take a few steps of iterative refinement.
//
  for ( j = 1; j <= 5; j++ )
  {
    cout << "\n";
    cout << "Iterative refinement step " << j << "\n";
    cout << "\n";
//
//  Improve the solution.
//
    job = 0;
    x_new = r8ge_sl_it ( N, a, a_lu, pivot, b, job, x );
//
//  Compute and print the residual.
//
    delete [] r;

    r = r8ge_res ( N, N, a, x_new, b );

    r8vec2_print_some ( N, x_new, r, 10, "  i, x, b-A*x" );

    for ( i = 0; i < N; i++ )
    {
      x[i] = x_new[i];
    }
    delete [] x_new;
  }

  delete [] a;
  delete [] r;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test41 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST41 tests R8GE_TRF, R8GE_TRS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define M 5
# define NRHS 1

  double a[N*N];
  double b[N*NRHS];
  int i;
  int info;
  int j;
  int pivot[N];
  double *x;

  cout << "\n";
  cout << "TEST41\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_TRF computes the LU factors,\n";
  cout << "  R8GE_TRS solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( i == j )
      {
        a[i+j*N] = 2.0;
      }
      else if ( i == j - 1 )
      {
        a[i+j*N] = - 1.0;
      }
      else if ( i == j + 1 )
      {
        a[i+j*N] = - 1.0;
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  info = r8ge_trf ( M, N, a, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST41 - Fatal error!\n";
    cout << "  R8GE_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  for ( i = 0; i < N-1; i++ )
  {
    b[i+0*N] = 0.0;
  }
  b[N-1+0*N] = ( double ) ( N + 1 );

  x = r8ge_trs ( N, NRHS, 'N', a, pivot, b );

  if ( x == NULL )
  {
    cout << "\n";
    cout << "TEST41 - Fatal error!\n";
    cout << "  R8GE_TRS returned an error condition!\n";
    return;
  }
  r8vec_print ( N, x, "  Solution:" );

  for ( i = 0; i < N-1; i++ )
  {
    b[i+0*N] = 0.0;
  }
  b[N-1+0*N] = ( double ) ( N + 1 );

  delete [] x;

  x = r8ge_trs ( N, NRHS, 'T', a, pivot, b );

  if ( x == NULL )
  {
    cout << "\n";
    cout << "TEST41 - Fatal error!\n";
    cout << "  R8GE_TRS returned an error condition!\n";
    return;
  }

  r8vec_print ( N, x, "  Solution to transposed system:" );

  delete [] x;

  return;
# undef M
# undef N
# undef NRHS
}
//****************************************************************************80

void test42 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST42 tests R8GE_NP_TRF, R8GE_NP_TRM, R8GE_NP_TRS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 10
# define N 10
# define NRHS 1

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST42\n";
  cout << "  For a matrix in general storage,\n";
  cout << "  R8GE_NP_TRF factors without pivoting,\n";
  cout << "  R8GE_NP_TRS solves factored systems.\n";
  cout << "  R8GE_NP_TRM computes A*X for factored A.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8ge_random ( M, N, &seed );
//
//  Set the desired solution.
//
  x = new double[N];
  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mxv ( M, N, a, x );
//
//  Factor the matrix.
//
  info = r8ge_np_trf ( M, N, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST42 - Fatal error!\n";
    cout << "  R8GE_NP_TRF declares the matrix is singular!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
//
//  Solve the linear system.
//
  delete [] x;
  x = r8ge_np_trs ( N, NRHS, 'N', a, b );
 
  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST42 - Fatal error!\n";
    cout << "  R8GE_TRS returned an error condition!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }
  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 0;
  delete [] b;
  b = r8ge_np_trm ( M, N, a, x, job );
//
//  Solve the system
//
  x = r8ge_np_trs ( N, NRHS, 'N', a, b );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST42 - Fatal error!\n";
    cout << "  R8GE_TRS returned an error condition!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  job = 1;
  delete [] b;
  b = r8ge_np_trm ( M, N, a, x, job );
//
//  Solve the system.
//
  delete [] x;
  x = r8ge_np_trs ( N, NRHS, 'T', a, b );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST42 - Fatal error!\n";
    cout << "  R8GE_TRS returned an error condition!\n";
    cout << "  The value of INFO is " << info << "\n";
    return;
  }

  r8vec_print ( N, x, "  Solution of transposed system:" );

  delete [] a;
  delete [] b;
  delete [] x;

  return;
# undef M
# undef N
# undef NRHS
}
//****************************************************************************80

void test422 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST422 tests R8CC_GET, R8CC_IJK, R8CC_INC, R8CC_KIJ, R8CC_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int i;
  int j;
  int k;
  int rowind[NZ_NUM] = {
    1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;
  int test;
  int test_num = 20;
  double value;

  cout << "\n";
  cout << "TEST422\n";
  cout << "  For a matrix in the R8CC format,\n";
  cout << "  (double precision Harwell-Boeing Unsymmetric Assembled)\n";
  cout << "  R8CC_GET gets an entry;\n";
  cout << "  R8CC_IJK gets K from (I,J)\n";
  cout << "  R8CC_INC increments an entry;\n";
  cout << "  R8CC_KIJ gets (I,J) from K;\n";
  cout << "  R8CC_SET sets an entry;\n";
  cout << "\n";
  cout << "  Matrix rows M    = " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Nonzeros NZ_NUM  = " << NZ_NUM << "\n";

  i4vec_print ( N+1, colptr, "  The COLPTR vector:" );

  i4vec_print ( NZ_NUM, rowind, "  The ROWIND vector:" );
//
//  Initialize the matrix to random values.
//
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, 
    "  The initial R8CC matrix:" );

  cout << "\n";
  cout << "  R8CC_IJK locates some (I,J) entries.\n";
  cout << "\n";
  cout << "         I         J         K\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i4_uniform ( 1, M, &seed );
    j = i4_uniform ( 1, N, &seed );
    k = r8cc_ijk ( M, N, NZ_NUM, colptr, rowind, i, j );
    cout << "  " << setw(8) << i
         << "  " << setw(8) << j
         << "  " << setw(8) << k << "\n";
  }

  cout << "\n";
  cout << "  R8CC_KIJ locates some K entries.\n";
  cout << "\n";
  cout << "         K         I         J\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    cout << "  " << setw(8) << k
         << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";
  }

  cout << "\n";
  cout << "  R8CC_SET sets 10 entries at random.\n";
  cout << "\n";
  cout << "         I         J         K      NEW_VALUE\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = 100.0 + ( double ) ( test );
    r8cc_set ( M, N, NZ_NUM, colptr, rowind, a, i, j, value );
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << j
         << "  " << setw(8)  << k 
         << "  " << setw(14) << value << "\n";
  }

  cout << "\n";
  cout << "  R8CC_INC increments 10 entries at random.\n";
  cout << "\n";
  cout << "         I         J         K       NEW_VALUE\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = 20.0 + ( double ) ( test );
    r8cc_inc ( M, N, NZ_NUM, colptr, rowind, a, i, j, value );
    value = r8cc_get ( M, N, NZ_NUM, colptr, rowind, a, i, j );
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << j
         << "  " << setw(8)  << k 
         << "  " << setw(14) << value << "\n";
  }

  cout << "\n";
  cout << "  R8CC_GET retrieves 10 entries.\n";
  cout << "\n";
  cout << "         I         J         K      VALUE\n";
  cout << "\n";

  for ( test = 1; test <= 10; test++ )
  {
    k = i4_uniform ( 1, NZ_NUM, &seed );
    r8cc_kij ( M, N, NZ_NUM, colptr, rowind, k, &i, &j );
    value = r8cc_get ( M, N, NZ_NUM, colptr, rowind, a, i, j );
    cout << "  " << setw(8)  << i
         << "  " << setw(8)  << j
         << "  " << setw(8)  << k 
         << "  " << setw(14) << value << "\n";
  }

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The final R8CC matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test423 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST423 tests R8CC_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };

  cout << "\n";
  cout << "TEST423\n";
  cout << "  R8CC_INDICATOR sets up an SHBUA indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M    = " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Nonzeros NZ_NUM  = " << NZ_NUM << "\n";

  a = r8cc_indicator ( M, N, NZ_NUM, colptr, rowind );

  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The R8CC indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test425 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST425 tests R8CC_MXV, R8CC_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  double *bm;
  double *bn;
  double *c;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int i;
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;
  double xn[N];
  double xm[M];

  cout << "\n";
  cout << "TEST425\n";
  cout << "  R8CC_MXV multiplies an SHBUA matrix by a vector;\n";
  cout << "  R8CC_VXM multiplies a vector by an SHBUA matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M    = " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Nonzeros NZ_NUM  = " << NZ_NUM << "\n";
//
//  Set the matrix.
//
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );
//
//  Make an R8GE copy.
//
  c = r8cc_to_r8ge ( M, N, NZ_NUM, colptr, rowind, a );
//
//  Print the R8GE copy.
//
  r8ge_print ( N, N, c, "  The R8CC matrix, in R8GE form:" );
//
//  Compute A * xn = bm.
//
  xn[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    xn[i] = 0.0;
  }
  xn[N-1] = -1.0;

  r8vec_print ( N, xn, "  The vector xn:" );

  bm = r8cc_mxv ( M, N, NZ_NUM, colptr, rowind, a, xn );

  r8vec_print ( M, bm, "  The product A * xn:" );
//
//  Compute xm * A = bn.
//
  xm[0] = 1.0;
  for ( i = 1; i < M-1; i++ )
  {
    xm[i] = 0.0;
  }
  xm[M-1] = -1.0;

  bn = r8cc_vxm ( M, N, NZ_NUM, colptr, rowind, a, xm );

  r8vec_print ( N, bn, "  The product A' * xm:" );

  delete [] a;
  delete [] bm;
  delete [] bn;
  delete [] c;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test426 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST426 tests R8CC_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5
# define NZ_NUM 12

  double *a;
  int colptr[N+1] = { 1, 4, 6, 8, 10, 13 };
  int rowind[NZ_NUM] = { 1, 2, 4, 1, 2, 3, 5, 4, 5, 1, 2, 5 };
  int seed = 123456789;

  cout << "\n";
  cout << "TEST426\n";
  cout << "  R8CC_PRINT prints an R8CC matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M    = " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Nonzeros NZ_NUM  = " << NZ_NUM << "\n";
//
//  Set the matrix.
//
  a = r8cc_random ( M, N, NZ_NUM, colptr, rowind, &seed );
//
//  Print the matrix.
//
  r8cc_print ( M, N, NZ_NUM, colptr, rowind, a, "  The R8CC matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test428 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST428 tests R8LT_INDICATOR;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 6
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST428\n";
  cout << "  For a matrix in lower triangular storage,\n";
  cout << "  R8LT_INDICATOR sets up an indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  a = r8lt_indicator ( M, N );

  r8lt_print ( M, N, a, "  The R8LT indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test43 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST43 tests R8LT_SL, R8LT_MXV, R8LT_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int j;
  int job;
  double *x;

  cout << "\n";
  cout << "TEST43\n";
  cout << "  For a matrix in lower triangular storage,\n";
  cout << "  R8LT_SL solves systems;\n";
  cout << "  R8LT_MXV computes matrix-vector products;\n";
  cout << "  R8LT_VXM computes vector-matrix products;\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      if ( j <= i )
      {
        a[i+j*N] = ( double ) ( j + 1 );
      }
      else
      {
        a[i+j*N] = 0.0;
      }
    }
  }

  r8lt_print ( N, N, a, "  The lower triangular matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8lt_mxv ( N, N, a, x );
    }
    else
    {
      b = r8lt_vxm ( N, N, a, x );
    }
//
//  Solve the linear system.
//
    delete [] x;
    x = r8lt_sl ( N, a, b, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to the transposed system:" );
    }
    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test44 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST44 tests R8LT_DET, R8LT_INVERSE, R8LT_MXM, R8LT_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  double *c;
  double det;
  int i;
  int j;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST44\n";
  cout << "  For a matrix in lower triangular storage,\n";
  cout << "  R8LT_DET computes the determinant.\n";
  cout << "  R8LT_INVERSE computes the inverse.\n";
  cout << "  R8LT_MXM computes matrix products.\n";
  cout << "  R8LT_RANDOM sets a random value.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8lt_random ( N, N, &seed );

  r8lt_print ( N, N, a, "  Matrix A:" );
//
//  Compute the determinant.
//
  det = r8lt_det ( N, a );

  cout << "\n";
  cout << "  Determinant is " << det << "\n";
//
//  Compute the inverse matrix.
//
  b = r8lt_inverse ( N, a );

  r8lt_print ( N, N, b, "  Inverse matrix B:" );
//
//  Check.
//
  c = r8lt_mxm ( N, a, b );

  r8lt_print ( N, N, c, "  Product A * B:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test443 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST443 tests R8NCF_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//   John Burkardt
//
{
# define M 7
# define N 5
# define NZ_NUM 15

  double *a;
  int rowcol[2*NZ_NUM] = {
    1, 1, 
    2, 2,
    3, 3,
    4, 4, 
    5, 5, 
    2, 1, 
    5, 1, 
    1, 2, 
    5, 2,
    1, 4, 
    2, 4, 
    3, 4,
    4, 5, 
    4, 6, 
    1, 7 };

  cout << "\n";
  cout << "TEST443\n";
  cout << "  R8NCF_INDICATOR sets up a R8NCF indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";
  cout << "  Matrix nonzeros =  " << NZ_NUM << "\n";

  a = r8ncf_indicator ( M, N, NZ_NUM, rowcol );

  r8ncf_print ( M, N, NZ_NUM, rowcol, a, "  The R8NCF indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test445 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST445 tests R8PBL_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9
# define MU 3

  double *a;
  int i;
  int j;

  cout << "\n";
  cout << "TEST445\n";
  cout << "  R8PBL_INDICATOR sets up a R8PBL indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Bandwidth MU = " << MU << "\n";

  a = r8pbl_indicator ( N, MU );

  r8pbl_print ( N, MU, a, "  The R8PBL indicator matrix:" );

  delete [] a;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test45 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST45 tests R8PBU_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 50
# define MU 1

  double a[(MU+1)*N];
  double *b;
  int i;
  int j;
  double *r;
  double res_max;
  double *x;
  double *x_init;

  cout << "\n";
  cout << "TEST45\n";
  cout << "  R8PBU_CG applies the conjugate gradient method\n";
  cout << "    to a symmetric positive definite banded\n";
  cout << "    linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix values.
//
  a[0+0*(MU+1)] = 0.0;
  for ( j = 1; j < N; j++ )
  {
    a[0+j*(MU+1)] = -1.0;
  }
  for ( j = 0; j < N; j++ )
  {
    a[1+j*(MU+1)] = 2.0;
  }

  r8pbu_print_some ( N, MU, a, 1, 1, 10, 10, "The symmetric banded matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the right hand side.
//
  b = r8pbu_mxv ( N, MU, a, x );
//
//  Set the approximate solution.
//
  x_init = new double[N];
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = 1.0;
  }
//
//  Call the conjugate gradient method.
//
  delete [] x;

  x = r8pbu_cg ( N, MU, a, b, x_init );

  r8vec_print_some ( N, x, 10, "  Solution:" );
//
//  Compute the residual, A*x-b
//
  r = r8pbu_mxv ( N, MU, a, x );

  res_max = 0.0;
  for ( i = 0; i < N; i++ )
  { 
    res_max = r8_max ( res_max, r8_abs ( r[i] - b[i] ) );
  }

  cout << "\n";
  cout << "  Maximum residual = " << res_max << "\n";

  delete [] r;
  delete [] x;
  delete [] x_init;
 
  return;
# undef MU
# undef N
}
//****************************************************************************80

void test46 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST46 tests R8PBU_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define MU 3

  double *a;
  double *a2;
  double *a_lu;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  cout << " \n";  
  cout << "TEST46\n";
  cout << "  R8PBU_DET, determinant of a positive definite\n";
  cout << "    symmetric banded matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8pbu_random ( N, MU, &seed );

  r8pbu_print ( N, MU, a, "  The R8PBU matrix:" );
//
//  Copy the matrix into a general array.
//
  a2 = r8pbu_to_r8ge ( N, MU, a );
//
//  Factor the matrix.
//
  a_lu = r8pbu_fa ( N, MU, a );
  r8pbu_print ( N, MU, a_lu, "  The factored R8PBU matrix:" );
//
//  Compute the determinant.
//
  det = r8pbu_det ( N, MU, a_lu );

  cout << "\n";
  cout << "  R8PBU_DET computes the determinant = " << det << "\n";
//
//  Factor the general matrix.
//
  info = r8ge_fa ( N, a2, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, a2, pivot );

  cout << "  R8GE_DET computes the determinant =  " << det << "\n";

  delete [] a;
  delete [] a2;
  delete [] a_lu;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test47 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST47 tests R8PBU_FA, R8PBU_RANDOM, R8PBU_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 50
# define MU 1

  double *a;
  double *a_lu;
  double *b;
  int i;
  int info;
  int j;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST47\n";
  cout << "  For a banded positive definite symmetric matrix,\n";
  cout << "  R8PBU_FA computes the LU factors.\n";
  cout << "  R8PBU_RANDOM sets a random value.\n";
  cout << "  R8PBU_SL solves a linear system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
  cout << "\n";
//
//  Set the matrix values.
//
  a = r8pbu_random ( N, MU, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the right hand side.
//
  b = r8pbu_mxv ( N, MU, a, x );
//
//  Factor the matrix.
//
  a_lu = r8pbu_fa ( N, MU, a );
//
//  Solve the linear system.
//
  delete [] x;
  x = r8pbu_sl ( N, MU, a_lu, b );

  r8vec_print_some ( N, x, 10, "  Solution:" );

  delete a;
  delete a_lu;
  delete b; 
  delete x;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test48 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST48 tests R8PBU_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define MU 3

  double *a;
  double *a_lu;
  double *b;
  double *b2;
  int i;
  int info;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST48\n";
  cout << "  R8PBU_ML computes A*x \n";
  cout << "    where A has been factored by R8PBU_FA.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";
//
//  Set the matrix.
//
  a = r8pbu_random ( N, MU, &seed );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8pbu_mxv ( N, MU, a, x );
//
//  Factor the matrix.
//
  a_lu = r8pbu_fa ( N, MU, a );
//
//  Now multiply the factored matrix times solution to get right hand side again.
//
  b2 = r8pbu_ml ( N, MU, a_lu, x );

  r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );

  delete [] a;
  delete [] a_lu;
  delete [] b;
  delete [] b2;
  delete [] x;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test485 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST485 tests R8PBU_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9
# define MU 3

  double *a;

  cout << " \n";
  cout << "TEST485\n";
  cout << "  R8PBU_INDICATOR sets up an SPBU indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Bandwidth MU = " << MU << "\n";

  a = r8pbu_indicator ( N, MU );

  r8pbu_print ( N, MU, a, "  The R8PBU indicator matrix:" );

  delete [] a;

  return;
# undef MU
# undef N
}
//****************************************************************************80

void test49 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST49 tests R8PBU_SOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 50
# define MU 1

  double a[(MU+1)*N];
  double *b;
  double *b2;
  double eps;
  double err;
  int i;
  int itchk;
  int itmax;
  int j;
  int k;
  double omega;
  double pi = 3.14159265;
  double t;
  double *x;
  double x_init[N];

  cout << "\n";
  cout << "TEST49\n";
  cout << "  R8PBU_SOR, SOR routine for iterative\n";
  cout << "    solution of A*x=b.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Upper bandwidth MU = " << MU << "\n";

  for ( k = 1; k <= 3; k++ )
  {
    if ( k == 1 )
    {
      omega = 0.25;
    }
    else if ( k == 2 )
    {
      omega = 0.75;
    }
    else
    {
      omega = 1.00;
    }
//
//  Set matrix values.
//
    a[0+0*(MU+1)] = 0.0;
    for ( j = 1; j < N; j++ )
    {
      a[0+j*(MU+1)] = -1.0;
    }

    for ( j = 0; j < N; j++ )
    {
      a[1+j*(MU+1)] = 2.0;
    }
//
//  Set the desired solution.
//
    x = new double[N];
    for ( i = 0; i < N; i++ )
    {
      t = pi * ( double ) ( i ) / ( double ) ( N - 1 );
      x[i] = sin ( t );
    }
//
//  Compute the right hand side.
//
    b = r8pbu_mxv ( N, MU, a, x );
//
//  Set the initial solution estimate.
//
    for ( i = 0; i < N; i++ )
    {
      x_init[i] = 1.0;
    }
 
    itchk = 1;
    itmax = 8000;
    eps = 0.0001;

    delete [] x;
    x = r8pbu_sor ( N, MU, a, b, eps, itchk, itmax, omega, x_init );
//
//  Compute residual, A*x-b
//
    b2 = r8pbu_mxv ( N, MU, a, x );
 
    err = 0.0;
    for ( i = 0; i < N; i++ )
    {
      err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
    }
 
    cout << "\n";
    cout << "SOR iteration.\n";
    cout << "\n";
    cout << "  Relaxation factor OMEGA = " << omega << "\n";

    r8vec_print_some ( N, x, 10, "  Solution:" );

    cout << "\n";
    cout << "  Maximum error = " << err << "\n";

    delete [] b;
    delete [] b2;
    delete [] x;
  }
 
  return;
# undef MU
# undef N
}
//****************************************************************************80

void test50 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST50 tests R8PO_FA, R8PO_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int info;
  int j;
  double *r;
  double *x;

  cout << "\n";
  cout << "TEST50\n";
  cout << "  R8PO_FA factors a positive definite symmetric\n";
  cout << "    linear system,\n";
  cout << "  R8PO_SL solves a factored system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
  b = r8po_mxv ( N, a, x );
//
//  Factor the matrix.
//
  r = r8po_fa ( N, a );

  if ( r == NULL )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8PO_FA declares the matrix is singular!\n";
    return;
  }
//
//  Solve the linear system.
//
  delete [] x;
  x = r8po_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );
//
//  Set the desired solution.
//
  delete [] x;
  x = new double[N];

  for ( i = 0; i < N; i++ )
  {
    x[i] = 1.0;
  }
//
//  Compute the corresponding right hand side, using the factored matrix.
//
  delete [] b;
  b = r8po_ml ( N, r, x );
//
//  Solve the linear system.
//
  delete [] x;
  x = r8po_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );

  delete [] b;
  delete [] r;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test505 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST505 tests R8PO_FA;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N*N];
  int i;
  int j;
  int k;
  double *r;

  cout << "\n";
  cout << "TEST505\n";
  cout << "  R8PO_FA factors a positive definite symmetric\n";
  cout << "    linear system,\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }

  r8po_print ( N, a, "  The matrix A:" );
//
//  Factor the matrix.
//
  r = r8po_fa ( N, a );

  if ( r == NULL )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8PO_FA declares the matrix is singular!\n";
    return;
  }
  r8ut_print ( N, N, r, "  The factor R (an R8UT matrix):" );
//
//  Compute the product R' * R.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = 0.0;
      for ( k = 0; k < N; k++ )
      {
        a[i+j*N] = a[i+j*N] + r[k+i*N] * r[k+j*N];
      }
    }
  }

  r8ge_print ( N, N, a, "  The product R' * R:" );

  delete [] r;

  return;
# undef N
}
//****************************************************************************80

void test51 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST51 tests R8PO_DET, R8PO_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N];
  double *a_inv;
  double *c;
  double det;
  int i;
  int info;
  int j;
  int k;
  double *r;

  cout << "\n";
  cout << "TEST51\n";
  cout << "  For a symmetric positive definite matrix\n";
  cout << "    factored by R8PO_FA,\n";
  cout << "  R8PO_DET computes the determinant;\n";
  cout << "  R8PO_INVERSE computes the inverse.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  for ( i = 0; i < N; i++ )
  {
    for ( j = 0; j < N; j++ )
    {
      a[i+j*N] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }

  r8po_print ( N, a, "  Matrix A:" );
//
//  Factor the matrix.
//
  r = r8po_fa ( N, a );
//
//  Compute the determinant.
//
  det = r8po_det ( N, r );
 
  cout << "\n";
  cout << "  Matrix determinant = " << det << "\n";
//
//  Compute the inverse.
//
  a_inv = r8po_inverse ( N, r );

  r8po_print ( N, a_inv, "  Inverse matrix A_INV:" );
//
//  Check.
//
  c = r8po_mxm ( N, a, a_inv );

  r8po_print ( N, c, "  Product A * A_INV:" );

  delete [] a_inv;
  delete [] c;
  delete [] r;

  return;
# undef N
}
//****************************************************************************80

void test515 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST515 tests R8PO_FA, R8PO_SL, R8PO_ML, R8PO_MXV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int info;
  int seed = 123456789;
  double *r;
  double *x;

  cout << "\n";
  cout << "TEST515\n";
  cout << "  For a positive definite symmetric matrix,\n";
  cout << "  R8PO_FA computes the Cholesky factor.\n";
  cout << "  R8PO_SL solves a factored linear system.\n";
  cout << "  R8PO_MXV multiplies unfactored A * x\n";
  cout << "  R8PO_ML multiplies factored A * x\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8po_random ( N, &seed );

  r8po_print ( N, a, "  The matrix A:" );
//
//  Compute the desired right hand side.
//
  x = r8vec_indicator_new ( N );

  b = r8po_mxv ( N, a, x );

  r8vec_print ( N, b, "  Right hand side, computed by R8PO_MXV" );
//
//  Factor the matrix.
//
  r = r8po_fa ( N, a );
//
//  Solve the linear system.
//
  delete [] x;
  x = r8po_sl ( N, r, b );

  r8vec_print ( N, x, "  Solution (should be 1,2,3...)" );
//
//  Recompute the desired right hand side.
//  Since A has been factored, we have to call R8PO_ML now.
//
  delete [] x;
  x = r8vec_indicator_new ( N );

  delete [] b;
  b = r8po_ml ( N, a, x );

  r8vec_print ( N, b, "  Right hand side, computed by R8PO_ML" );
//
//  Solve the linear system.
//
  delete [] x;
  x = r8po_sl ( N, a, b );

  r8vec_print ( N, x, "  Solution (should be 1,2,3...)" );

  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x;

  return;
}
//****************************************************************************80

void test517 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST517 tests R8PO_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST517\n";
  cout << "  R8PO_INDICATOR sets up an R8PO indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8po_indicator ( N );

  r8po_print ( N, a, "  The R8PO indicator matrix:" );
 
  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test52 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST52 tests R8PO_RANDOM, R8PO_TO_R8GE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST52\n";
  cout << "  R8PO_RANDOM computes a random positive definite\n";
  cout << "    symmetric matrix.\n";
  cout << "  R8PO_TO_R8GE converts an R8PO matrix to R8GE format.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8po_random ( N, &seed );

  r8po_print ( N, a, "  The random R8PO matrix:" );
 
  r8ge_print ( N, N, a, "  The random R8PO matrix (printed by R8GE_PRINT):" );

  b = r8po_to_r8ge ( N, a );

  r8ge_print ( N, N, b, "  The random R8GE matrix (printed by R8GE_PRINT):" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test525 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST525 tests R8PP_FA, R8PP_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int j;
  double *r;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST525\n";
  cout << "  R8PP_FA factors an R8PP system,\n";
  cout << "  R8PP_SL solves an R8PP system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8pp_random ( N, &seed );

  r8pp_print ( N, a, "  The R8PP matrix:" );
//
//  Set the desired solution.
//
  x = r8vec_indicator_new ( N );

  r8vec_print ( N, x, "  The desired solution:" );
//
//  Compute the corresponding right hand side.
//
  b = r8pp_mxv ( N, a, x );

  r8vec_print ( N, b, "  The right hand side:" );
//
//  Factor the matrix.
//
  r = r8pp_fa ( N, a );

  if ( r == NULL )
  {
    cout << "\n";
    cout << "  Fatal error!\n";
    cout << "  R8PP_FA declares the matrix is singular!\n";
    return;
  }
  cout << "\n";
  cout << "  The R8PP matrix has been factored.\n";
//
//  Solve the linear system.
//
  delete [] x;
  x = r8pp_sl ( N, r, b );
 
  r8vec_print ( N, x, "  Solution:" );

  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test527 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST527 tests R8PP_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST527\n";
  cout << "  R8PP_INDICATOR sets up an R8PP indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8pp_indicator ( N );

  r8pp_print ( N, a, "  The R8PP indicator matrix:" );
 
  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test53 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST53 tests R8PP_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST53\n";
  cout << "  R8PP_RANDOM, compute a random positive definite\n";
  cout << "    symmetric packed matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8pp_random ( N, &seed );

  r8pp_print ( N, a, "  The matrix (printed by R8PP_PRINT):" );
 
  b = r8pp_to_r8ge ( N, a );

  r8ge_print ( N, N, b, "  The random R8PP matrix (printed by R8GE_PRINT):" );

  delete [] a;
  delete [] b;

  return;
# undef N
}
//****************************************************************************80

void test534 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST534 tests R8S3_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define NZ_NUM ( 3 * N - 2 )

  double a[NZ_NUM];
  int col[NZ_NUM];
  int i;
  int isym;
  int j;
  int k;
  string output_file = "r8s3_matrix.txt";
  int row[NZ_NUM];

  cout << "\n";
  cout << "TEST534\n";
  cout << "  For a R8S3 matrix,\n";
  cout << "  R8S3_WRITE writes the matrix to a file.\n";
  cout << "\n";
  cout << "  Matrix order N =         " << N << "\n";
  cout << "  Matrix nonzeros NZ_NUM = " << NZ_NUM << "\n";

  isym = 0;
//
//  Set the matrix values.
//
  k = 0;
  for ( i = 1; i <= N; i++ )
  {

    j = i;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 2; i <= N; i++ )
  {
    j = i - 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 1; i <= N-1; i++ )
  {
    j = i + 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  r8s3_print_some ( N, N, NZ_NUM, isym, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8S3 matrix:" );

  r8s3_write ( N, NZ_NUM, isym, row, col, a, output_file );

  cout << "  R8S3_WRITE wrote the matrix data to \"" << output_file << "\".\n";

  return;
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test535 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST535 tests R8S3_READ, R8S3_READ_SIZE..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int *col;
  int isym;
  int n;
  int nz_num;
  string input_file = "r8s3_matrix.txt";
  int *row;

  cout << "\n";
  cout << "TEST535\n";
  cout << "  For a R8S3 matrix,\n";
  cout << "  R8S3_READ reads a matrix from a file.\n";
  cout << "  R8S3_READ_SIZE reads the sizes of the matrix from a file.\n";

  r8s3_read_size ( input_file, &n, &nz_num );

  cout << "\n";
  cout << "  R8S3_READ_SIZE reports matrix size data:\n";
  cout << "\n";

  cout << "\n";
  cout << "  Matrix order N =         " << n << "\n";
  cout << "  Matrix nonzeros NZ_NUM = " << nz_num << "\n";

  row = new int[nz_num];
  col = new int[nz_num];
  a = new double[nz_num];

  r8s3_read ( input_file, n, nz_num, row, col, a );

  isym = 0;

  r8s3_print_some ( n, n, nz_num, isym, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8S3 matrix:" );

  cout << "\n";
  cout << "  Deleting the matrix data file \"" << input_file << "\"\n";

  file_delete ( input_file );

  delete [] a;
  delete [] col;
  delete [] row;

  return;
}
//****************************************************************************80

void test54 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST54 tests R8SD_CG.
//
//  Discussion:
//
//    NX and NY are the number of grid points in the X and Y directions.
//    N is the number of unknowns.
//    NDIAG is the number of nonzero diagonals we will store.  We only
//    store the main diagonal, and the superdiagonals.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define NDIAG 3
# define NX 10
# define NY 10
# define N ( NX * NY )

  double a[N*NDIAG];
  double *b;
  double *b2;
  double err;
  int i;
  int j;
  int k;
  int offset[NDIAG] = { 0, 1, NX };
  double *x;
  double x_init[N];

  cout << "\n";
  cout << "TEST54\n";
  cout << "  R8SD_CG applies the conjugate gradient method\n";
  cout << "    to a symmetric positive definite linear\n";
  cout << "    system stored by diagonals.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Matrix diagonals NDIAG = " << NDIAG << "\n";
  cout << "\n";
//
//  OFFSET tells us where the nonzero diagonals are.  It does this
//  by recording how "high" or to the right the diagonals are from
//  the main diagonal.
//
//
//  Now we compute the numbers that go into the diagonals.  For this
//  problem, we could simply store a column of 4's, and two columns of
//  -1's, but I wanted to go through the motions of thinking about the
//  value of each entry.  "K" counts the row of the original matrix
//  that we are working on.
//
  k = 0;

  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
//
//  Central
//
      a[k+0*N] = 4.0;
//
//  East ( = West )
//
      if ( i == NX )
      {
        a[k+1*N] = 0.0;
      }
      else
      {
        a[k+1*N] = -1.0;
      }
//
//  North ( = South )
//
      if ( j == NY )
      {
        a[k+2*N] = 0.0;
      }
      else
      {
        a[k+2*N] = -1.0;
      }
      k = k + 1;
    }
  }
//
//  Print some of the matrix.
//
  r8sd_print_some ( N, NDIAG, offset, a, 1, 1, 10, 10, 
    "  First 10 rows and columns of matrix." );
//
//  Set the desired solution.
//
  x = new double[N];
  k = 0;
  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
      x[k] = ( double ) ( 10 * i + j );
      k = k + 1;
    }
  }
//
//  Compute the corresponding right hand side.
//
  b = r8sd_mxv ( N, NDIAG, offset, a, x );

  r8vec_print_some ( N, b, 10, "  Right hand side:" );
//
//  Set X to zero so no one accuses us of cheating.
//
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = 0.0;
  }
//
//  Call the conjugate gradient method.
//
  delete [] x;
  x = r8sd_cg ( N, NDIAG, offset, a, b, x_init );
//
//  Compute the residual, A*x-b
//
  b2 = r8sd_mxv ( N, NDIAG, offset, a, x );
 
  err = 0.0;
  for ( i = 0; i < N; i++ )
  {
    err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
  }
 
  r8vec_print_some ( N, x, 10, "  Solution:" );

  cout << "\n";
  cout << "  Maximum residual = " << err << "\n";
//
//  Note that if we're not satisfied with the solution, we can
//  call again, using the computed X as our starting estimate.
//
//  Call the conjugate gradient method AGAIN.
//
  for ( i = 0; i < N; i++ )
  {
    x_init[i] = x[i];
  }
  delete [] x;

  x = r8sd_cg ( N, NDIAG, offset, a, b, x_init );
//
//  Compute the residual, A*x-b
//
  delete [] b2;
  b2 = r8sd_mxv ( N, NDIAG, offset, a, x );
 
  err = 0.0;
  for ( i = 0; i < N; i++ )
  {
    err = r8_max ( err, r8_abs ( b2[i] - b[i] ) );
  }
 
  r8vec_print_some ( N, x, 10, "  Second attempt at solution:" );

  cout << "\n";
  cout << "  Maximum residual of second attempt = " << err << "\n";

  delete [] b;
  delete [] b2;
  delete [] x;

  return;
# undef N
# undef NDIAG
# undef NX
# undef NY
}
//****************************************************************************80

void test55 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST55 tests R8SD_CG.
//
//  Discussion:
//
//    This is a sample demonstration of how to compute some eigenvalues
//    and corresponding eigenvectors of a matrix.  The matrix is the
//    discretized Laplacian operator, which can be stored by diagonals,
//    and handled by the conjugate gradient method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define MAXVEC 3
# define NDIAG 3
# define NX 10
# define NY 10
# define N ( NX * NY )
# define PI 3.141592653589

  double a[N*NDIAG];
  double del;
  double dot;
  double eval;
  int i;
  int iter;
  int i4vec;
  int j;
  int k;
  double lambda;
  double lambda_old;
  double lamvec[MAXVEC];
  double norm;
  int nvec;
  int offset[NDIAG] = { 0, 1, NX };
  double vec[N*MAXVEC];
  double x[N];
  double *xnew;

  cout << "\n";
  cout << "TEST55\n";
  cout << "  R8SD_CG is used for linear equation solving\n";
  cout << "    in a demonstration of inverse iteration to\n";
  cout << "    compute eigenvalues and eigenvectors of a\n";
  cout << "    symmetric matrix stored by diagonals.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Matrix diagonals NDIAG = " << NDIAG << "\n";
  cout << "\n";
  cout << "  Here are 25 of the smallest eigenvalues:\n";
  cout << "\n";
  cout << "  I, J, eigenvalue(I,J):\n";
  cout << "\n";

  for ( i = 1; i <= i4_min ( 5, NX ); i++ )
  {
    for ( j = 1; j <= i4_min ( 5, NY ); j++ )
    {
      eval = 4.0 
        - 2.0 * cos ( ( double ) ( i ) * PI / ( double ) ( NX + 1 ) ) 
        - 2.0 * cos ( ( double ) ( j ) * PI / ( double ) ( NY + 1 ) );
      cout << setw(6) <<  i    << "  "
           << setw(6) <<  j    << "  "
           << setw(12) << eval << "\n";
    }
  }
//
//  Now we compute the numbers that go into the diagonals.  For this
//  problem, we could simply store a column of 4's, and two columns of
//  -1's, but I wanted to go through the motions of thinking about the
//  value of each entry.  "K" counts the row of the original matrix
//  that we are working on.
//
  k = 0;
  for ( j = 1; j <= NY; j++ )
  {
    for ( i = 1; i <= NX; i++ )
    {
//
//  Central
//
      a[k+0*N] = 4.0;
//
//  East ( = West )
//
      if ( i == NX )
      {
        a[k+1*N] = 0.0;
      }
      else
      {
        a[k+1*N] = -1.0;
      }
//
//  North ( = South )
//
      if ( j == NY )
      {
        a[k+2*N] = 0.0;
      }
      else
      {
        a[k+2*N] = -1.0;
      }
      k = k + 1;
    }
  }

  nvec = 0;
//
//  Set the starting eigenvector and eigenvalue estimates.
//
  for ( ; ; )
  {
    cout << "\n";

    lambda = 0.0;

    k = 0;
    for ( j = 1; j <= NY; j++ )
    {
      for ( i = 1; i <= NX; i++ )
      {
        x[k] = 1.0;
        k = k + 1;
      }
    }
//
//  Remove any components of previous eigenvectors.
//
    for ( i4vec = 0; i4vec < nvec; i4vec++ )
    {
      dot = 0.0;
      for ( i = 0; i < N; i++ )
      {
        dot = dot + x[i] * vec[i+i4vec*N];
      }
      for ( i = 0; i < N; i++ )
      {
        x[i] = x[i] - dot * vec[i+i4vec*N];
      }
    }

    xnew = new double[N];
    for ( i = 0; i < N; i++ )
    {
      xnew[i] = x[i];
    }
//
//  Iterate
//
    for ( iter = 1; iter <= 40; iter++ )
    {
      norm = 0.0;
      for ( i = 0; i < N; i++ )
      {
        norm = norm + xnew[i] * xnew[i];
      }
      norm = sqrt ( norm );

      for ( i = 0; i < N; i++ )
      {
        xnew[i] = xnew[i] / norm;
      }

      lambda_old = lambda;
      lambda = 1.0 / norm;
//
//  Check for convergence.
//
      if ( 1 < iter )
      {
        del = r8_abs ( lambda - lambda_old );

        if ( del < 0.000001 )
        {
          cout << "\n";
          cout << "Lambda estimate = " << lambda << "\n";
          cout << "Converged on step " << iter << "\n";
          break;
        }
      }
//
//  Call the conjugate gradient method, solving A * XNEW = X.
//
      for ( i = 0; i < N; i++ )
      {
        x[i] = xnew[i];
      }
      delete xnew;

      xnew = r8sd_cg ( N, NDIAG, offset, a, x, x );

      for ( i4vec = 0; i4vec < nvec; i4vec++ )
      {
        dot = 0.0;
        for ( i = 0; i < N; i++ )
        {
          dot = dot + xnew[i] * vec[i+i4vec*N];
        }
        for ( i = 0; i < N; i++ )
        {
          xnew[i] = xnew[i] - dot * vec[i+i4vec*N];
        }
      }

    }

    lamvec[nvec] = lambda;
    for ( i = 0; i < N; i++ )
    {
      vec[i+nvec*N] = xnew[i];
    }
    nvec = nvec + 1;

    if ( MAXVEC <= nvec )
    {
      break;
    }

  }

  delete xnew;

  return;
# undef MAXVEC
# undef N
# undef NDIAG
# undef NX
# undef NY
# undef PI
}
//****************************************************************************80

void test555 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST555 tests R8SD_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define NDIAG 3

  double *a;
  int offset[NDIAG] = { 0, 1, 3 };

  cout << "\n";
  cout << "TEST555\n";
  cout << "  R8SD_INDICATOR sets up an R8SD indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
  cout << "  Matrix diagonals NDIAG = " << NDIAG << "\n";

  a = r8sd_indicator ( N, NDIAG, offset );

  r8sd_print ( N, NDIAG, offset, a, "  The R8SD indicator matrix:" );

  delete [] a;

  return;
# undef N
# undef NDIAG
}
//****************************************************************************80

void test56 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST56 tests R8SM_ML.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 7

  double a[M*N];
  double *b;
  double *b2;
  int info;
  int i;
  int job;
  int pivot[N];
  int seed = 123456789;
  double u[M];
  double v[N];
  double *x;

  cout << "\n";
  cout << "TEST56\n";
  cout << "  R8SM_ML computes A*x or A'*X\n";
  cout << "    where A is a Sherman Morrison matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the matrix.
//
    r8sm_random ( M, N, a, u, v, &seed );

    r8sm_print ( M, N, a, u, v, "  The Sherman Morrison matrix:" );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8sm_mxv ( M, N, a, u, v, x );
    }
    else
    {
      b = r8sm_vxm ( M, N, a, u, v, x );
    }
//
//  Factor the matrix.
//
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "  Fatal error!\n";
      cout << "  R8GE_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      return;
    }
//
//  Now multiply factored matrix times solution to get right hand side again.
//
    b2 = r8sm_ml ( N, a, u, v, pivot, x, job );

    if ( job == 0 )
    {
      r8vec2_print_some ( N, b, b2, 10, "  A*x and PLU*x" );
    }
    else
    {
      r8vec2_print_some ( N, b, b2, 10, "  A'*x and (PLU)'*x" );
    }

    delete [] b;
    delete [] b2;
    delete [] x;
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test57 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST57 tests R8SM_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 5

  double a[M*N];
  double *b;
  int i;
  int ierror;
  int info;
  int job;
  int pivot[N];
  int seed = 123456789;
  double u[M];
  double v[N];
  double *x;

  cout << "\n";
  cout << "TEST57\n";
  cout << "  R8SM_SL implements the Sherman-Morrison method \n";
  cout << "    for solving a perturbed linear system.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  for ( job = 1; 0 <= job; job-- )
  {
//
//  Set the matrix.
//
    r8sm_random ( M, N, a, u, v, &seed );

    r8sm_print ( M, N, a, u, v, "  The Sherman-Morrison matrix A:" );
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8sm_mxv ( M, N, a, u, v, x );
    }
    else
    {
      b = r8sm_vxm ( M, N, a, u, v, x );
    }

    r8vec_print ( N, b, "  The right hand side vector B:" );
//
//  Factor the matrix.
//
    info = r8ge_fa ( N, a, pivot );

    if ( info != 0 )
    {
      cout << "\n";
      cout << "  Fatal error!\n";
      cout << "  R8GE_FA declares the matrix is singular!\n";
      cout << "  The value of INFO is " << info << "\n";
      continue;
    }
//
//  Solve the linear system.
//
    b = r8sm_sl ( N, a, u, v, b, pivot, job );
 
    if ( job == 0 )
    {
      r8vec_print ( N, b, "  Solution to A * X = B:" );
    }
    else
    {
      r8vec_print ( N, b, "  Solution to A' * X = B:" );
    }

    delete [] b;
  }

  return;
# undef M
# undef N
}
//****************************************************************************80

void test5705 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST5705 tests R8SP_IJ_TO_K.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NZ_NUM 10

  bool check;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int i;
  int j;
  int k;
  int m = 7;
  int n = 5;
  int nz_num = NZ_NUM;
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };

  cout << "\n";
  cout << "TEST5705\n";
  cout << "  R8SP_IJ_TO_K returns the R8SP index of (I,J).\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << m << "\n";
  cout << "  Matrix columns N = " << n << "\n";
  cout << "  Matrix nonzeros =  " << nz_num << "\n";

  check = r8sp_check ( m, n, nz_num, row, col );

  if ( !check )
  {
    cout << "\n";
    cout << "R8SP_CHECK - Error!\n";
    cout << "  The matrix is not in the proper sorted format.\n";
    return;
  }

  cout << "\n";
  cout << "         I         J         K\n";
  cout << "\n";
  for ( i = 1; i <= m; i++ )
  {
    for ( j = 1; j <= n; j++ )
    {
      k = r8sp_ij_to_k ( nz_num, row, col, i, j );

      cout << "  " << setw(8) << i
           << "  " << setw(8) << j
           << "  " << setw(8) << k << "\n";
    }
  }

  return;
# undef NZ_NUM
}
//****************************************************************************80

void test571 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST571 tests R8SP_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 5
# define NZ_NUM 10

  double *a;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };

  cout << "\n";
  cout << "TEST571\n";
  cout << "  R8SP_INDICATOR sets up an R8SP indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M =          " << M << "\n";
  cout << "  Matrix columns N =       " << N << "\n";
  cout << "  Matrix nonzeros NZ_NUM = " << NZ_NUM << "\n";

  a = r8sp_indicator ( M, N, NZ_NUM, row, col );

  r8sp_print ( M, N, NZ_NUM, row, col, a, "  The R8SP indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test572 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST572 tests R8SP_MXV, R8SP_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 7
# define N 5
# define NZ_NUM 10

  double *a;
  double *b;
  double *c;
  int col[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1 };
  int i;
  int row[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7 };
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST572\n";
  cout << "  R8SP_MXV multiplies an R8SP matrix by a vector;\n";
  cout << "  R8SP_VXM multiplies a vector by an R8SP matrix;\n";
  cout << "\n";
  cout << "  Matrix rows M =          " << M << "\n";
  cout << "  Matrix columns N =       " << N << "\n";
  cout << "  Matrix nonzeros NZ_NUM = " << NZ_NUM << "\n";
//
//  Set the matrix.
//
  a = r8sp_random ( M, N, NZ_NUM, row, col, &seed );
//
//  Make an R8GE copy.
//
  c = r8sp_to_r8ge ( M, N, NZ_NUM, row, col, a );
//
//  Print the R8GE copy.
//
  r8ge_print ( M, N, c, "  The R8SP matrix, in R8GE form:" );

  x = new double[N];

  x[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    x[i] = 0.0;
  }
  x[N-1] = -1.0;

  r8vec_print ( N, x, "  The vector x:" );

  b = r8sp_mxv ( M, N, NZ_NUM, row, col, a, x );

  r8vec_print ( M, b, "  The product A * x:" );

  delete [] x;
  x = new double[M];

  x[0] = 1.0;
  for ( i = 1; i < M-1; i++ )
  {
    x[i] = 0.0;
  }
  x[M-1] = -1.0;

  r8vec_print ( M, x, "  The vector x:" );

  delete [] b;

  b = r8sp_vxm ( M, N, NZ_NUM, row, col, a, x );

  r8vec_print ( N, b, "  The product A' * x:" );

  delete [] a;
  delete [] b;
  delete [] c;
  delete [] x;

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test5722 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST5722 tests R8SP_PRINT.
//
//  Discussion:
//
//    Because MATLAB seems to allow a R8SP matrix to store the same index
//    several times, presumably with the matrix entry being the SUM of
//    these occurrences, I modified R8SP_PRINT to handle this situation
//    (I hope).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 5
# define N 7
# define NZ_NUM 12

  double a[NZ_NUM] = {
    21.0,  51.0, 12.0, 52.0, 14.0, 
    24.0,  34.0, 45.0, 46.0, 17.0, 
   100.0, 200.0 };
  int col[NZ_NUM] = { 1, 1, 2, 2, 4, 4, 4, 5, 6, 7, 2, 4 };
  int row[NZ_NUM] = { 2, 5, 1, 5, 1, 2, 3, 4, 4, 1, 1, 3 };


  cout << "\n";
  cout << "TEST5722\n";
  cout << "  R8SP_PRINT prints a R8SP matrix;\n";
  cout << "  In this example, we have listed several matrix\n";
  cout << "  locations TWICE.  R8SP_PRINT should compute the\n";
  cout << "  sum of these values.\n";
  cout << "\n";
  cout << "  In particular, we want A(1,2) = 112 and A(3,4) = 234;\n";
  cout << "\n";
  cout << "  Matrix rows M =          " << M << "\n";
  cout << "  Matrix columns N =       " << N << "\n";
  cout << "  Matrix nonzeros NZ_NUM = " << NZ_NUM << "\n";

  r8sp_print ( M, N, NZ_NUM, row, col, a, "  The R8SP matrix:" );

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test5724 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST5724 tests R8SP_WRITE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 100
# define N 100
# define NZ_NUM ( 3 * N - 2 )

  double a[NZ_NUM];
  int col[NZ_NUM];
  int i;
  int j;
  int k;
  string output_file = "r8sp_matrix.txt";
  int row[NZ_NUM];

  cout << "\n";
  cout << "TEST5724\n";
  cout << "  For a R8SP matrix,\n";
  cout << "  R8SP_WRITE writes the matrix to a file.\n";
  cout << "\n";
  cout << "  Matrix number of rows M =    " << M << "\n";
  cout << "  Matrix number of columns N = " << N << "\n";
  cout << "  Matrix nonzeros NZ_NUM =     " << NZ_NUM << "\n";
//
//  Set the matrix values.
//
  k = 0;
  for ( i = 1; i <= N; i++ )
  {

    j = i;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 2; i <= N; i++ )
  {
    j = i - 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  for ( i = 1; i <= N-1; i++ )
  {
    j = i + 1;
    row[k] = i;
    col[k] = j;
    a[k] = ( double ) ( 100 * i + j );
    k = k + 1;
  }

  r8sp_print_some ( M, N, NZ_NUM, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8SP matrix:" );

  r8sp_write ( M, N, NZ_NUM, row, col, a, output_file );

  cout << "  R8SP_WRITE wrote the matrix data to \"" << output_file << "\".\n";

  return;
# undef M
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void test5725 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST5725 tests R8SP_READ, R8SP_READ_SIZE..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int *col;
  int m;
  int n;
  int nz_num;
  string input_file = "r8sp_matrix.txt";
  int *row;

  cout << "\n";
  cout << "TEST5725\n";
  cout << "  For a R8SP matrix,\n";
  cout << "  R8SP_READ reads a matrix from a file.\n";
  cout << "  R8SP_READ_SIZE reads the sizes of the matrix from a file.\n";

  r8sp_read_size ( input_file, &m, &n, &nz_num );

  cout << "\n";
  cout << "  R8SP_READ_SIZE reports matrix size data:\n";
  cout << "\n";

  cout << "\n";
  cout << "  Matrix number of rows M =    " << m << "\n";
  cout << "  Matrix number of columns N = " << n << "\n";
  cout << "  Matrix nonzeros NZ_NUM =     " << nz_num << "\n";

  row = new int[nz_num];
  col = new int[nz_num];
  a = new double[nz_num];

  r8sp_read ( input_file, m, n, nz_num, row, col, a );

  r8sp_print_some ( m, n, nz_num, row, col, a, 1, 1, 
    10, 10, "  Initial 10x10 block of R8SP matrix:" );

  cout << "\n";
  cout << "  Deleting the matrix data file \"" << input_file << "\"\n";

  file_delete ( input_file );

  delete [] a;
  delete [] col;
  delete [] row;

  return;
}
//****************************************************************************80

void test573 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST573 tests R8SR_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define NZ 7

  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };

  cout << "\n";
  cout << "TEST573\n";
  cout << "  R8SR_INDICATOR sets up an R8SR indicator matrix;\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  r8sr_indicator ( N, NZ, row, col, diag, off );

  r8sr_print ( N, NZ, row, col, diag, off, "  The R8SR indicator matrix:" );

  return;
# undef N
# undef NZ
}
//****************************************************************************80

void test574 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST574 tests R8SR_MXV, R8SR_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define NZ 7

  double *b;
  double *c;
  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  int i;
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;
  double x[N];

  cout << "\n";
  cout << "TEST574\n";
  cout << "  R8SR_MXV multiplies an R8SR matrix by a vector;\n";
  cout << "  R8SR_VXM multiplies a vector by an R8SR matrix;\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
//
//  Make an R8GE copy.
//
  c = r8sr_to_r8ge ( N, NZ, row, col, diag, off );
//
//  Print the R8GE copy.
//
  r8ge_print ( N, N, c, "  The R8SR matrix, in R8GE form:" );

  x[0] = 1.0;
  for ( i = 1; i < N-1; i++ )
  {
    x[i] = 0.0;
  }
  x[N-1] = -1.0;

  r8vec_print ( N, x, "  The vector x:" );

  b = r8sr_mxv ( N, NZ, row, col, diag, off, x );

  r8vec_print ( N, b, "  The product A * x:" );

  delete [] b;
  b = r8sr_vxm ( N, NZ, row, col, diag, off, x );

  r8vec_print ( N, b, "  The product A' * x:" );

  delete [] b;
  delete [] c;

  return;
# undef N
# undef NZ
}
//****************************************************************************80

void test5745 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST5745 tests R8SR_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define NZ 7

  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;

  cout << "\n";
  cout << "TEST5745\n";
  cout << "  R8SR_PRINT prints an R8SR matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
//
//  Print the matrix.
//
  r8sr_print ( N, NZ, row, col, diag, off, "  The R8SR matrix:" );

  return;
# undef N
# undef NZ
}
//****************************************************************************80

void test575 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST575 tests R8SR_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define NZ 7

  double *b;
  int col[NZ] = { 2, 5, 5, 1, 1, 2, 3 };
  double diag[N];
  double off[NZ];
  int row[N+1] = { 1, 3, 4, 5, 6, 8 };
  int seed = 123456789;

  cout << "\n";
  cout << "TEST575\n";
  cout << "  R8SR_RANDOM randomizes an R8SR matrix\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  r8sr_random ( N, NZ, row, col, diag, off, &seed );
//
//  Make an R8GE copy.
//
  b = r8sr_to_r8ge ( N, NZ, row, col, diag, off );
//
//  Print the R8GE copy.
//
  r8ge_print ( N, N, b, "  The R8SR matrix, in R8GE form:" );

  delete [] b;

  return;
# undef N
# undef NZ
}
//****************************************************************************80

void test577 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST577 tests R8SS_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  double *a;
  int diag[N];
  int na;

  cout << "\n";
  cout << "TEST577\n";
  cout << "  For a symmetric skyline storage matrix,\n";
  cout << "  R8SS_INDICATOR sets up an indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8ss_indicator ( N, &na, diag );

  r8ss_print ( N, na, diag, a, "  The R8SS indicator matrix:" );

  delete [] a;

  return;  
# undef N
}
//****************************************************************************80

void test58 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST58 tests R8SS_MXV, R8SS_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9

  double a[(N*(N+1))/2];
  double *a2;
  double *b;
  double *b2;
  int diag[N];
  int i;
  int ij;
  int ilo;
  int j;
  int na;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST58\n";
  cout << "  For a symmetric skyline storage matrix,\n";
  cout << "  R8SS_MXV computes A*x,\n";
  cout << "  R8SS_PRINT prints it.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  r8ss_random ( N, &na, diag, a, &seed );

  cout << "\n";
  cout << "  Number of nonzero entries stored is " << na << "\n";
  cout << "\n";
  cout << "  Diagonal storage indices:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    cout << setw(6) << i + 1   << "  "
         << setw(6) << diag[i] << "\n";
  }
//
//  Replace the random entries by marker values.
//
  ij = 0;
  for ( j = 1; j <= N; j++ )
  {
    if ( j == 1 )
    {
      ilo = 1;
    }
    else
    {
      ilo = diag[j-2] - diag[j-1] + j + 1;
    }

    for ( i = ilo; i <= j; i++ )
    {
      a[ij] = ( double ) ( 10 * i + j );
      ij = ij + 1;
    }
  }
  r8ss_print ( N, na, diag, a, "  The R8SS matrix:" );
//
//  Copy the matrix into a general matrix.
//
  a2 = r8ss_to_r8ge ( N, na, diag, a );
//
//  Set the vector X.
//
  x = r8vec_indicator_new ( N );
//
//  Compute the product.
//
  b = r8ss_mxv ( N, na, diag, a, x );
//
//  Compute the product using the general matrix.
//
  b2 = r8ge_mxv ( N, N, a2, x );
//
//  Compare the results.
//
  r8vec2_print_some ( N, b, b2, 10, "  R8SS_MXV verse R8GE_MXV" );

  delete [] a2;
  delete [] b;
  delete [] b2;
  delete [] x;

  return;  
# undef N
}
//****************************************************************************80

void test581 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST581 tests R8STO_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double *a;

  cout << "\n";
  cout << "TEST581\n";
  cout << "  R8STO_INDICATOR sets up an R8STO indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8sto_indicator ( N );

  r8sto_print ( N, a, "  The R8STO indicator matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test583 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST583 tests R8STO_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N] = { 4.0, 2.0, 0.8 };
  double *a2;
  double *b;
  double *c;

  cout << "\n";
  cout << "TEST583\n";
  cout << "  R8STO_INVERSE computes the inverse of a positive\n";
  cout << "  definite symmetric Toeplitz matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix A:" );

  b = r8sto_inverse ( N, a );

  r8ge_print ( N, N, b, "  The inverse matrix B:" );

  a2 = r8sto_to_r8ge ( N, a );

  c = r8ge_mxm ( N, a2, b );

  r8ge_print ( N, N, c, "  The product C = A * B:" );

  delete [] a2;
  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test585 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST585 tests R8STO_MXV, R8STO_YW_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N];
  double *b;
  int i;
  int job;
  double r[N+1] = { 1.0, 0.5, 0.2, 0.1 };
  double *x;

  for ( i = 0; i < N; i++ )
  {
    a[i] = r[i];
  }

  cout << "\n";
  cout << " TEST585\n";
  cout << "  R8STO_YW_SL solves the Yule-Walker equations for a\n";
  cout << "  symmetric Toeplitz system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix:" );

  b = new double[N];
  for ( i = 0; i < N; i++ )
  {
    b[i] = -r[i+1];
  }

  r8vec_print ( N, b, "  The right hand side, B:" );

  for ( i = 0; i < N; i++ )
  {
    b[i] = -b[i];
  }

  x = r8sto_yw_sl ( N, b );

  r8vec_print ( N, x, "  The computed solution, X:" );

  delete [] b;
  b = r8sto_mxv ( N, a, x );

  r8vec_print ( N, b, "  The product A*x:" );

  delete [] b;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test587 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST587 tests R8STO_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N] = { 1.0, 0.5, 0.2 };
  double *ax;
  double b[N] = { 4.0, -1.0, 3.0 };
  double *x;

  cout << "\n";
  cout << "TEST587\n";
  cout << "  R8STO_SL solves a positive definite symmetric Toeplitz system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  r8sto_print ( N, a, "  The symmetric Toeplitz matrix A:" );

  r8vec_print ( N, b, "  The right hand side vector b:" );

  x = r8sto_sl ( N, a, b );

  r8vec_print ( N, x, "  The computed solution x:" );

  ax = r8sto_mxv ( N, a, x );

  r8vec_print ( N, ax, "  The product vector A * x:" );

  delete [] ax;
  delete [] x;

  return;
# undef N
}
//****************************************************************************80

void test589 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST589 tests R8TO_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;

  cout << "\n";
  cout << "TEST589\n";
  cout << "  R8TO_INDICATOR sets up an R8TO indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8to_indicator ( N );

  r8to_print ( N, a, "  The R8TO indicator matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test59 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST59 tests R8TO_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *b;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST59\n";
  cout << "  R8TO_SL solves a Toeplitz system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8to_random ( N, &seed );

  r8to_print ( N, a, "  The Toeplitz matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8to_mxv ( N, a, x );
    }
    else
    {
      b = r8to_vxm ( N, a, x );
    }
//
//  Solve the linear system.
//
    delete [] x;
    x = r8to_sl ( N, a, b, job );

    if ( job == 0 )
    {
      r8vec_print_some ( N, x, 10, "  Solution:" );
    }
    else
    {
      r8vec_print_some ( N, x, 10, "  Solution to transposed system:" );
    }
    delete [] b;
    delete [] x;
  }
  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void test60 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST60 tests R8UT_INVERSE, R8UT_DET, R8UT_MXM, R8UT_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  double *c;
  double det;
  int i;
  int j;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST60\n";
  cout << "  For an upper triangular matrix,\n";
  cout << "  R8UT_INVERSE computes the inverse.\n";
  cout << "  R8UT_DET computes the determinant.\n";
  cout << "  R8UT_MXM computes matrix products.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";

  a = r8ut_random ( N, N, &seed );

  r8ut_print ( N, N, a, "  The matrix A:" );
//
//  Compute the determinant.
//
  det = r8ut_det ( N, a );

  cout << "\n";
  cout << "  Determinant is " << det << "\n";
//
//  Compute the inverse matrix B.
//
  b = r8ut_inverse ( N, a );

  r8ut_print ( N, N, b, "  The inverse matrix B:" );
//
//  Check
//
  c = r8ut_mxm ( N, a, b );

  r8ut_print ( N, N, c, "  The product A * B:" );

  delete [] a;
  delete [] b;
  delete [] c;

  return;
# undef N
}
//****************************************************************************80

void test605 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST605 tests R8UT_INDICATOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define M 8
# define N 5

  double *a;

  cout << "\n";
  cout << "TEST605\n";
  cout << "  For an upper triangular matrix,\n";
  cout << "  R8UT_INDICATOR sets up an indicator matrix.\n";
  cout << "\n";
  cout << "  Matrix rows M =    " << M << "\n";
  cout << "  Matrix columns N = " << N << "\n";

  a = r8ut_indicator ( M, N );

  r8ut_print ( M, N, a, "  The R8UT indicator matrix:" );

  delete [] a;

  return;
# undef M
# undef N
}
//****************************************************************************80

void test61 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST61 tests R8UT_SL, R8UT_MXV, R8UT_VXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double a[N*N];
  double *b;
  int i;
  int j;
  int job;
  double *x;

  cout << "\n";
  cout << "TEST61\n";
  cout << "  For an upper triangular matrix,\n";
  cout << "  R8UT_SL solves systems;\n";
  cout << "  R8UT_MXV computes matrix-vector products.\n";
  cout << "  R8UT_VXM computes vector-matrix products.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << " \n";

  for ( i = 1; i <= N; i++ )
  {
    for ( j = 1; j <= N; j++ )
    {
      if ( i <= j )
      {
        a[i-1+(j-1)*N] = ( double ) j;
      }
      else
      {
        a[i-1+(j-1)*N] = 0.0;
      }
    }
  }

  r8ut_print ( N, N, a, "  The upper triangular matrix:" );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8ut_mxv ( N, N, a, x );
    }
    else
    {
      b = r8ut_vxm ( N, N, a, x );
    }
//
//  Solve the linear system.
//
    delete [] x;
    x = r8ut_sl ( N, a, b, job );

    if ( job == 0 )
    {
      r8vec_print ( N, x, "  Solution:" );
    }
    else
    {
      r8vec_print ( N, x, "  Solution to transposed system:" );
    }
    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void test62 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST62 tests R8VM_DET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  double *a;
  double *a2;
  double det;
  int info;
  int pivot[N];
  int seed = 123456789;

  cout << "\n";
  cout << "TEST62\n";
  cout << "  R8VM_DET, determinant of a Vandermonde matrix.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8vm_random ( N, N, &seed );

  r8vm_print ( N, N, a, "  The Vandermonde matrix:" );
//
//  Copy the matrix into a general array.
//
  a2 = r8vm_to_r8ge ( N, N, a );
//
//  Compute the determinant.
//
  det = r8vm_det ( N, a );

  cout << "\n";
  cout << "  R8VM_DET computes the determinant = " << det << "\n";
//
//  Factor the general matrix.
//
  info = r8ge_fa ( N, a2, pivot );
//
//  Compute the determinant.
//
  det = r8ge_det ( N, a2, pivot );

  cout << "  R8GE_DET computes the determinant = " << det << "\n";

  delete [] a;
  delete [] a2;

  return;
# undef N
}
//****************************************************************************80

void test63 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST63 tests R8VM_SL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *b;
  int i;
  int info;
  int job;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST63\n";
  cout << "  R8VM_SL solves a Vandermonde system.\n";
  cout << "\n";
  cout << "  Matrix order N = " << N << "\n";
//
//  Set the matrix.
//
  a = r8vm_random ( N, N, &seed );

  for ( job = 0; job <= 1; job++ )
  {
//
//  Set the desired solution.
//
    x = r8vec_indicator_new ( N );
//
//  Compute the corresponding right hand side.
//
    if ( job == 0 )
    {
      b = r8vm_mxv ( N, N, a, x );
    }
    else
    {
      b = r8vm_vxm ( N, N, a, x );
    }
//
//  Solve the linear system.
//
    delete [] x;
    x = r8vm_sl ( N, a, b, job, &info );

    if ( job == 0 )
    {
      r8vec_print_some ( N, x, 10, "  Solution:" );
    }
    else
    {
      r8vec_print_some ( N, x, 10, "  Solution to transposed system:" );
    }

    delete [] b;
    delete [] x;
  }

  delete [] a;

  return;
# undef N
}
