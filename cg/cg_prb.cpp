# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "cg.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test023 ( );
void test025 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CG_PRB.
//
//  Discussion:
//
//    CG_PRB tests the CG library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CG_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CG library.\n";

  test01 ( );
  test02 ( );
  test023 ( );
  test025 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << " \n";
  cout << "CG_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << " \n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests R8GE_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int n;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test R83_CG, applying CG to a full storage matrix.\n";

  n = 10;
//
//  Let A be the -1 2 -1 matrix.
//
  seed = 123456789;
  a = pds_random ( n, seed );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r8ge_mv ( n, n, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r8ge_cg ( n, a, b, x2 );
//
//  Compute the residual.
//
  r = r8ge_resid ( n, n, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R83_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int n;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test R83_CG, applying CG to an R83 matrix.\n";

  n = 10;
//
//  Let A be the -1 2 -1 matrix.
//
  a = r83_dif2 ( n, n );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r83_mv ( n, n, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r83_cg ( n, a, b, x2 );
//
//  Compute the residual.
//
  r = r83_resid ( n, n, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests R83S_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int n;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  Test R83S_CG, applying CG to an R83S matrix.\n";

  n = 10;
//
//  Let A be the -1 2 -1 matrix.
//
  a = r83s_dif2 ( n, n );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r83s_mv ( n, n, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r83s_cg ( n, a, b, x2 );
//
//  Compute the residual.
//
  r = r83s_resid ( n, n, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests R83T_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int n;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  Test R83T_CG, applying CG to an R83T matrix.\n";

  n = 10;
//
//  Let A be the -1 2 -1 matrix.
//
  a = r83t_dif2 ( n, n );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r83t_mv ( n, n, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r83t_cg ( n, a, b, x2 );
//
//  Compute the residual.
//
  r = r83t_resid ( n, n, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R8PBU_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int mu;
  int n;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test R8PBU_CG, applying CG to an R8PBU matrix.\n";

  n = 10;
  mu = 1;
//
//  Let A be the -1 2 -1 matrix.
//
  a = r8pbu_dif2 ( n, n, mu );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r8pbu_mv ( n, n, mu, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r8pbu_cg ( n, mu, a, b, x2 );
//
//  Compute the residual.
//
  r = r8pbu_resid ( n, n, mu, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests R8SD_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e_norm;
  int i;
  int n;
  int ndiag;
  int *offset;
  double *r;
  double r_norm;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test R8SD_CG, applying CG to an R8SD matrix.\n";

  n = 10;

  ndiag = 2;
  offset = new int[ndiag];
  offset[0] = 0;
  offset[1] = 1;
//
//  Let A be the -1 2 -1 matrix.
//
  a = r8sd_dif2 ( n, n, ndiag, offset );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r8sd_mv ( n, n, ndiag, offset, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r8sd_cg ( n, ndiag, offset, a, b, x2 );
//
//  Compute the residual.
//
  r = r8sd_resid ( n, n, ndiag, offset, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  delete [] offset;
  delete [] r;
  delete [] x1;
  delete [] x2;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests R8SP_CG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int *col;
  double e_norm;
  int i;
  int n;
  int nz_num;
  double *r;
  double r_norm;
  int *row;
  int seed;
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Test R8SP_CG, applying CG to an R8SP matrix.\n";

  n = 10;

  nz_num = 3 * n;
  row = ( int * ) malloc ( nz_num * sizeof ( int ) );
  col = ( int * ) malloc ( nz_num * sizeof ( int ) );
//
//  Let A be the -1 2 -1 matrix.
//
  a = r8sp_dif2 ( n, n, nz_num, row, col );
//
//  Choose a random solution.
//
  seed = 123456789;
  x1 = r8vec_uniform_01_new ( n, seed );
//
//  Compute the corresponding right hand side.
//
  b = r8sp_mv ( n, n, nz_num, row, col, a, x1 );
//
//  Call the CG routine.
//
  x2 = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x2[i] = 1.0;
  }
  r8sp_cg ( n, nz_num, row, col, a, b, x2 );
//
//  Compute the residual.
//
  r = r8sp_resid ( n, n, nz_num, row, col, a, x2, b );
  r_norm = r8vec_norm ( n, r );
//
//  Compute the error.
//
  e_norm = r8vec_diff_norm ( n, x1, x2 );
//
//  Report.
//
  cout << "\n";
  cout << "  Number of variables N = " << n << "\n";
  cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";
  cout << "  Norm of error ||x1-x2|| = " << e_norm << "\n";
//
//  Free memory.
//
  delete [] a;
  delete [] b;
  free ( col );
  delete [] r;
  free ( row );
  delete [] x1;
  delete [] x2;

  return;
}
