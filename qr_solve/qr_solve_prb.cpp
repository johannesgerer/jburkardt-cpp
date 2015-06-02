# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "qr_solve.hpp"
# include "test_ls.hpp"
# include "r8lib.hpp"

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
//    MAIN is the main program for QR_SOLVE_PRB.
//
//  Discussion:
//
//    QR_SOLVE_PRB tests the QR_SOLVE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "QR_SOLVE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the QR_SOLVE library.\n";
  cout << "  QR_SOLVE needs the R8LIB library.\n";
  cout << "  This test also needs the TEST_LS library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QR_SOLVE_PRB\n";
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
//    TEST01 tests NORMAL_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double b_norm;
  int flag;
  int i;
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  NORMAL_SOLVE is a function with a simple interface which\n";
  cout << "  solves a linear system A*x = b in the least squares sense.\n";
  cout << "  Compare a tabulated solution X1 to the NORMAL_SOLVE result X2.\n";
  cout << "\n";
  cout << "  NORMAL_SOLVE cannot be applied when N < M,\n";
  cout << "  or if the matrix does not have full column rank.\n";

  prob_num = p00_prob_num ( );

  cout << "\n";
  cout << "  Number of problems = " << prob_num << "\n";
  cout << "\n";
  cout << "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Get problem size.
//
    m = p00_m ( prob );
    n = p00_n ( prob );
//
//  Retrieve problem data.
//
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
//
//  Use NORMAL_SOLVE on the problem.
//
    x2 = normal_solve ( m, n, a, b, flag );

    if ( flag != 0 )
    {
      cout << "  " << setw(5) << prob
           << "  " << setw(4) << m
           << "  " << setw(4) << n
           << "  " << setw(12) << b_norm
           << "  " << "------------"
           << "  " << setw(12) << x1_norm
           << "  " << "------------"
           << "  " << setw(12) << r1_norm
           << "  " << "------------" << "\n";
    }
    else
    {
      x2_norm = r8vec_norm ( n, x2 );
      r2 = r8mat_mv_new ( m, n, a, x2 );
      for ( i = 0; i < m; i++ )
      {
        r2[i] = r2[i] - b[i];
      }
      r2_norm = r8vec_norm ( m, r2 );
//
//  Compare tabulated and computed solutions.
//
      x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
//
//  Report results for this problem.
//
      cout << "  " << setw(5) << prob
           << "  " << setw(4) << m
           << "  " << setw(4) << n
           << "  " << setw(12) << b_norm
           << "  " << setw(12) << x_diff_norm
           << "  " << setw(12) << x1_norm
           << "  " << setw(12) << x2_norm
           << "  " << setw(12) << r1_norm
           << "  " << setw(12) << r2_norm << "\n";

      delete [] r2;
      delete [] x2;
    }
//
//  Deallocate memory.
//
    delete [] a;
    delete [] b;
    delete [] r1;
    delete [] x1;
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests QR_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
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
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  QR_SOLVE is a function with a simple interface which\n";
  cout << "  solves a linear system A*x = b in the least squares sense.\n";
  cout << "  Compare a tabulated solution X1 to the QR_SOLVE result X2.\n";

  prob_num = p00_prob_num ( );

  cout << "\n";
  cout << "  Number of problems = " << prob_num << "\n";
  cout << "\n";
  cout << "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Get problem size.
//
    m = p00_m ( prob );
    n = p00_n ( prob );
//
//  Retrieve problem data.
//
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
//
//  Use QR_SOLVE on the problem.
//
    x2 = qr_solve ( m, n, a, b );

    x2_norm = r8vec_norm ( n, x2 );
    r2 = r8mat_mv_new ( m, n, a, x2 );
    for ( i = 0; i < m; i++ )
    {
      r2[i] = r2[i] - b[i];
    }
    r2_norm = r8vec_norm ( m, r2 );
//
//  Compare tabulated and computed solutions.
//
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
//
//  Report results for this problem.
//
    cout << "  " << setw(5) << prob
         << "  " << setw(4) << m
         << "  " << setw(4) << n
         << "  " << setw(12) << b_norm
         << "  " << setw(12) << x_diff_norm
         << "  " << setw(12) << x1_norm
         << "  " << setw(12) << x2_norm
         << "  " << setw(12) << r1_norm
         << "  " << setw(12) << r2_norm << "\n";
//
//  Deallocate memory.
//
    delete [] a;
    delete [] b;
    delete [] r1;
    delete [] r2;
    delete [] x1;
    delete [] x2;
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests SVD_SOLVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
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
  int m;
  int n;
  int prob;
  int prob_num;
  double *r1;
  double r1_norm;
  double *r2;
  double r2_norm;
  double x_diff_norm;
  double *x1;
  double x1_norm;
  double *x2;
  double x2_norm;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  SVD_SOLVE is a function with a simple interface which\n";
  cout << "  solves a linear system A*x = b in the least squares sense.\n";
  cout << "  Compare a tabulated solution X1 to the SVD_SOLVE result X2.\n";

  prob_num = p00_prob_num ( );

  cout << "\n";
  cout << "  Number of problems = " << prob_num << "\n";
  cout << "\n";
  cout << "  Index     M     N     ||B||         ||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||\n";
  cout << "\n";

  for ( prob = 1; prob <= prob_num; prob++ )
  {
//
//  Get problem size.
//
    m = p00_m ( prob );
    n = p00_n ( prob );
//
//  Retrieve problem data.
//
    a = p00_a ( prob, m, n );
    b = p00_b ( prob, m );
    x1 = p00_x ( prob, n );

    b_norm = r8vec_norm ( m, b );
    x1_norm = r8vec_norm ( n, x1 );
    r1 = r8mat_mv_new ( m, n, a, x1 );
    for ( i = 0; i < m; i++ )
    {
      r1[i] = r1[i] - b[i];
    }
    r1_norm = r8vec_norm ( m, r1 );
//
//  Use SVD_SOLVE on the problem.
//
    x2 = svd_solve ( m, n, a, b );

    x2_norm = r8vec_norm ( n, x2 );
    r2 = r8mat_mv_new ( m, n, a, x2 );
    for ( i = 0; i < m; i++ )
    {
      r2[i] = r2[i] - b[i];
    }
    r2_norm = r8vec_norm ( m, r2 );
//
//  Compare tabulated and computed solutions.
//
    x_diff_norm = r8vec_norm_affine ( n, x1, x2 );
//
//  Report results for this problem.
//
    cout << "  " << setw(5) << prob
         << "  " << setw(4) << m
         << "  " << setw(4) << n
         << "  " << setw(12) << b_norm
         << "  " << setw(12) << x_diff_norm
         << "  " << setw(12) << x1_norm
         << "  " << setw(12) << x2_norm
         << "  " << setw(12) << r1_norm
         << "  " << setw(12) << r2_norm << "\n";
//
//  Deallocate memory.
//
    delete [] a;
    delete [] b;
    delete [] r1;
    delete [] r2;
    delete [] x1;
    delete [] x2;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DQRLS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double b[5] = { 1.0, 2.3, 4.6, 3.1, 1.2 };
  int i;
  int ind;
  int itask;
  int j;
  int *jpvt;
  int kr;
  int m = 5;
  int n = 3;
  double *qraux;
  double tol;
  double *x;

  a = new double[m*n];
  jpvt = new int[n];
  qraux = new double[n];
  x = new double[n];
//
//  Set up least-squares problem
//  quadratic model, equally-spaced points
//
  cout << "\n";
  cout << "TEST04\n";
  cout << "  DQRLS solves a linear system A*x = b in the least squares sense.\n";

  for ( i = 0; i < m; i++ )
  {
    a[i+0*m] = 1.0;
    for ( j = 1; j < n; j++ )
    {
      a[i+j*m] = a[i+(j-1)*m] * ( double ) ( i + 1 );
    }
  }

  tol = 1.0E-06;

  r8mat_print ( m, n, a, "  Coefficient matrix A:" );

  r8vec_print ( m, b, "  Right hand side b:" );
//
//  Solve least-squares problem
//
  itask = 1;
  ind = dqrls ( a, m, m, n, tol, kr, b, x, b, jpvt, qraux, itask );
//
//  Print results
//
  cout << "\n";
  cout << "  Error code = " << ind << "\n";
  cout << "  Estimated matrix rank = " << kr << "\n";

  r8vec_print ( n, x, "  Least squares solution x:" );

  r8vec_print ( m, b, "  Residuals A*x-b" );

  delete [] a;
  delete [] jpvt;
  delete [] qraux;
  delete [] x;

  return;
}
