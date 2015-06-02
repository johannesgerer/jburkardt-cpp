# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "wishart.hpp"
# include "pdflib.hpp"
# include "rnglib.hpp"

int main ( );
void wishart_test01 ( );
void wishart_test02 ( );
void wishart_test03 ( );
void wishart_test04 ( );
void wishart_test05 ( );
void wishart_test06 ( );
void wishart_test07 ( );
void wishart_test08 ( );
void wishart_test09 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for WISHART_PRB.
//
//  Discussion:
//
//    WISHART_PRB tests the WISHART library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "WISHART_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the WISHART library.\n";

  wishart_test01 ( );
  wishart_test02 ( );
  wishart_test03 ( );
  wishart_test04 ( );
  wishart_test05 ( );
  wishart_test06 ( );
  wishart_test07 ( );
  wishart_test08 ( );
  wishart_test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "WISHART_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void wishart_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST01 demonstrates the unit Wishart sampling function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
  int rot_num;
  double *v;
  double *w;
//
//  Initialize the RNGLIB library.
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST01:\n";
  cout << "  We can compute sample unit Wishart matrices by:\n";
  cout << "    W = wishart_unit_sample ( n, df );\n";
//
//  Set the parameters and call.
//
  n = 5;
  df = 8;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 8 ):" );
  delete [] w;
//
//  Calling again yields a new matrix.
//
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 8 ):" );
  delete [] w;
//
//  Reduce DF
//
  n = 5;
  df = 5;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 5, 5 ):" );
  delete [] w;
//
//  Try a smaller matrix.
//
  n = 3;
  df = 5;
  w = wishart_unit_sample ( n, df );
  r8mat_print ( n, n, w, "  wishart_unit_sample ( 3, 5 ):" );
//
//  What is the eigendecomposition of the matrix?
//
  it_max = 50;
  v = new double[n*n];
  lambda = new double[n];

  jacobi_eigenvalue ( n, w, it_max, v, lambda, it_num, rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
//
//  Free memory.
//
  delete [] lambda;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST02 demonstrates the unit Bartlett sampling function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
  int rot_num;
  double *t;
  double *v;
  double *w;
//
//   Initialize the RNGLIB library.
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST02:\n";
  cout << "  We can compute sample unit Bartlett matrices by:\n";
  cout << "    T = bartlett_unit_sample ( n, df );\n";
//
//   Set the parameters and call.
//
  n = 5;
  df = 8;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 8 ):" );
  delete [] t;
//
//   Calling again yields a new matrix.
//
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 8 ):" );
  delete [] t;
//
//   Reduce DF.
//
  n = 5;
  df = 5;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 5, 5 ):" );
  delete [] t;
//
//   Try a smaller matrix.
//
  n = 3;
  df = 5;
  t = bartlett_unit_sample ( n, df );
  r8mat_print ( n, n, t, "  bartlett_unit_sample ( 3, 5 ):" );
//
//   What is the eigendecomposition of the matrix T' * T?
//
  w = r8mat_mtm_new ( n, n, n, t, t );

  it_max = 50;
  v = new double[n*n];
  lambda = new double[n];

  jacobi_eigenvalue ( n, w, it_max, v, lambda, it_num, rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
//
//  Free memory.
//
  delete [] lambda;
  delete [] t;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST03 compares the unit Wishart and Bartlett sample matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double diff;
  int n;
  double *t;
  double *tt;
  double *w;
//
//   Initialize the RNGLIB library.
//   Normally, we would do this just once, here at the beginning.
//   In this example, however, we really want to do it just before
//   we call each of the sampling routines, so that they both access
//   the same set of random numbers...
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST03:\n";
  cout << "  Verify that, if using the same set of random numbers,\n";
  cout << "    W = T' * T,\n";
  cout << "  where\n";
  cout << "    W = wishart_unit_sample ( n, df );\n";
  cout << "    T = bartlett_unit_sample ( n, df );\n";
//
//   Set the parameters.
//
  n = 5;
  df = 8;
//
//   Initialize the random number package and compute W.
//
  initialize ( );
  w = wishart_unit_sample ( n, df );
//
//   Initialize the random number package again, and compute T.
//
  initialize ( );
  t = bartlett_unit_sample ( n, df );
//
//   Compute T' * T.
//
  tt = r8mat_mtm_new ( n, n, n, t, t );
//
//   Compare T'T to W.
//
  diff = r8mat_norm_fro_affine ( n, n, w, tt );
  cout << "\n";
  cout << "  Frobenius norm of error is " << diff << "\n";
//
//  Free memory.
//
  delete [] t;
  delete [] tt;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST04 demonstrates the Wishart sampling function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  int i;
  int it_max;
  int it_num;
  int j;
  double *lambda;
  int n;
  int rot_num;
//
//  Note that R is an upper triangular matrix,
//  whose entries here are listed in column major order.
//
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  double *sigma;
  double sigma_diag[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  double *v;
  double *w;
//
//   Initialize the RNGLIB library.
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST04:\n";
  cout << "  We can compute sample Wishart matrices by:\n";
  cout << "    W = wishart_sample ( n, df, sigma );\n";
//
//   Set the parameters and call.
//
  n = 5;
  df = 8;
  sigma = r8mat_identity_new ( n );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, Identity ):" );
  delete [] w;
//
//   Calling again yields a new matrix.
//
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, Identity ):" );
  delete [] sigma;
  delete [] w;
//
//   Try a diagonal matrix.
//
  sigma = r8mat_diagonal_new ( n, sigma_diag );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 5, 8, diag(1,2,3,4,5) ):" );
  delete [] sigma;
  delete [] w;
//
//   Try a smaller matrix.  Sigma must be positive definite symmetric.
//
  n = 3;
  df = 3;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Set covariance SIGMA:" );
  w = wishart_sample ( n, df, sigma );
  r8mat_print ( n, n, w, "  wishart_sample ( 3, 3, sigma ):" );
//
//   What is the eigendecomposition of this matrix?
//
  it_max = 50;
  v = new double[n*n];
  lambda = new double[n];

  jacobi_eigenvalue ( n, w, it_max, v, lambda, it_num, rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
//
//  Free memory.
//
  delete [] lambda;
  delete [] sigma;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST05 demonstrates the Bartlett sampling function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  int it_max;
  int it_num;
  double *lambda;
  int n;
//
//  Note that R is an upper triangular matrix,
//  whose entries here are listed in column major order.
//
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  int rot_num;
  double *sigma;
  double sigma_diag[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };
  double *t;
  double *v;
  double *w;
//
//   Initialize the RNGLIB library.
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST05:\n";
  cout << "  We can compute sample Bartlett matrices by:\n";
  cout << "    T = bartlett_sample ( n, df, sigma );\n";
//
//   Set the parameters and call.
//
  n = 5;
  df = 8;
  sigma = r8mat_identity_new ( n );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, Identity ):" );
  delete [] t;
//
//   Calling again yields a new matrix.
//
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, Identity ):" );
  delete [] sigma;
  delete [] t;
//
//   Try a diagonal matrix.
//
  sigma = r8mat_diagonal_new ( n, sigma_diag );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 5, 8, diag(1,2,3,4,5) ):" );
  delete [] sigma;
  delete [] t;
//
//   Try a smaller matrix.
//
  n = 3;
  df = 3;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Set covariance SIGMA:" );
  t = bartlett_sample ( n, df, sigma );
  r8mat_print ( n, n, t, "  bartlett_sample ( 3, 3, sigma ):" );
//
//   What is the eigendecomposition of T' * T?
//
  w = r8mat_mtm_new ( n, n, n, t, t );
  it_max = 50;
  v = new double[n*n];
  lambda = new double[n];

  jacobi_eigenvalue ( n, w, it_max, v, lambda, it_num, rot_num );
  r8mat_print ( n, n, v, "  Eigenvectors of previous matrix:" );
  r8vec_print ( n, lambda, "  Eigenvalues of previous matrix:" );
//
//  Free memory.
//
  delete [] lambda;
  delete [] sigma;
  delete [] t;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST06 compares the Wishart and Bartlett sample matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double diff;
  int n;
//
//  Note that R is an upper triangular matrix,
//  whose entries here are listed in column major order.
//
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  double *sigma;
  double *t;
  double *tt;
  double *w;
//
//   Initialize the RNGLIB library.
//   Normally, we would do this just once, here at the beginning.
//   In this example, however, we really want to do it just before
//   we call each of the sampling routines, so that they both access
//   the same set of random numbers...
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST06:\n";
  cout << "  Verify that, if using the same set of random numbers,\n";
  cout << "    W = T'' * T,\n";
  cout << "  where\n";
  cout << "    W = wishart_sample ( n, df, sigma );\n";
  cout << "    T = bartlett_sample ( n, df, sigma );\n";
//
//   Set the parameters.
//
  n = 3;
  df = 5;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Covariance SIGMA:" );
//
//   Initialize the random number package and compute W.
//
  initialize ( );
  w = wishart_sample ( n, df, sigma );
//
//   Initialize the random number package again, and compute T.
//
  initialize ( );
  t = bartlett_sample ( n, df, sigma );
//
//   Compute T' * T.
//
  tt = r8mat_mtm_new ( n, n, n, t, t );
//
//   Compare T'T to W.
//
  diff = r8mat_norm_fro_affine ( n, n, w, tt );
  cout << "\n";
  cout << "  Frobenius norm of error is " << diff << "\n";
//
//  Free memory.
//
  delete [] sigma;
  delete [] t;
  delete [] tt;
  delete [] w;

  return;
}
//****************************************************************************80

void wishart_test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST07 demonstrates a property of the Wishart distribution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double diff;
  double divisor;
  int i;
  int n;
//
//  Note that R is an upper triangular matrix,
//  whose entries here are listed in column major order.
//
  double r[3*3] = { 
    5.0, 0.0, 0.0,
    1.0, 4.0, 0.0,
    3.0, 2.0, 6.0 };
  int sample_num;
  double *sigma;
  double *w;
  double *w_average;
//
//   Initialize the RNGLIB library.
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST07:\n";
  cout << "  For given values of N, DF, SIGMA, the random\n";
  cout << "  matrices from the Wishart distribution:\n";
  cout << "    W = wishart_sample ( n, df, sigma );\n";
  cout << "  should have mean DF * SIGMA.\n";
//
//   Set the parameters.
//
  n = 3;
  cout << "  Fix N = " << n << "\n";
  df = 5;
  cout << "  Fix DF = " << df << "\n";
  sigma = r8mat_mtm_new ( n, n, n, r, r );
  r8mat_print ( n, n, sigma, "  Fix covariance SIGMA:" );
//
//   Sample many times and average.
//
  sample_num = 1000;
  w_average = r8mat_zero_new ( n, n );
  for ( i = 1; i <= sample_num; i++ )
  {
    w = wishart_sample ( n, df, sigma );
    r8mat_add ( n, n, w, w_average );
    delete [] w;
  }
  divisor = ( double ) sample_num;
  r8mat_divide ( n, n, divisor, w_average );
//
//   Compare SIGMA and W_SAMPLE / DF.
//
  divisor = ( double ) df;
  r8mat_divide ( n, n, divisor, w_average );

  r8mat_print ( n, n, w_average, "  W_Average / DF: " );

  diff = r8mat_norm_fro_affine ( n, n, sigma, w_average );
  cout << "\n";
  cout << "  Frobenius norm of SIGMA-W_average/DF = " << diff << "\n";
//
//  Free memory.
//
  delete [] sigma;
  delete [] w_average;

  return;
}
//****************************************************************************80

void wishart_test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST08 samples the unit Wishart and unit Wishart inverse matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double diff;
  int i;
  double *ident;
  int j;
  double *m;
  int n;
  double *w;
  double *wm;
//
//   Initialize the RNGLIB library.
//   Normally, we would do this just once, here at the beginning.
//   In this example, however, we really want to do it just before
//   we call each of the sampling routines, so that they both access
//   the same set of random numbers...
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST08:\n";
  cout << "  Verify that, if using the same set of random numbers,\n";
  cout << "    inverse(W) = M,\n";
  cout << "  where\n";
  cout << "    W = wishart_unit_sample ( n, df );\n";
  cout << "    M = wishart_unit_sample_inverse ( n, df );\n";
//
//   Set the parameters.
//
  n = 5;
  df = 8;
//
//   Initialize the random number package and compute W.
//
  initialize ( );
  w = wishart_unit_sample ( n, df );
//
//   Initialize the random number package again, and compute T.
//
  initialize ( );
  m = wishart_unit_sample_inverse ( n, df );
//
//   Compute W * M.
//
  wm = r8mat_mm_new ( n, n, n, w, m );
//
//   Compare M * W to I.
//
  ident = r8mat_identity_new ( n );
  diff = r8mat_norm_fro_affine ( n, n, wm, ident );
  cout << "\n";
  cout << "  Frobenius norm of error is " << diff << "\n";
//
//  Free memory.
//
  delete [] ident;
  delete [] m;
  delete [] w;
  delete [] wm;

  return;
}
//****************************************************************************80

void wishart_test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    WISHART_TEST09 samples the Wishart and Wishart inverse matrices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int df;
  double diff;
  int i;
  double *ident;
  int j;
  double *m;
  int n;
//
//  Note that R is an upper triangular matrix,
//  whose entries here are listed in column major order.
//
  double r[5*5] = { 
    3.0, 0.0, 0.0, 0.0, 0.0, 
    1.0, 7.0, 0.0, 0.0, 0.0, 
    1.0, 1.0, 5.0, 0.0, 0.0, 
    1.0, 2.0, 1.0, 4.0, 0.0, 
    1.0, 3.0, 3.0, 2.0, 6.0 };
  double *sigma;
  double *w;
  double *wm;
//
//   Initialize the RNGLIB library.
//   Normally, we would do this just once, here at the beginning.
//   In this example, however, we really want to do it just before
//   we call each of the sampling routines, so that they both access
//   the same set of random numbers...
//
  initialize ( );

  cout << "\n";
  cout << "WISHART_TEST09:\n";
  cout << "  Verify that, if using the same set of random numbers,\n";
  cout << "    inverse(W) = M,\n";
  cout << "  where\n";
  cout << "    W = wishart_sample ( n, df, sigma );\n";
  cout << "    M = wishart_sample_inverse ( n, df, sigma );\n";
//
//   Set the parameters.
//
  n = 5;
  df = 8;
  sigma = r8mat_mtm_new ( n, n, n, r, r );
//
//   Initialize the random number package and compute W.
//
  initialize ( );
  w = wishart_sample ( n, df, sigma );
//
//   Initialize the random number package again, and compute T.
//
  initialize ( );
  m = wishart_sample_inverse ( n, df, sigma );
//
//   Compute W * M.
//
  wm = r8mat_mm_new ( n, n, n, w, m );
//
//   Compare M * W to I.
//
  ident = r8mat_identity_new ( n );
  diff = r8mat_norm_fro_affine ( n, n, wm, ident );
  cout << "\n";
  cout << "  Frobenius norm of error is " << diff << "\n";
//
//  Free memory.
//
  delete [] ident;
  delete [] m;
  delete [] sigma;
  delete [] w;
  delete [] wm;

  return;
}
