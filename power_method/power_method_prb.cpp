# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "power_method.hpp"

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
//    MAIN is the main program for POWER_METHOD_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 20008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POWER_METHOD_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the POWER_METHOD library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POWER_METHOD_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 uses POWER_METHOD on the Fibonacci2 matrix.
//
//  Discussion:
//
//    This matrix, despite having a single dominant eigenvalue, will generally
//    converge only very slowly under the power method.  This has to do with
//    the fact that the matrix has only 3 eigenvectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 July 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double cos_x1x2;
  double ctime;
  double ctime1;
  double ctime2;
  int i;
  int it_max;
  int it_num;
  double lambda;
  int n = 50;
  double norm;
  double phi;
  int seed;
  double sin_x1x2;
  double tol;
  double *x;
  double *x2;

  a = fibonacci2 ( n );

  seed = 123456789;
  x = r8vec_uniform_01 ( n, &seed );

  it_max = 300;
  tol = 0.000001;

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use POWER_METHOD on the Fibonacci2 matrix.\n";
  cout << "\n";
  cout << "  Matrix order N       = " << n << "\n";
  cout << "  Maximum iterations   = " << it_max << "\n";
  cout << "  Error tolerance      = " << tol << "\n";

  ctime1 = cpu_time ( );

  power_method ( n, a, x, it_max, tol, &lambda, &it_num );

  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  cout << "\n";
  cout << "  Number of iterations = " << it_num << "\n";
  cout << "  CPU time             = " << ctime << "\n";
  cout << "  Estimated eigenvalue = " << setprecision(14) << lambda << "\n";
  cout << "  Correct value        = " << setprecision(14) << phi << "\n";
  cout << "  Error                = " << r8_abs ( lambda - phi ) << "\n";
//
//  X2 is the exact eigenvector.
//
  x2 = new double[n];

  x2[0] = 1.0;
  for ( i = 1; i < n; i++ )
  {
    x2[i] = phi * x2[i-1];
  }
  norm = r8vec_norm_l2 ( n, x2 );
  for ( i = 0; i < n; i++ )
  {
    x2[i] = x2[i] / norm;
  }
//
//  The sine of the angle between X and X2 is a measure of error.
//
  cos_x1x2 = r8vec_dot ( n, x, x2 );
  sin_x1x2 = sqrt ( ( 1.0 - cos_x1x2 ) * ( 1.0 + cos_x1x2 ) );

  cout << "\n";
  cout << "  Sine of angle between true and estimated vectors = " << sin_x1x2 << "\n";

  delete [] a;
  delete [] x;
  delete [] x2;

  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses POWER_METHOD2 on the Fibonacci2 matrix.
//
//  Discussion:
//
//    This matrix, despite having a single dominant eigenvalue, will generally
//    converge only very slowly under the power method.  This has to do with
//    the fact that the matrix has only 3 eigenvectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double ctime;
  double ctime1;
  double ctime2;
  int i;
  int it_max;
  int it_num;
  complex <double> lambda;
  int n = 50;
  double phi;
  int seed;
  double tol;
  complex <double> *v;
  double *x;

  a = fibonacci2 ( n );
  v = new complex<double> [n];

  seed = 123456789;
  x = r8vec_uniform_01 ( n, &seed );

  it_max = 300;
  tol = 0.000001;

  phi = ( 1.0 + sqrt ( 5.0 ) ) / 2.0;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use POWER_METHOD2 on the Fibonacci2 matrix.\n";
  cout << "\n";
  cout << "  Matrix order N       = " << n << "\n";
  cout << "  Maximum iterations   = " << it_max << "\n";
  cout << "  Error tolerance      = " << tol << "\n";

  ctime1 = cpu_time ( );

  power_method2 ( n, a, x, it_max, tol, &lambda, v, &it_num );

  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  cout << "\n";
  cout << "  Number of iterations = " << it_num << "\n";
  cout << "  CPU time             = " << ctime << "\n";
  cout << "  Estimated eigenvalue = " 
       << "  " << setprecision(14) << real ( lambda )
       << "  " << setprecision(14) << imag ( lambda ) << "\n";
  cout << "  Correct value        = " << setprecision(14) << phi << "\n";
  cout << "  Error                = " << abs ( lambda - phi ) << "\n";

  delete [] a;
  delete [] v;
  delete [] x;

  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses POWER_METHOD2 on the TRIS matrix.
//
//  Discussion:
//
//    This matrix, despite having a single dominant eigenvalue, will generally
//    converge only very slowly under the power method.  This has to do with
//    the fact that the matrix has only 3 eigenvectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double alpha;
  double beta;
  double ctime;
  double ctime1;
  double ctime2;
  double gamma;
  int i;
  int it_max;
  int it_num;
  complex <double> lambda;
  complex <double> lambda_max;
  complex <double> *lambda_vec;
  int n = 50;
  int seed;
  double tol;
  complex <double> *v;
  double *x;

  alpha = -1.0;
  beta = 10.0;
  gamma = 8.0;

  a = tris ( n, n, alpha, beta, gamma );

  v = new complex<double> [n];

  seed = 123456789;
  x = r8vec_uniform_01 ( n, &seed );

  it_max = 4000;
  tol = 0.000001;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Use POWER_METHOD2 on the TRIS (tridiagonal scalar) matrix.\n";
  cout << "\n";
  cout << "  Matrix order N         = " << n << "\n";
  cout << "  Maximum iterations     = " << it_max << "\n";
  cout << "  Error tolerance        = " << tol << "\n";

  ctime1 = cpu_time ( );

  power_method2 ( n, a, x, it_max, tol, &lambda, v, &it_num );

  ctime2 = cpu_time ( );
  ctime = ctime2 - ctime1;

  cout << "\n";
  cout << "  Number of iterations   = " << it_num << "\n";
  cout << "  CPU time               = " << ctime << "\n";
  cout << "  Estimated eigenvalue   = " 
       << setprecision(14) << real ( lambda )
       << "  " << setprecision(14) << imag ( lambda ) << "\n";

  lambda_vec = tris_eigenvalues ( n, alpha, beta, gamma );

  lambda_max = lambda_vec[0];
  for ( i = 1; i < n; i++ )
  {
    if ( abs ( lambda_max ) < abs ( lambda_vec[i] ) )
    {
      lambda_max = lambda_vec[i];
    }
  }

  cout << "  Correct max eigenvalue = " 
       << setprecision(14) << real ( lambda_max )
       << "  " << setprecision(14) << imag ( lambda_max ) << "\n";

  cout << "  Error                  = " << abs ( lambda - lambda_max ) << "\n";

  delete [] a;
  delete [] lambda_vec;
  delete [] v;
  delete [] x;

  return;
}
