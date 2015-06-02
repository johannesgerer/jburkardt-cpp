# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "jacobi_eigenvalue.hpp"

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
//    MAIN is the main program for JACOBI_EIGENVALUE_PRB.
//
//  Discussion:
//
//    JACOBI_EIGENVALUE_PRB tests the JACOBI_EIGENVALUE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "JACOBI_EIGENVALUE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the JACOBI_EIGENVALUE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "JACOBI_EIGENVALUE_PRB\n";
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
//    TEST01 uses a 4x4 test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N] = {
      4.0,  -30.0,    60.0,   -35.0, 
    -30.0,  300.0,  -675.0,   420.0, 
     60.0, -675.0,  1620.0, -1050.0, 
    -35.0,  420.0, -1050.0,   700.0 };
  double d[N];
  double error_frobenius;
  int it_max;
  int it_num;
  int n = N;
  int rot_num;
  double v[N*N];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  For a symmetric matrix A,\n";
  cout << "  JACOBI_EIGENVALUE computes the eigenvalues D\n";
  cout << "  and eigenvectors V so that A * V = D * V.\n";

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num );

  cout << "\n";
  cout << "  Number of iterations = " << it_num << "\n";
  cout << "  Number of rotations  = " << rot_num << "\n";

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
//
//  Compute eigentest.
//
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  cout << "\n";
  cout << "  Frobenius norm error in eigensystem A*V-D*V = " 
       << error_frobenius << "\n";

  return;
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 uses a 4x4 test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
{
# define N 4

  double a[N*N] = {
    4.0, 0.0, 0.0, 0.0, 
    0.0, 1.0, 0.0, 0.0, 
    0.0, 0.0, 3.0, 0.0, 
    0.0, 0.0, 0.0, 2.0 };
  double d[N];
  double error_frobenius;
  int it_max;
  int it_num;
  int n = N;
  int rot_num;
  double v[N*N];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For a symmetric matrix A,\n";
  cout << "  JACOBI_EIGENVALUE computes the eigenvalues D\n";
  cout << "  and eigenvectors V so that A * V = D * V.\n";
  cout << "\n";
  cout << "As a sanity check, input a diagonal matrix.\n";

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num );

  cout << "\n";
  cout << "  Number of iterations = " << it_num << "\n";
  cout << "  Number of rotations  = " << rot_num << "\n";

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
//
//  Compute eigentest.
//
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  cout << "\n";
  cout << "  Frobenius norm error in eigensystem A*V-D*V = " 
       << error_frobenius << "\n";

  return;
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 uses a 5x5 test matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 July 2013
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double a[N*N];
  double d[N];
  double error_frobenius;
  int i;
  int it_max;
  int it_num;
  int j;
  int n = N;
  int rot_num;
  double v[N*N];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For a symmetric matrix A,\n";
  cout << "  JACOBI_EIGENVALUE computes the eigenvalues D\n";
  cout << "  and eigenvectors V so that A * V = D * V.\n";
  cout << "\n";
  cout << "  Use the discretized second derivative matrix.\n";

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[i+j*n] = -2.0;
      }
      else if ( i == j + 1 || i == j - 1 )
      {
        a[i+j*n] = 1.0;
      }
      else
      {
        a[i+j*n] = 0.0;
      }
    }
  }

  r8mat_print ( n, n, a, "  Input matrix A:" );

  it_max = 100;

  jacobi_eigenvalue ( n, a, it_max, v, d, it_num, rot_num );

  cout << "\n";
  cout << "  Number of iterations = " << it_num << "\n";
  cout << "  Number of rotations  = " << rot_num << "\n";

  r8vec_print ( n, d, "  Eigenvalues D:" );

  r8mat_print ( n, n, v, "  Eigenvector matrix V:" );
//
//  Compute eigentest.
//
  error_frobenius = r8mat_is_eigen_right ( n, n, a, v, d );
  cout << "\n";
  cout << "  Frobenius norm error in eigensystem A*V-D*V = " 
       << error_frobenius << "\n";

  return;
# undef N
}
