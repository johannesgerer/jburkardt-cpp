# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "haar.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HAAR_PRB.
//
//  Discussion:
//
//    HAAR_PRB tests the HAAR library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "HAAR_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the HAAR library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HAAR_PRB\n";
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
//    TEST01 tests HAAR_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  HAAR_1D computes the Haar transform of a vector.\n";
//
//  Random data.
//
  n = 16;
  seed = 123456789;
  u = r8vec_uniform_01_new ( n, &seed );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  cout << "\n";
  cout << "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }
  delete [] u;
  delete [] v;
  delete [] w;
//
//  Constant signal.
//
  n = 8;
  u = r8vec_ones_new ( n );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  cout << "\n";
  cout << "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }
  delete [] u;
  delete [] v;
  delete [] w;
//
//  Linear signal.
//
  n = 16;
  u = r8vec_linspace_new ( n, 1.0, ( double ) n );
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  cout << "\n";
  cout << "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }
  delete [] u;
  delete [] v;
  delete [] w;
//
//  Quadratic data.
//
  n = 8;
  u = ( double * ) malloc ( n * sizeof ( double ) );
  u[0] = 25.0;
  u[1] = 16.0;
  u[2] = 9.0;
  u[3] = 4.0;
  u[4] = 1.0;
  u[5] = 0.0;
  u[6] = 1.0;
  u[7] = 4.0;
  v = r8vec_copy_new ( n, u );

  haar_1d ( n, v );

  w = r8vec_copy_new ( n, v );
  haar_1d_inverse ( n, w );

  cout << "\n";
  cout << "   i      U(i)        H(U)(i)  Hinv(H(U))(i)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(10) << u[i]
         << "  " << setw(10) << v[i]
         << "  " << setw(10) << w[i] << "\n";
  }
  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HAAR_2D and HAAR_2D_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m = 16;
  int n = 4;
  int seed;
  double *u;
  double *v;
  double *w;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HAAR_2D computes the Haar transform of an array.\n";
  cout << "  HAAR_2D_INVERSE inverts the transform.\n";
//
//  Demonstrate successful inversion.
//
  seed = 123456789;
  u = r8mat_uniform_01_new ( m, n, &seed );

  r8mat_print ( m, n, u, "  Input array U:" );

  v = r8mat_copy_new ( m, n, u );

  haar_2d ( m, n, v );

  r8mat_print ( m, n, v, "  Transformed array V:" );

  w = r8mat_copy_new ( m, n, v );

  haar_2d_inverse ( m, n, w );

  r8mat_print ( m, n, w, "  Recovered array W:" );

  delete [] u;
  delete [] v;
  delete [] w;

  return;
}
