# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_triangle_monte_carlo.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SPHERE_TRIANGLE_MONTE_CARLO_PRB.
//
//  Discussion:
//
//    SPHERE_TRIANGLE_MONTE_CARLO_PRB tests the SPHERE_TRIANGLE_MONTE_CARLO 
//    library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_TRIANGLE_MONTE_CARLO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_TRIANGLE_MONTE_CARLO library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_TRIANGLE_MONTE_CARLO_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//*****************************************************************************/

void test01 ( )

//*****************************************************************************/
//
//  Purpose:
//
//    TEST01 uses SPHERE_TRIANGLE_SAMPLE_01 with an increasing number of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define M 3

  double area;
  int e[M];
  int e_test[M*7] = {
    0, 0, 0, 
    2, 0, 0, 
    0, 2, 0, 
    0, 0, 2, 
    4, 0, 0, 
    2, 2, 0, 
    0, 0, 4 };
  double error;
  int i;
  int j;
  int k;
  int m = M;
  int n;
  double r8_pi = 3.1415926535897932384626434;
  double result;
  int seed;
  double shrink;
  double v1[M];
  double v2[M];
  double v3[M];
  double wc[M];
  double *w1;
  double *w2;
  double *w3;
  double *value;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Estimate monomial integrals over a sphere triangle\n";
  cout << "  using the Monte Carlo method.\n";

  seed = 123456789;
//
//  Choose three points at random to define a spherical triangle.
//
  w1 = sphere01_sample ( 1, seed );
  w2 = sphere01_sample ( 1, seed );
  w3 = sphere01_sample ( 1, seed );

  for ( i = 0; i < m; i++ )
  {
    wc[i] = ( w1[i] + w2[i] + w3[i] ) / 3.0;
  }
  r8vec_normalize ( m, wc );
//
//  Shrink triangle by factor F.
//
  shrink = 2.0;

  for ( k = 1; k <= 3; k++ )
  {
    shrink = shrink / 2.0;

    for ( i = 0; i < m; i++ )
    {
      v1[i] = wc[i] + shrink * ( w1[i] - wc[i] );
      v2[i] = wc[i] + shrink * ( w2[i] - wc[i] );
      v3[i] = wc[i] + shrink * ( w3[i] - wc[i] );
    }
    r8vec_normalize ( m, v1 );
    r8vec_normalize ( m, v2 );
    r8vec_normalize ( m, v3 );

    area = sphere01_triangle_vertices_to_area ( v1, v2, v3 );

    cout << "\n";
    cout << "  Vertices of random spherical triangle\n";
    cout << "  with shrink factor = " << shrink << "\n";
    cout << "  and area = " << area << "\n";
    cout << "\n";
    r8vec_transpose_print ( m, v1, "  V1:" );
    r8vec_transpose_print ( m, v2, "  V2:" );
    r8vec_transpose_print ( m, v3, "  V3:" );
//
//  Estimate integrals.
//
    cout << "\n";
    cout << "         N        1              X^2             Y^2";
    cout << "             Z^2             X^4           X^2Y^2           Z^4\n";
    cout << "\n";

    n = 1;

    while ( n <= 4 * 65536 )
    {
      x = sphere01_triangle_sample ( n, v1, v2, v3, seed );

      cout << "  " << setw(8) << n;
      for ( j = 0; j < 7; j++ )
      {
        for ( i = 0; i < m; i++ )
        {
          e[i] = e_test[i+j*m];
        }
        value = monomial_value ( m, n, e, x );

        result = area * r8vec_sum ( n, value ) / ( double ) ( n );
        cout << "  " << setw(14) << result;
      }

      cout << "\n";

      delete [] value;
      delete [] x;

      n = 2 * n;
    }
  }

  delete [] w1;
  delete [] w2;
  delete [] w3;

  return;
# undef M
}
