# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "sphere_stereograph.hpp"

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
//    MAIN is the main program for SPHERE_STEREOGRAPH_PRB.
//
//  Discussion:
//
//    SPHERE_STEREOGRAPH_PRB tests the SPHERE_STEREOGRAPH library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SPHERE_STEREOGRAPH_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the SPHERE_STEREOGRAPH library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SPHERE_STEREOGRAPH_PRB\n";
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
//    TEST01 checks that the two functions are inverses.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  double dif;
  int i;
  int j;
  int m;
  int n;
  double *p;
  double *p1;
  double *p2;
  double *q;
  double *q1;
  double *q2;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SPHERE_STEREOGRAPH maps from sphere to plane.\n";
  cout << "  SPHERE_STEREOGRAPH_INVERSE is the inverse map.\n";
  cout << "  Check that these two functions are inverses.\n";

  m = 3;
  n = 100;
//
//  Check #1.
//
  seed = 123456789;
  p1 = uniform_on_sphere01_map ( m, n, &seed );
  q  = sphere_stereograph ( m, n, p1 );
  p2 = sphere_stereograph_inverse ( m, n, q );
  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( p1[i+j*m] - p2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from sphere to plane to sphere.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] p1;
  delete [] p2;
  delete [] q;
//
//  Check #2.
//
  q1 = r8mat_uniform_01_new ( m, n, &seed );
  for ( j = 0; j < n; j++ )
  {
    q1[m-1+j*m] = 1.0;
  }

  p = sphere_stereograph_inverse ( m, n, q1 );
  q2 = sphere_stereograph ( m, n, p );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( q1[i+j*m] - q2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from plane to sphere to plane.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] p;
  delete [] q1;
  delete [] q2;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 checks the generalized mapping.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *center;
  double dif;
  double *focus;
  int i;
  int j;
  int m;
  int n;
  double *p1;
  double *p2;
  double *q1;
  double *q2;
  int seed;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  SPHERE_STEREOGRAPH standard mapping from sphere to plane.\n";
  cout << "  SPHERE_STEREOGRAPH2 generalized mapping:\n";
  cout << "  (focus and center are arbitrary)\n";
  cout << "  Check that these two functions can agree.\n";

  m = 3;
  n = 100;

  focus = new double[m];
  for ( i = 0; i < m - 1; i++ )
  {
    focus[i] = 0.0;
  }
  focus[m-1] = -1.0;

  center = new double[m];
  for ( i = 0; i < m; i++ )
  {
    center[i] = 0.0;
  }
//
//  Check #1.
//
  seed = 123456789;
  p1 = uniform_on_sphere01_map ( m, n, &seed );

  q1 = sphere_stereograph ( m, n, p1 );

  q2 = sphere_stereograph2 ( m, n, p1, focus, center );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( q1[i+j*m] - q2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from sphere to plane.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] p1;
  delete [] q1;
  delete [] q2;
//
//  Check #2.
//
  q1 = r8mat_uniform_01_new ( m, n, &seed );
  for ( j = 0; j < n; j++ )
  {
    q1[m-1+j*m] = 1.0;
  }

  p1 = sphere_stereograph_inverse ( m, n, q1 );

  p2 = sphere_stereograph2_inverse ( m, n, q1, focus, center );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( p1[i+j*m] - p2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from plane to sphere.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] p1;
  delete [] p2;
  delete [] q1;

  delete [] center;
  delete [] focus;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 checks that the two functions are inverses.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *alpha;
  double *beta;
  double *center;
  double dif;
  double *focus;
  int i;
  int j;
  int m;
  int n;
  double *normal;
  double *p;
  double *p1;
  double *p2;
  double *pq;
  double *pr;
  double *q;
  double *q1;
  double *q2;
  double r;
  int seed;
  double *tang;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  SPHERE_STEREOGRAPH2 maps from sphere to plane.\n";
  cout << "  SPHERE_STEREOGRAPH2_INVERSE is the inverse map.\n";
  cout << "  Check that these two functions are inverses.\n";

  m = 3;
  n = 100;
  seed = 123456789;

  focus = r8vec_uniform_01_new ( m, &seed );
  center = r8vec_uniform_01_new ( m, &seed );
  r = r8vec_norm_affine ( m, focus, center );

  cout << "\n";
  cout << "  Using radius = " << r << "\n";
  r8vec_transpose_print ( m, center, "  Center:" );
  r8vec_transpose_print ( m, focus, "  Focus:" );
//
//  Check #1.
//
  p1 = uniform_on_sphere01_map ( m, n, &seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
       p1[i+j*m] = center[i] + r * p1[i+j*m];
    }
  }
  q  = sphere_stereograph2 ( m, n, p1, focus, center );

  p2 = sphere_stereograph2_inverse ( m, n, q, focus, center );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( p1[i+j*m] - p2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from sphere to plane to sphere.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] p1;
  delete [] p2;
  delete [] q;
//
//  Check #2.
//  We have to work hard to get random points on the plane, since
//  all we know to begin with is the point of tangency and the normal.
//
  tang = new double[m];
  for ( i = 0; i < m; i++ )
  {
    tang[i] = 2.0 * center[i] - focus[i];
  }
  normal = new double[m];
  for ( i = 0; i < m; i++ )
  {
    normal[i] = center[i] - focus[i];
  }

  pr = new double[m];
  pq = new double[m];
  plane_normal_basis_3d ( tang, normal, pr, pq );

  q1 = new double[m*n];

  alpha = r8vec_uniform_01_new ( n, &seed );
  beta = r8vec_uniform_01_new ( n, &seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      q1[i+j*m] = tang[i] + pr[i] * alpha[j] + pq[i] * beta[j];
    }
  }
  p = sphere_stereograph2_inverse ( m, n, q1, focus, center );
  q2 = sphere_stereograph2 ( m, n, p, focus, center );

  dif = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dif = dif + pow ( q1[i+j*m] - q2[i+j*m], 2 );
    }
  }
  dif = sqrt ( dif ) / ( double ) ( n );
  cout << "\n";
  cout << "  Map points from plane to sphere to plane.\n";
  cout << "  RMS difference for " << n << " points was " << dif << "\n";

  delete [] alpha;
  delete [] beta;
  delete [] normal;
  delete [] p;
  delete [] pq;
  delete [] pr;
  delete [] q1;
  delete [] q2;
  delete [] tang;

  return;
}
