# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "lagrange_interp_nd.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
double *f_sinr ( int m, int n, double x[] );
double *f_poly352 ( int m, int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LAGRANGE_INTERP_ND_PRB.
//
//  Discussion:
//
//    LAGRANGE_INTERP_ND_PRB tests the LAGRANGE_INTERP_ND library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LAGRANGE_INTERP_ND_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the LAGRANGE_INTERP_ND library.\n";
//
//  Use the interface that passes in the orders directly.
//
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Repeat tests 1, 2, 3, and 4,
//  using the interface that passes in the orders indirectly,
//  based on the "level".
//
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Experiment with anisotropic orders.
//
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LAGRANGE_INTERP_ND_PRB:\n";
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
//    TEST01 interpolates in 1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int j;
  int m;
  int *n_1d;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Interpolate in 1D, using orders.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE evaluates it.\n";

  m = 1;

  n_1d = new int[m];
  for ( i = 0; i < m; i++ )
  {
    n_1d[i] = 5;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size ( m, n_1d );
  xd = lagrange_interp_nd_grid ( m, n_1d, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] n_1d;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80


void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 interpolates in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int j;
  int m;
  int *n_1d;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Interpolate in 2D, using orders.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE evaluates it.\n";

  m = 2;

  n_1d = new int[m];
  for ( i = 0; i < m; i++ )
  {
    n_1d[i] = 5;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size ( m, n_1d );
  xd = lagrange_interp_nd_grid ( m, n_1d, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] n_1d;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80


void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 interpolates in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int j;
  int m;
  int *n_1d;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Interpolate in 3D, using orders.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE evaluates it.\n";

  m = 3;

  n_1d = new int[m];
  for ( i = 0; i < m; i++ )
  {
    n_1d[i] = 5;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size ( m, n_1d );
  xd = lagrange_interp_nd_grid ( m, n_1d, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] n_1d;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 interpolates in 3D, using increasing resolution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e;
  int i;
  int j;
  int l;
  int m;
  int *n_1d;
  int nd;
  int ni;
  int order;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Interpolate in 3D, using orders.\n";
  cout << "  Use a sequence of increasing orders.\n";

  m = 3;

  n_1d = new int[m];

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  ni = 20;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );

  cout << "\n";
  cout << "  Level     Order   Average Error\n";
  cout << "\n";

  for ( l = 1; l <= 5; l++ )
  {
    if ( l == 0 )
    {
      order = 1;
    }
    else
    {
      order = i4_power ( 2, l ) + 1;
    }

    for ( i = 0; i < m; i++ )
    {
      n_1d[i] = order;
    }
 
    nd = lagrange_interp_nd_size ( m, n_1d );
    xd = lagrange_interp_nd_grid ( m, n_1d, a, b, nd );
    zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
    zi = lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi );

    e = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

    cout << "  " << setw(5) << l
         << "  " << setw(5) << nd
         << "  " << setw(10) << e << "\n";

    delete [] xd;
    delete [] zd;
    delete [] zi;
  }

  delete [] a;
  delete [] b;
  delete [] n_1d;
  delete [] xi;
  delete [] ze;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 repeats test 1 using levels.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int *ind;
  int j;
  int m;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST05:\n";
  cout << "  Repeat test #1, using levels.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE2 evaluates it.\n";

  m = 1;

  ind = new int[m];
  for ( i = 0; i < m; i++ )
  {
    ind[i] = 2;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size2 ( m, ind );
  xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] ind;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 repeats test 2 using levels.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int *ind;
  int j;
  int m;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST06:\n";
  cout << "  Repeat test #2, using levels.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE2 evaluates it.\n";

  m = 2;

  ind = new int[m];
  for ( i = 0; i < m; i++ )
  {
    ind[i] = 2;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size2 ( m, ind );
  xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] ind;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 repeats test 3 using levels.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int *ind;
  int j;
  int m;
  int nd;
  int ni;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST07:\n";
  cout << "  Repeat test #3,  using levels.\n";
  cout << "  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.\n";
  cout << "  LAGRANGE_INTERP_ND_VALUE2 evaluates it.\n";

  m = 3;

  ind = new int[m];
  for ( i = 0; i < m; i++ )
  {
    ind[i] = 2;
  }

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  nd = lagrange_interp_nd_size2 ( m, ind );
  xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
  zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
  cout << "\n";
  cout << "         Zinterp          Zexact      Error\n";
  cout << "\n";

  ni = 5;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );
  zi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );

  for ( j = 0; j < ni; j++ )
  {
    cout << "  " << setw(14) << zi[j]
         << "  " << setw(14) << ze[j]
         << "  " << setw(10) << r8_abs ( zi[j] - ze[j] ) << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] ind;
  delete [] xd;
  delete [] xi;
  delete [] zd;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 repeats test 4 using levels.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e;
  int i;
  int *ind;
  int j;
  int l;
  int m;
  int nd;
  int ni;
  int order;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST08:\n";
  cout << "  Interpolate in 3D, using levels.\n";
  cout << "  Use a sequence of increasing levels.\n";

  m = 3;

  ind = new int[m];

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  ni = 20;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_sinr ( m, ni, xi );

  cout << "\n";
  cout << "  Level     Order   Average Error\n";
  cout << "\n";

  for ( l = 0; l <= 5; l++ )
  {
    for ( i = 0; i < m; i++ )
    {
      ind[i] = l;
    }

    nd = lagrange_interp_nd_size2 ( m, ind );
    xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
    zd = f_sinr ( m, nd, xd );
//
//  Evaluate.
//
    zi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );

    e = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

    cout << "  " << setw(5) << l
         << "  " << setw(5) << nd
         << "  " << setw(10) << e << "\n";

    delete [] xd;
    delete [] zd;
    delete [] zi;
  }

  delete [] a;
  delete [] b;
  delete [] ind;
  delete [] xi;
  delete [] ze;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 interpolates in 3D, using anisotropic resolution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double e;
  int i;
  int j;
  int l;
  int m;
  int *n_1d;
  int nd;
  int ni;
  int order;
  int seed;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;

  cout << "\n";
  cout << "TEST09:\n";
  cout << "  Interpolate in 3D, using orders.\n";
  cout << "  Use a sequence of increasing orders.\n";
  cout << "  Use anisotropic resolution.\n";
  cout << "  The interpoland is a polynomial of degrees 3, 5, 2\n";
  cout << "  so our orders need to be at least 4, 6, 3 to match it.\n";

  m = 3;

  n_1d = new int[m];

  a = new double[m];
  b = new double[m];
  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }

  ni = 20;
  seed = 123456789;
  xi = r8mat_uniform_01_new ( m, ni, seed );
  ze = f_poly352 ( m, ni, xi );

  cout << "\n";
  cout << "  Level     Orders   Average Error\n";
  cout << "\n";

  for ( l = 0; l <= 10; l++ )
  {
    if ( l == 0 )
    {
      n_1d[0] = 1;
      n_1d[1] = 1;
      n_1d[2] = 1;
    }
    else if ( l == 1 )
    {
      n_1d[0] = 2;
      n_1d[1] = 1;
      n_1d[2] = 1;
    }
    else if ( l == 2 )
    {
      n_1d[0] = 1;
      n_1d[1] = 2;
      n_1d[2] = 1;
    }
    else if ( l == 3 )
    {
      n_1d[0] = 1;
      n_1d[1] = 1;
      n_1d[2] = 2;
    }
    else if ( l == 4 )
    {
      n_1d[0] = 4;
      n_1d[1] = 2;
      n_1d[2] = 2;
    }
    else if ( l == 5 )
    {
      n_1d[0] = 2;
      n_1d[1] = 4;
      n_1d[2] = 2;
    }
    else if ( l == 6 )
    {
      n_1d[0] = 2;
      n_1d[1] = 2;
      n_1d[2] = 4;
    }
    else if ( l == 8 )
    {
      n_1d[0] = 6;
      n_1d[1] = 4;
      n_1d[2] = 4;
    }
    else if ( l == 9 )
    {
      n_1d[0] = 4;
      n_1d[1] = 6;
      n_1d[2] = 4;
    }
    else if ( l == 10 )
    {
      n_1d[0] = 4;
      n_1d[1] = 4;
      n_1d[2] = 6;
    }
  
    nd = lagrange_interp_nd_size ( m, n_1d );
    xd = lagrange_interp_nd_grid ( m, n_1d, a, b, nd );
    zd = f_poly352 ( m, nd, xd );
//
//  Evaluate.
//
    zi = lagrange_interp_nd_value ( m, n_1d, a, b, nd, zd, ni, xi );

    e = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

    cout << "  " << setw(5) << l
         << "  " << setw(5) << n_1d[0]
         << "  " << setw(5) << n_1d[1]
         << "  " << setw(5) << n_1d[2]
         << "  " << setw(10) << e << "\n";

    delete [] xd;
    delete [] zd;
    delete [] zi;
  }

  delete [] a;
  delete [] b;
  delete [] n_1d;
  delete [] xi;
  delete [] ze;

  return;
}
//****************************************************************************80

double *f_sinr ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the points.
//
//    Output, double F_SINR[N], the value of the function at each point.
//
{
  int i;
  int j;
  double r;
  double *z;

  z = new double[n];

  for ( j = 0; j < n; j++ )
  {
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r + pow ( x[i+j*m], 2 );
    }
    r = sqrt ( r );
    z[j] = sin ( r );
  }
  return z;
}
//****************************************************************************80

double *f_poly352 ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F_POLY253 is a scalar function of a 3-dimensional argument, to be interpolated.
//
//  Discussion:
//
//    The polynomial is of maximum degrees 3, 5, and 2, in the three variables.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the points.
//
//    Output, double F_POLY253[N], the value of the function at each point.
//
{
  int j;
  double *z;

  z = new double[n];

  for ( j = 0; j < n; j++ )
  {
    z[j] = 
        1.0 
      + pow ( x[0+j*m], 2 ) * pow ( x[1+j*m], 5 ) * pow ( x[2+j*m], 2 ) 
      +       x[0+j*m]      * pow ( x[1+j*m], 2 ) *       x[2+j*m]
      + pow ( x[0+j*m], 3 );
  }

  return z;
}
