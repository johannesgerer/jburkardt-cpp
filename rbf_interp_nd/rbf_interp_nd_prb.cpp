# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "rbf_interp_nd.hpp"
# include "r8lib.hpp"

int main ( );
void rbf_interp_nd_test01 ( );
void rbf_interp_nd_test02 ( );
void rbf_interp_nd_test03 ( );
void rbf_interp_nd_test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_ND_TEST tests RBF_INTERP_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "RBF_INTERP_ND_TEST:\n";
  cout << "  C++ version\n";
  cout << "  Test the RBF_INTERP_ND library.\n";
  cout << "  The R8LIB library is also needed.\n";

  rbf_interp_nd_test01 ( );
  rbf_interp_nd_test02 ( );
  rbf_interp_nd_test03 ( );
  rbf_interp_nd_test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "RBF_INTERP_ND_TEST:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void rbf_interp_nd_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_ND_TEST01 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double app_error;
  double b;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  cout << "\n";
  cout << "RBF_INTERP_ND_TEST01:\n";
  cout << "  RBF_WEIGHT computes weights for RBF interpolation.\n";
  cout << "  RBF_INTERP_ND evaluates the RBF interpolant.\n";
  cout << "  Use the multiquadratic basis function PHI1(R).\n";

  a = 0.0;
  b = 2.0;

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = new double[m*nd];

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  r8mat_transpose_print ( m, nd, xd, "  The product points:" );

  r0 = ( b - a ) / ( double ) ( n1d );

  cout << "\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  fd = new double[nd];
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  r8vec_print ( nd, fd, "  Function data:" );

  w = rbf_weight ( m, nd, xd, r0, phi1, fd );

  r8vec_print ( nd, w, "  Weight vector:" );
//
//  #1: Interpolation test.  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] fi;
  delete [] xi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi1, w, ni, xi );

  fe = new double[ni];
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per 1000 samples = " << app_error << "\n";

  delete [] fd;
  delete [] fe;
  delete [] fi;
  delete [] w;
  delete [] x1d;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80

void rbf_interp_nd_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_ND_TEST02 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double app_error;
  double b;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  cout << "\n";
  cout << "RBF_INTERP_ND_TEST02:\n";
  cout << "  RBF_WEIGHT computes weights for RBF interpolation.\n";
  cout << "  RBF_INTERP_ND evaluates the RBF interpolant.\n";
  cout << "  Use the inverse multiquadratic basis function PHI2(R).\n";

  a = 0.0;
  b = 2.0;

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = new double[m*nd];

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  r8mat_transpose_print ( m, nd, xd, "  The product points:" );

  r0 = ( b - a ) / ( double ) ( n1d );

  cout << "\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  fd = new double[nd];
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  r8vec_print ( nd, fd, "  Function data:" );

  w = rbf_weight ( m, nd, xd, r0, phi2, fd );

  r8vec_print ( nd, w, "  Weight vector:" );
//
//  #1: Interpolation test.  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] fi;
  delete [] xi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi2, w, ni, xi );

  fe = new double[ni];
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per 1000 samples = " << app_error << "\n";

  delete [] fd;
  delete [] fe;
  delete [] fi;
  delete [] w;
  delete [] x1d;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80

void rbf_interp_nd_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_ND_TEST03 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double app_error;
  double b;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  cout << "\n";
  cout << "RBF_INTERP_ND_TEST03:\n";
  cout << "  RBF_WEIGHT computes weights for RBF interpolation.\n";
  cout << "  RBF_INTERP_ND evaluates the RBF interpolant.\n";
  cout << "  Use the thin-plate spline basis function PHI3(R).\n";

  a = 0.0;
  b = 2.0;

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = new double[m*nd];

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  r8mat_transpose_print ( m, nd, xd, "  The product points:" );

  r0 = ( b - a ) / ( double ) ( n1d );

  cout << "\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  fd = new double[nd];
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  r8vec_print ( nd, fd, "  Function data:" );

  w = rbf_weight ( m, nd, xd, r0, phi3, fd );

  r8vec_print ( nd, w, "  Weight vector:" );
//
//  #1: Interpolation test.  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] fi;
  delete [] xi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi3, w, ni, xi );

  fe = new double[ni];
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per 1000 samples = " << app_error << "\n";

  delete [] fd;
  delete [] fe;
  delete [] fi;
  delete [] w;
  delete [] x1d;
  delete [] xd;
  delete [] xi;

  return;
}
//****************************************************************************80

void rbf_interp_nd_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    RBF_INTERP_ND_TEST04 tests RBF_WEIGHTS and RBF_INTERP_ND with PHI4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double app_error;
  double b;
  double *fd;
  double *fe;
  double *fi;
  int i;
  double int_error;
  int j;
  int m = 2;
  int n1d = 5;
  int nd;
  int ni;
  double r0;
  int seed;
  double *w;
  double *x1d;
  double *xd;
  double *xi;

  cout << "\n";
  cout << "RBF_INTERP_ND_TEST04:\n";
  cout << "  RBF_WEIGHT computes weights for RBF interpolation.\n";
  cout << "  RBF_INTERP_ND evaluates the RBF interpolant.\n";
  cout << "  Use the gaussian basis function PHI4(R).\n";

  a = 0.0;
  b = 2.0;

  x1d = r8vec_linspace_new ( n1d, a, b );
  nd = i4_power ( n1d, m );
  xd = new double[m*nd];

  for ( i = 0; i < m; i++ )
  {
    r8vec_direct_product ( i, n1d, x1d, m, nd, xd );
  }

  r8mat_transpose_print ( m, nd, xd, "  The product points:" );

  r0 = ( b - a ) / ( double ) ( n1d );

  cout << "\n";
  cout << "  Scale factor R0 = " << r0 << "\n";

  fd = new double[nd];
  for ( j = 0; j < nd; j++ )
  {
    fd[j] = xd[0+j*m] * xd[1+j*m] * exp ( - xd[0+j*m] * xd[1+j*m] );
  }
  r8vec_print ( nd, fd, "  Function data:" );

  w = rbf_weight ( m, nd, xd, r0, phi4, fd );

  r8vec_print ( nd, w, "  Weight vector:" );
//
//  #1: Interpolation test.  Does interpolant match function at interpolation points?
//
  ni = nd;
  xi = r8mat_copy_new ( m, ni, xd );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi );

  int_error = r8vec_norm_affine ( nd, fd, fi ) / ( double ) ( nd );

  cout << "\n";
  cout << "  L2 interpolation error averaged per interpolant node = " << int_error << "\n";

  delete [] fi;
  delete [] xi;
//
//  #2: Approximation test.  Estimate the integral (f-interp(f))^2.
//
  ni = 1000;
  seed = 123456789;

  xi = r8mat_uniform_ab_new ( m, ni, a, b, seed );

  fi = rbf_interp_nd ( m, nd, xd, r0, phi4, w, ni, xi );

  fe = new double[ni];
  for ( j = 0; j < ni; j++ )
  {
    fe[j] = xi[0+j*m] * xi[1+j*m] * exp ( - xi[0+j*m] * xi[1+j*m] );
  }

  app_error = pow ( b - a, m ) * r8vec_norm_affine ( ni, fi, fe ) / ( double ) ( ni );

  cout << "\n";
  cout << "  L2 approximation error averaged per 1000 samples = " << app_error << "\n";

  delete [] fd;
  delete [] fe;
  delete [] fi;
  delete [] w;
  delete [] x1d;
  delete [] xd;
  delete [] xi;

  return;
}
