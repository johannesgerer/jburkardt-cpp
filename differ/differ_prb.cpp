# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "differ.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DIFFER_PRB.
//
//  Discussion:
//
//    DIFFER_PRB tests the DIFFER library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "DIFFER_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the DIFFER library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "DIFFER_PRB:\n";
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
//    TEST01 tests DIFFER_MATRIX.
//
//  Discussion:
//
//    DIFFER_MATRIX computes a modified Vandermonde matrix A1.
//
//    The solution of a system A1 * X1 = B is related to the solution
//    of the system A2 * X2 = B, where A2 is the standard Vandermonde
//    matrix, simply by X2(I) = X1(I) * A(I,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int i;
  int info;
  int job;
  int n = 4;
  double stencil[4] = { 2.5, 3.3, -1.3, 0.5 };
  double x[4] = { 1.0, 2.0, 3.0, 4.0 };
  double *x1;
  double *x2;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Demonstrate that the DIFFER matrix is 'really'\n";
  cout << "  a Vandermonde matrix.\n";

  a = differ_matrix ( n, stencil );
  r8mat_print ( n, n, a, "  Stencil matrix:" );
  b = r8mat_mv_new ( n, n, a, x );
//
//  Set up and solve the DIFFER system.
//
  a = differ_matrix ( n, stencil );
  x1 = r8mat_fs_new ( n, a, b );

  r8vec_print ( n, x1, "  Solution of DIFFER system:" );
//
//  R8VM_SL solves the related Vandermonde system.
//  A simple transformation gives us the solution to the DIFFER system.
//
  job = 0;
  x2 = r8vm_sl_new ( n, stencil, b, job, info );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "TEST01 - Warning!\n";
    cerr << "  VANDERMONDE system is singular.\n";
    exit ( 1 );
  }

  r8vec_print ( n, x2, "  Solution of VANDERMONDE system:" );

  for ( i = 0; i < n; i++ )
  {
    x2[i] = x2[i] / stencil[i];
  }
  r8vec_print ( n, x2, "  Transformed solution of VANDERMONDE system:" );

  delete [] a;
  delete [] b;
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
//    TEST02 tests DIFFER_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double err;
  int n;
  int n_max = 8;
  int seed;
  int test;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  DIFFER_INVERSE returns the inverse of a DIFFER matrix;\n";
  cout << "\n";
  cout << "   N    Inverse error\n";

  seed = 123456789;

  for ( n = 2; n <= n_max; n++ )
  {
    cout << "\n";

    for ( test = 1; test <= 5; test++ )
    {
      x = r8vec_uniform_01_new ( n, seed );
      a = differ_matrix ( n, x );
      b = differ_inverse ( n, x );
      err = inverse_error ( n, a, b );
      cout << "  " << setw(2)  << n
           << "  " << setw(14) <<  err << "\n";
      delete [] a;
      delete [] b;
      delete [] x;
    }
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DIFFER_MATRIX.
//
//  Discussion:
//
//    Reproduce a specific example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  double *c;
  double df;
  double dfdx;
  double dx;
  int i;
  int n = 4;
  int order;
  double stencil[4] = { -3.0, -2.0, -1.0, 1.0 };
  double x0;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Reproduce a specific example.\n";
//
//  Compute the coefficients for a specific stencil.
//
  b = new double[n];
  for ( i = 0; i < n; i++ )
  {
    b[i] = 0.0;
  }
  order = 1;
  b[order-1] = 1.0;
  a = differ_matrix ( n, stencil );
  c = r8mat_fs_new ( n, a, b );

  r8vec_print ( n, c, "  Solution of DIFFER system:" );
//
//  Use the coefficients C to estimate the first derivative of EXP(X)
//  at X0, using a spacing of DX = 0.1.
//
  x0 = 1.3;
  dx = 0.1;
  df = 0.0;
  for ( i = 0; i < n; i++ )
  {
    df = df + c[i] * ( exp ( x0 + stencil[i] * dx ) - exp ( x0 ) );
  }
  dfdx = df / dx;

  cout << "\n";
  cout << "  DFDX =         " << dfdx << "\n";
  cout << "  d exp(x) /dx = " << exp ( x0 ) << "\n";

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests DIFFER_FORWARD, DIFFER_BACKWARD, DIFFER_CENTRAL.
//
//  Discussion:
//
//    Evaluate the coefficients for uniformly spaced finite difference
//    approximations of derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double h;
  string label;
  int n;
  int o;
  int p;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  DIFFER_FORWARD,\n";
  cout << "  DIFFER_BACKWARD, and\n";
  cout << "  DIFFER_CENTRAL produce coefficients for difference\n";
  cout << "  approximations of the O-th derivative,\n";
  cout << "  with error of order H^P, for a uniform spacing of H.\n";

  h = 1.0;
  cout << "\n";
  cout << "  Use a spacing of H = " << h << " for all examples.\n";
//
//  Forward difference approximation to the third derivative with error of O(h).
//
  o = 3;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_forward ( h, o, p, c, x );
  label = "  Forward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Backward difference approximation to the third derivative with error of O(h).
//
  o = 3;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_backward ( h, o, p, c, x );
  label = "  Backward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Central difference approximation to the third derivative with error of O(h^2).
//
  o = 3;
  p = 2;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_central ( h, o, p, c, x );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Central difference approximation to the third derivative with error of O(h^4).
//
  o = 3;
  p = 4;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_central ( h, o, p, c, x );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Forward difference approximation to the fourth derivative with error of O(h).
//
  o = 4;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_forward ( h, o, p, c, x );
  label = "  Forward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Backward difference approximation to the fourth derivative with error of O(h).
//
  o = 4;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_backward ( h, o, p, c, x );
  label = "  Backward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//   Central difference approximation to the fourth derivative with error of O(h^3).
//
  o = 4;
  p = 3;
  n = o + p;
  c = new double[n];
  x = new double[n];
  differ_central ( h, o, p, c, x );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests DIFFER_STENCIL.
//
//  Discussion:
//
//    Evaluate the coefficients for uniformly spaced finite difference
//    approximations of derivatives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double h;
  int i;
  string label;
  int n;
  int o;
  int p;
  double *x;
  double x0;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  DIFFER_STENCIL produces coefficients for difference\n";
  cout << "  approximations of the O-th derivative,\n";
  cout << "  using arbitrarily spaced data, with maximum spacing H\n";
  cout << "  with error of order H^P.\n";
//
//  Let X0 = 1.0.
//
  x0 = 0.0;
  h = 1.0;
  cout << "\n";
  cout << "  Use a spacing of H = " << h << " for all examples.\n";
//
//  Forward difference approximation to the third derivative with error of O(h).
//
  o = 3;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Forward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Backward difference approximation to the third derivative with error of O(h).
//
  o = 3;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 - n ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Backward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Central difference approximation to the third derivative with error of O(h^2).
//
  o = 3;
  p = 2;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Central difference approximation to the third derivative with error of O(h^4).
//
  o = 3;
  p = 4;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Forward difference approximation to the fourth derivative with error of O(h).
//
  o = 4;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Forward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//  Backward difference approximation to the fourth derivative with error of O(h).
//
  o = 4;
  p = 1;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( i + 1 - n ) * h;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Backward difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;
//
//   Central difference approximation to the fourth derivative with error of O(h^3).
//
  o = 4;
  p = 3;
  n = o + p;
  c = new double[n];
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) * h / 2.0;
  }
  differ_stencil ( x0, o, p, x, c );
  label = "  Central difference coefficients, O = " + i4_to_string ( o )
    + ", P = " + i4_to_string ( p );
  r8vec2_print ( n, x, c, label );
  delete [] c;
  delete [] x;

  return;
}
