# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "bernstein_polynomial.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BERNSTEIN_POLYNOMIAL_PRB.
//
//  Discussion:
//
//    BERNSTEIN_POLYNOMIAL_PRB tests the BERNSTEIN_POLYNOMIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "BERNSTEIN_POLYNOMIAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BERNSTEIN_POLYNOMIAL library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BERNSTEIN_POLYNOMIAL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << " \n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BERNSTEIN_POLY_01 and BERNSTEIN_POLY_01_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double *bvec;
  int k;
  int n;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BERNSTEIN_POLY_01 evaluates the Bernstein polynomials\n";
  cout << "  based on the interval [0,1].\n";
  cout << "  BERNSTEIN_POLY_01_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "     N     K     X       Exact         BP01(N,K)(X)\n";
  cout << "\n";

  n_data = 0;

  while ( true )
  {
    bernstein_poly_01_values ( &n_data, &n, &k, &x, &b );

    if ( n_data == 0 )
    {
      break;
    }

    bvec = bernstein_poly_01 ( n, x );

    cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << x
         << "  " << setw(14) << b
         << "  " << setw(14) << bvec[k] << "\n";

    delete [] bvec;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BERNSTEIN_POLY_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *bern;
  int k;
  int n = 10;
  double x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  BERNSTEIN_POLY_AB evaluates Bernstein polynomials over an\n";
  cout << "  arbitrary interval [A,B].\n";
  cout << "\n";
  cout << "  Here, we demonstrate that \n";
  cout << "    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)\n";
  cout << "  provided only that\n";
  cout << "    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).\n";

  x = 0.3;
  a = 0.0;
  b = 1.0;
  bern = bernstein_poly_ab ( n, a, b, x );
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n";
  for ( k = 0; k <= n; k++ )
  {
    cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;
 
  x = 1.3;
  a = 1.0;
  b = 2.0;
  bern = bernstein_poly_ab ( n, a, b, x );
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n"; 
  for ( k = 0; k <= n; k++ )
  {
    cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;

  x = 2.6;
  a = 2.0;
  b = 4.0;
  bern = bernstein_poly_ab ( n, a, b, x );
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n";
 
  for ( k = 0; k <= n; k++ )
  {
    cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests the Partition-of-Unity property.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *bvec;
  int n;
  int n_data;
  int seed;
  double x;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  BERNSTEIN_POLY_01 evaluates the Bernstein polynomials\n";
  cout << "  based on the interval [0,1].\n";
  cout << "\n";
  cout << "  Here we test the partition of unity property.\n";
  cout << "\n";
  cout << "     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)\n";
  cout << "\n";

  seed = 123456789;

  for ( n = 0; n <= 10; n++ )
  {
    x = r8_uniform_01 ( &seed );

    bvec = bernstein_poly_01 ( n, x );

    cout << "  " << setw(4) << n
         << "  " << setw(7) << x
         << "  " << setw(14) << r8vec_sum ( n + 1, bvec ) << "\n";

    delete [] bvec;
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BERNSTEIN_POLY_AB_APPROX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error_max;
  int i;
  int maxdata = 20;
  int ndata;
  int nsample;
  int nval = 501;
  double *xdata;
  double *xval;
  double *ydata;
  double *yval;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  BERNSTEIN_POLY_AB_APPROX evaluates the Bernstein polynomial\n";
  cout << "  approximant to a function F(X).\n";

  a = 1.0;
  b = 3.0;

  cout << "\n";
  cout << "     N      Max Error\n";
  cout << "\n";

  for ( ndata = 0; ndata <= maxdata; ndata++ )
  {
//
//  Generate data values.
//
    xdata = new double[ndata+1];
    ydata = new double[ndata+1];
    for ( i = 0; i <= ndata; i++)
    {
      if ( ndata == 0 )
      {
        xdata[i] = 0.5 * ( a + b );
      }
      else
      {
        xdata[i] = ( ( double ) ( ndata - i ) * a   
                   + ( double ) (         i ) * b ) 
                   / ( double ) ( ndata     );
      }
      ydata[i] = sin ( xdata[i] );
    }
//
//  Compare the true function and the approximant.
//
    xval = r8vec_linspace_new ( nval, a, b );

    error_max = 0.0;

    yval = bernstein_poly_ab_approx ( ndata, a, b, ydata, nval, xval );

    error_max = 0.0;
    for ( i = 0; i < nval; i++ )
    {
      error_max = r8_max ( error_max, fabs ( yval[i] - sin ( xval[i] ) ) );
    }
    cout << "  " << setw(4) << ndata
         << "  " << setw(14) << error_max << "\n";

    delete [] xdata;
    delete [] xval;
    delete [] ydata;
    delete [] yval;
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests BERNSTEIN_MATRIX and BERNSTEIN_MATRIX_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double a_norm_frobenius;
  double *b;
  double b_norm_frobenius;
  double *c;
  double error_norm_frobenius;
  int n;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  BERNSTEIN_MATRIX returns a matrix A which transforms a\n";
  cout << "  polynomial coefficient vector from the power basis to\n";
  cout << "  the Bernstein basis.\n";
  cout << "  BERNSTEIN_MATRIX_INVERSE computes the inverse B.\n";
  cout << "\n";
  cout << "     N     ||A||            ||B||      ||I-A*B||\n";
  cout << "\n";

  for ( n = 5; n <= 15; n++ )
  {
    a = bernstein_matrix ( n );
    a_norm_frobenius = r8mat_norm_fro ( n, n, a );

    b = bernstein_matrix_inverse ( n );
    b_norm_frobenius = r8mat_norm_fro ( n, n, b );

    c = r8mat_mm_new ( n, n, n, a, b );
    error_norm_frobenius = r8mat_is_identity ( n, c );

    cout << "  " << setw(4) << n
         << "  " << setw(14) << a_norm_frobenius
         << "  " << setw(14) << b_norm_frobenius
         << "  " << setw(14) << error_norm_frobenius << "\n";

    delete [] a;
    delete [] b;
    delete [] c;
  }
  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 uses BERNSTEIN_MATRIX_INVERSE to describe Bernstein polynomials.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2014
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *ax;
  int i;
  int k;
  int n;
  double *x;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  BERNSTEIN_MATRIX returns a matrix which\n";
  cout << "  transforms a polynomial coefficient vector\n";
  cout << "  from the the Bernstein basis to the power basis.\n";
  cout << "  We can use this to get explicit values of the\n";
  cout << "  4-th degree Bernstein polynomial coefficients as\n";
  cout << "\n";
  cout << "    B(4,K)(X) = C4 * x^4\n";
  cout << "              + C3 * x^3\n";
  cout << "              + C2 * x^2\n";
  cout << "              + C1 * x\n";
  cout << "              + C0 * 1\n";

  n = 5;
  cout << "\n";
  cout << "     K       C4           C3            C2";
  cout << "            C1             C0\n";
  cout << "\n";

  a = bernstein_matrix ( n );
  x = new double[n];

  for ( k = 0; k < n; k++ )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = 0.0;
    }
    x[k] = 1.0;

    ax = r8mat_mv_new ( n, n, a, x );

    cout << "  " << setw(4) << k << "  ";
    for ( i = 0; i < n; i++ )
    {
      cout << "  " << setw(14) << ax[i];
    }
    cout << "\n";
  }

  delete [] a;
  delete [] ax;
  delete [] x;

  return;
}

