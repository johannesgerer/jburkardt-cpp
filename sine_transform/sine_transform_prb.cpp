# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "sine_transform.hpp"

int main ( );
void sine_transform_test01 ( );
void sine_transform_test02 ( );
void sine_transform_test03 ( );
void sine_transform_test04 ( );
void sine_transform_test05 ( );
double cosine_sum ( double x );
double poly5 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SINE_TRANSFORM_PRB.
//
//  Discussion:
//
//    SINE_TRANSFORM_PRB tests SINE_TRANSFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SINE_TRANSFORM_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the SINE_TRANSFORM library.\n";

  sine_transform_test01 ( );
  sine_transform_test02 ( );
  sine_transform_test03 ( );
  sine_transform_test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SINE_TRANSFORM_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void sine_transform_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_TEST01 demonstrates that the transform is its own inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int n = 10;
  int seed;
  double *r;
  double *s;
  double *t;

  seed = 123456789;

  cout << "\n";
  cout << "SINE_TRANSFORM_TEST01:\n";
  cout << "  SINE_TRANSFORM_DATA does a sine transform of data\n";
  cout << "  defined by a vector.\n";
  cout << "\n";
  cout << "  Demonstrate that the transform is its own inverse.\n";
  cout << "  Let R be a random N vector.\n";
  cout << "  Let S be the transform of D.\n";
  cout << "  Let T be the transform of E.\n";
  cout << "  Then R and T will be equal.\n";

  r = r8vec_uniform_01_new ( n, &seed );
  s = sine_transform_data ( n, r );
  t = sine_transform_data ( n, s );

  cout << "\n";
  cout << "     I      R(I)        S(I)        T(I)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << r[i]
         << "  " << setw(10) << s[i]
         << "  " << setw(10) << t[i] << "\n";
  }

  delete [] r;
  delete [] s;
  delete [] t;

  return;
}
//****************************************************************************80

void sine_transform_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_TEST02 uses the functional form of the routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *f1;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 9;
  double *s;
  double *x;

  a = 1.0;
  b = 3.0;
//
//  Evenly spaced points between A and B, but omitting
//  A and B themselves.
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }

  cout << "\n";
  cout << "SINE_TRANSFORM_TEST02:\n";
  cout << "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n";
  cout << "  defined by a function F(X) evaluated at equally spaced\n";
  cout << "  points in an interval [A,B].\n";
  cout << "\n";
  cout << "  Demonstrate that the transform is its own inverse.\n";
  cout << "  Let X(0:N+1) be N+2 equally spaced points in [A,B].\n";
  cout << "  Let S be the transform of F(X(1:N)).\n";
  cout << "  Let F1 be the linear interpolant of (A,F(A)), (B,F(B)).\n";
  cout << "  Let F2 be the transform of S.\n";
  cout << "  Then F(X(1:N)) = F1(X(1:N)) + F2(1:N).\n";

  s = sine_transform_function ( n, a, b, poly5 );

  fa = poly5 ( a );
  fb = poly5 ( b );

  f1 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f1[i] = ( ( b - x[i]     ) * fa   
            + (     x[i] - a ) * fb ) 
            / ( b        - a );
  }

  f2 = sine_transform_data ( n, s );

  cout << "\n";
  cout << "     I      X(I)      F(X(I))       S           F1          F2          F1+F2\n";
  cout << "\n";
  cout << "  " << setw(4) << 0
       << "  " << setw(10) << a
       << "  " << setw(10) << poly5 ( a )
       << "  " << setw(10) << 0.0
       << "  " << setw(10) << fa
       << "  " << setw(10) << 0.0
       << "  " << setw(10) << fa << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i + 1
         << "  " << setw(10) << x[i]
         << "  " << setw(10) << poly5 ( x[i] )
         << "  " << setw(10) << s[i]
         << "  " << setw(10) << f1[i]
         << "  " << setw(10) << f2[i]
         << "  " << setw(10) << f1[i] + f2[i] << "\n";
  }

  cout << "  " << setw(4) << n + 1
       << "  " << setw(10) << b
       << "  " << setw(10) << poly5 ( b )
       << "  " << setw(10) << 0.0
       << "  " << setw(10) << fb
       << "  " << setw(10) << 0.0
       << "  " << setw(10) << fb << "\n";

  delete [] f1;
  delete [] f2;
  delete [] s;
  delete [] x;

  return;
}
//****************************************************************************80

void sine_transform_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_TEST03 evaluates the sine transform interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 9;
  int n2 = 1 + 2 * ( n + 1 );
  double *s;
  double *x;
  double *x2;

  cout << "\n";
  cout << "SINE_TRANSFORM_TEST03:\n";
  cout << "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n";
  cout << "  defined by a function F(X) evaluated at N equally spaced\n";
  cout << "  points in an interval [A,B].\n";
  cout << "  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.\n";
  cout << "\n";
  cout << "  The interpolant will be 0 at the 0th and (N+1)-th points.\n";
  cout << "  It equals the function at points 1 through N.\n";
  cout << "  In between, it can approximate smooth functions,\n";
  cout << "  and the approximation improves with N.\n";
//
//  N determines the number of data points, indexed by 1 to N.  
//  However, we essentially have N+2 data points, indexed 0 to N+1,
//  with the data value being 0 at the first and last auxilliary points.
//
  a = 1.0;
  b = 4.0;
//
//  Evenly spaced points between A and B, but omitting
//  A and B themselves.
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }
//
//  Determine the interpolant coefficients.
//
  s = sine_transform_function ( n, a, b, poly5 );

  cout << "\n";
  cout << "     I      X(I)      F(X(I))        S(I)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << x[i]
         << "  " << setw(10) << poly5 ( x[i] )
         << "  " << setw(10) << s[i] << "\n";
  }
//
//  Evaluate the interpolant.
//
  fa = poly5 ( a );
  fb = poly5 ( b );
//
//  Evenly spaced points between A and B, including A and B,
//  and twice the density of the previous set of points.
//
  x2 = new double[n2];

  for ( i = 0; i < n2; i++ )
  {
    x2[i] = ( ( double ) ( n2 - i - 1 ) * a   
            + ( double ) (      i     ) * b ) 
            / ( double ) ( n2     - 1 );
  }

  f2 = sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2 );

  cout << "\n";
  cout << "     I      X            F(X)        FHAT(X)\n";
  cout << "\n";

  for ( i = 0; i < n2; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << x2[i]
         << "  " << setw(10) << poly5 ( x2[i] )
         << "  " << setw(10) << f2[i] << "\n";
  }

  delete [] f2;
  delete [] s;
  delete [] x;
  delete [] x2;

  return;
}
//****************************************************************************80

void sine_transform_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TRANSFORM_TEST04 evaluates the sine transform interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double *f2;
  double fa;
  double fb;
  int i;
  int n = 15;
  int n2 = 1 + 5 * ( n + 1 );
  double *s;
  double *x;
  double *x2;

  cout << "\n";
  cout << "SINE_TRANSFORM_TEST04:\n";
  cout << "  SINE_TRANSFORM_FUNCTION does a sine transform of data\n";
  cout << "  defined by a function F(X) evaluated at N equally spaced\n";
  cout << "  points in an interval [A,B].\n";
  cout << "  SINE_TRANSFORM_INTERPOLANT evaluates the interpolant.\n";
  cout << "\n";
  cout << "  The interpolant will be 0 at the 0th and (N+1)-th points.\n";
  cout << "  It equals the function at points 1 through N.\n";
  cout << "  In between, it can approximate smooth functions,\n";
  cout << "  and the approximation improves with N.\n";
//
//  N determines the number of data points, indexed by 1 to N.  
//  However, we essentially have N+2 data points, indexed 0 to N+1,
//  with the data value being 0 at the first and last auxilliary points.
//
  a = 0.0;
  b = 7.0;
//
//  Evenly spaced points between A and B, but omitting
//  A and B themselves.
//
  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * a   
           + ( double ) (     i + 1 ) * b ) 
           / ( double ) ( n     + 1 );
  }
//
//  Determine the interpolant coefficients.
//
  s = sine_transform_function ( n, a, b, cosine_sum );
//
//  Evaluate the interpolant.
//
  fa = cosine_sum ( a );
  fb = cosine_sum ( b );
//
//  Evenly spaced points between A and B, including A and B,
//  and twice the density of the previous set of points.
//
  x2 = new double[n2];

  for ( i = 0; i < n2; i++ )
  {
    x2[i] = ( ( double ) ( n2 - i - 1 ) * a   
            + ( double ) (      i     ) * b ) 
            / ( double ) ( n2     - 1 );
  }

  f2 = sine_transform_interpolant ( n, a, b, fa, fb, s, n2, x2 );

  cout << "\n";
  cout << "  Expect exact agreement every 5th sample.\n";
  cout << "\n";

  for ( i = 0; i < n2; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(10) << x2[i]
         << "  " << setw(10) << cosine_sum ( x2[i] )
         << "  " << setw(10) << f2[i] << "\n";
  }

  delete [] f2;
  delete [] s;
  delete [] x;
  delete [] x2;

  return;
}
//****************************************************************************80

double cosine_sum ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    COSINE_SUM evaluates a function which is a cosine sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double COSINE_SUM, the value.
//
{
  double value;

  value =   cos (       x ) 
    + 5.0 * cos ( 1.6 * x ) 
    - 2.0 * cos ( 2.0 * x ) 
    + 5.0 * cos ( 4.5 * x ) 
    + 7.0 * cos ( 9.0 * x );

  return value;
}
//****************************************************************************80

double poly5 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    POLY5 evaluates a particular fifth-degree polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double POLY5, the value of the polynomial at X.
//
{
  double value;

  value = ( x - 0.1 ) * 
          ( x - 0.2 ) * 
          ( x - 0.4 ) * 
          ( x - 2.1 ) * 
          ( x - 3.0 );

  return value;
}
