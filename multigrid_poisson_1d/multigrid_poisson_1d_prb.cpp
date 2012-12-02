# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>>

using namespace std;

# include "multigrid_poisson_1d.hpp"

int main ( );

void test01_mono ( );
void test01_multi ( );
void test02_mono ( );
void test02_multi ( );
double exact1 ( double x );
double force1 ( double x );
double exact2 ( double x );
double force2 ( double x );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MULTIGRID_POISSON_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "MULTIGRID_POISSON_1D:\n";
  cout << "  C++ version\n";
  cout << "  Test the MULTIGRID_POISSON_1D multigrid library.\n";

  test01_mono ( );
  test01_multi ( );
  test02_mono ( );
  test02_multi ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "MULTIGRID_POISSON_1D:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01_mono ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    TEST01_MONO tests MONOGRID_POISSON_1D on test case 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double x;

  cout << "\n";
  cout << "TEST01_MONO\n";
  cout << "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n";
  cout << "  using the Gauss-Seidel method.\n";
  cout << "\n";
  cout << "  -u''(x) = 1, for 0 < x < 1\n";
  cout << "  u(0) = u(1) = 0.\n";
  cout << "  Solution is u(x) = ( -x^2 + x ) / 2\n";

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = new double[n+1];

    cout << "\n";
    cout << "  Mesh index K = " << k << "\n";
    cout << "  Number of intervals N=2^K = " << n << "\n";
    cout << "  Number of nodes = 2^K+1 =   " << n + 1 << "\n";

    monogrid_poisson_1d ( n, force1, exact1, it_num, u );

    cout << "\n";
    cout << "     I        X(I)      U(I)         U Exact(X(I))\n";
    cout << "\n";
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      cout << "  " << setw(4) << i
           << "  " << setw(10) << x
           << "  " << setw(14) << u[i]
           << "  " << setw(14) << exact1 ( x ) << "\n";
    }

    cout << "\n";

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      difmax = r8_max ( difmax, r8_abs ( u[i] - exact1 ( x ) ) );
    } 
    cout << "  Maximum error = " << difmax << "\n";
    cout << "  Number of iterations = " << it_num << "\n";

    delete [] u;
  }
  return;
}
//****************************************************************************80

void test01_multi ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    TEST01_MULTI tests MULTIGRID_POISSON_1D on test case 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double x;

  cout << "\n";
  cout << "TEST01_MULTI\n";
  cout << "  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n";
  cout << "  using the multigrid method.\n";
  cout << "\n";
  cout << "  -u''(x) = 1, for 0 < x < 1\n";
  cout << "  u(0) = u(1) = 0.\n";
  cout << "  Solution is u(x) = ( -x^2 + x ) / 2\n";

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = new double[n+1];

    cout << "\n";
    cout << "  Mesh index K = " << k << "\n";
    cout << "  Number of intervals N=2^K = " << n << "\n";
    cout << "  Number of nodes = 2^K+1 =   " << n + 1 << "\n";

    multigrid_poisson_1d ( n, force1, exact1, it_num, u );

    cout << "\n";
    cout << "     I        X(I)      U(I)         U Exact(X(I))\n";
    cout << "\n";
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      cout << "  " << setw(4) << i
           << "  " << setw(10) << x
           << "  " << setw(14) << u[i]
           << "  " << setw(14) << exact1 ( x ) << "\n";
    }

    cout << "\n";

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      difmax = r8_max ( difmax, r8_abs ( u[i] - exact1 ( x ) ) );
    } 
    cout << "  Maximum error = " << difmax << "\n";
    cout << "  Number of iterations = " << it_num << "\n";

    delete [] u;
  }
  return;
}
//****************************************************************************80

double exact1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT1 evaluates the exact solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT1, the value of the exact solution at X.
//
{
  double value;

  value = 0.5 * ( - x * x + x );

  return value;
}
//****************************************************************************80

double force1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FORCE1 evaluates the forcing function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double FORCE1, the value of the forcing function at X.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

void test02_mono ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    TEST02_MONO tests MONOGRID_POISSON_1D on test case 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double x;

  cout << "\n";
  cout << "TEST02_MONO\n";
  cout << "  MONOGRID_POISSON_1D solves a 1D Poisson BVP\n";
  cout << "  using the Gauss-Seidel method.\n";
  cout << "\n";
  cout << "  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n";
  cout << "  u(0) = u(1) = 0.\n";
  cout << "  Solution is u(x) = x * (x-1) * exp(x)\n";

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = new double[n+1];

    cout << "\n";
    cout << "  Mesh index K = " << k << "\n";
    cout << "  Number of intervals N=2^K = " << n << "\n";
    cout << "  Number of nodes = 2^K+1 =   " << n + 1 << "\n";

    monogrid_poisson_1d ( n, force2, exact2, it_num, u );

    cout << "\n";
    cout << "     I        X(I)      U(I)         U Exact(X(I))\n";
    cout << "\n";
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      cout << "  " << setw(4) << i
           << "  " << setw(10) << x
           << "  " << setw(14) << u[i]
           << "  " << setw(14) << exact2 ( x ) << "\n";
    }

    cout << "\n";

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      difmax = r8_max ( difmax, r8_abs ( u[i] - exact2 ( x ) ) );
    } 
    cout << "  Maximum error = " << difmax << "\n";
    cout << "  Number of iterations = " << it_num << "\n";

    delete [] u;
  }
  return;
}
//****************************************************************************80

void test02_multi ( ) 

//****************************************************************************80
//
//  Purpose:
//
//    TEST02_MULTI tests MULTIGRID_POISSON_1D on test case 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
{
  double difmax;
  int i;
  int it_num;
  int k;
  int n;
  double *u;
  double x;

  cout << "\n";
  cout << "TEST02_MULTI\n";
  cout << "  MULTIGRID_POISSON_1D solves a 1D Poisson BVP\n";
  cout << "  using the multigrid method.\n";
  cout << "\n";
  cout << "  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1\n";
  cout << "  u(0) = u(1) = 0.\n";
  cout << "  Solution is u(x) = x * (x-1) * exp(x)\n";

  for ( k = 5; k <= 5; k++ )
  {
    n = i4_power ( 2, k );

    u = new double[n+1];

    cout << "\n";
    cout << "  Mesh index K = " << k << "\n";
    cout << "  Number of intervals N=2^K = " << n << "\n";
    cout << "  Number of nodes = 2^K+1 =   " << n + 1 << "\n";

    multigrid_poisson_1d ( n, force2, exact2, it_num, u );

    cout << "\n";
    cout << "     I        X(I)      U(I)         U Exact(X(I))\n";
    cout << "\n";
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      cout << "  " << setw(4) << i
           << "  " << setw(10) << x
           << "  " << setw(14) << u[i]
           << "  " << setw(14) << exact2 ( x ) << "\n";
    }

    cout << "\n";

    difmax = 0.0;
    for ( i = 0; i < n + 1; i++ )
    {
      x = ( double ) ( i ) / ( double ) ( n );
      difmax = r8_max ( difmax, r8_abs ( u[i] - exact2 ( x ) ) );
    } 
    cout << "  Maximum error = " << difmax << "\n";
    cout << "  Number of iterations = " << it_num << "\n";

    delete [] u;
  }
  return;
}
//****************************************************************************80

double exact2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT2 evaluates the exact solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double EXACT2, the value of the exact solution at X.
//
{
  double value;

  value = x * ( x - 1.0 ) * exp ( x );

  return value;
}
//****************************************************************************80

double force2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FORCE2 evaluates the forcing function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Hager,
//    Applied Numerical Linear Algebra,
//    Prentice-Hall, 1988,
//    ISBN13: 978-0130412942,
//    LC: QA184.H33.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Output, double FORCE2, the value of the forcing function at X.
//
{
  double value;

  value = - x * ( x + 3.0 ) * exp ( x );

  return value;
}
