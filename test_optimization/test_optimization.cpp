# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_optimization.hpp"

//****************************************************************************80

void p00_ab ( int problem, int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_AB evaluates the limits of the optimization region for any problem.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, int M, the number of variables.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  if ( problem == 1 )
  {
    p01_ab ( m, a, b );
  }
  else if ( problem == 2 )
  {
    p02_ab ( m, a, b );
  }
  else if ( problem == 3 )
  {
    p03_ab ( m, a, b );
  }
  else if ( problem == 4 )
  {
    p04_ab ( m, a, b );
  }
  else if ( problem == 5 )
  {
    p05_ab ( m, a, b );
  }
  else if ( problem == 6 )
  {
    p06_ab ( m, a, b );
  }
  else if ( problem == 7 )
  {
    p07_ab ( m, a, b );
  }
  else if ( problem == 8 )
  {
    p08_ab ( m, a, b );
  }
  else if ( problem == 9 )
  {
    p09_ab ( m, a, b );
  }
  else if ( problem == 10 )
  {
    p10_ab ( m, a, b );
  }
  else if ( problem == 11 )
  {
    p11_ab ( m, a, b );
  }
  else if ( problem == 12 )
  {
    p12_ab ( m, a, b );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_AB - Fatal error!\n";
    cerr << "  Illegal value of PROBLEM = " << problem << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double *p00_compass_search ( int problem, int m, double x0[], double delta_tol, 
  double delta_init, int k_max, double *fx, int *k )

//****************************************************************************80
//
//  Purpose:
//
//    P00_COMPASS_SEARCH carries out a direct search minimization algorithm.
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
//  Reference:
//
//    Tamara Kolda, Robert Michael Lewis, Virginia Torczon,
//    Optimization by Direct Search: New Perspectives on Some Classical 
//    and Modern Methods,
//    SIAM Review,
//    Volume 45, Number 3, 2003, pages 385-482. 
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, int M, the number of variables.
//
//    Input, double X0[M], a starting estimate for the minimizer.
//
//    Input, double DELTA_TOL, the smallest step size that is allowed.
//
//    Input, double DELTA_INIT, the starting stepsize.  
//
//    Input, int K_MAX, the maximum number of steps allowed.
//
//    Output, double COMPASS_SEARCH[M], the estimated minimizer.
//
//    Output, double *FX, the function value at X.
//
//    Output, int *K, the number of steps taken.
//
{
  int decrease;
  double delta;
  double *fxd;
  int i;
  int ii;
  int n = 1;
  double s;
  double *x;
  double *xd;

  *k = 0;
  x = new double[m];
  xd = new double[m];
  r8vec_copy ( m, x0, x );
  fxd = p00_f ( problem, m, n, x );
  *fx = fxd[0];
  delete [] fxd;

  if ( delta_tol <= 0 )
  {
    cerr << "\n";
    cerr << "P00_COMPASS_SEARCH - Fatal error!\n";
    cerr << "  DELTA_TOL <= 0.0.\n";
    cerr << "  DELTA_TOL = " << delta_tol << "\n";
    exit ( 1 );
  }

  if ( delta_init <= delta_tol )
  {
    cerr << "\n";
    cerr << "P00_COMPASS_SEARCH - Fatal error!\n";
    cerr << "  DELTA_INIT < DELTA_TOL.\n";
    cerr << "  DELTA_INIT = " << delta_init << "\n";
    cerr << "  DELTA_TOL = " << delta_tol << "\n";
    exit ( 1 );
  }

  delta = delta_init;

  while ( *k < k_max )
  {
    *k = *k + 1;
//
//  For each coordinate direction I, seek a lower function value
//  by increasing or decreasing X(I) by DELTA.
//
    decrease = 0;
    s = + 1.0;
    i = 0;

    for ( ii = 1; ii <= 2 * m; ii++ )
    {
      r8vec_copy ( m, x, xd );
      xd[i] = xd[i] + s * delta;
      fxd = p00_f ( problem, m, n, xd );
//
//  As soon as a decrease is noticed, accept the new point.
//
      if ( fxd[0] < *fx )
      {
        r8vec_copy ( m, xd, x );
        *fx = fxd[0];
        decrease = 1;
        break;
      }
      delete [] fxd;

      s = - s;
      if ( s == + 1.0 )
      {
        i = i + 1;
      }
    }
//
//  If no decrease occurred, reduce DELTA.
//
    if ( !decrease )
    {
      delta = delta / 2.0;
      if ( delta < delta_tol )
      {
        break;
      }
    }
  }
  delete [] xd;

  return x;
}
//****************************************************************************80

double *p00_f ( int problem, int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F evaluates the objective function for any problem.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int PROBLEM, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the argument of the objective function.
//
//    Output, double P00_F[N], the objective function evaluated at
//    each argument.
//
{
  double *f;

  if ( problem == 1 )
  {
    f = p01_f ( m, n, x );
  }
  else if ( problem == 2 )
  {
    f = p02_f ( m, n, x );
  }
  else if ( problem == 3 )
  {
    f = p03_f ( m, n, x );
  }
  else if ( problem == 4 )
  {
    f = p04_f ( m, n, x );
  }
  else if ( problem == 5 )
  {
    f = p05_f ( m, n, x );
  }
  else if ( problem == 6 )
  {
    f = p06_f ( m, n, x );
  }
  else if ( problem == 7 )
  {
    f = p07_f ( m, n, x );
  }
  else if ( problem == 8 )
  {
    f = p08_f ( m, n, x );
  }
  else if ( problem == 9 )
  {
    f = p09_f ( m, n, x );
  }
  else if ( problem == 10 )
  {
    f = p10_f ( m, n, x );
  }
  else if ( problem == 11 )
  {
    f = p11_f ( m, n, x );
  }
  else if ( problem == 12 )
  {
    f = p12_f ( m, n, x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_F - Fatal error!\n";
    cerr << "  Illegal value of PROBLEM = " << problem << "\n";
    exit ( 1 );
  }
  return f;
}
//****************************************************************************80

int p00_problem_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROBLEM_NUM returns the number of problems available.
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
//   Output, int P00_PROBLEM_NUM, the number of problems available.
//
{
  int problem_num;

  problem_num = 12;

  return problem_num;
}
//****************************************************************************80

double *p00_sol ( int problem, int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P00_SOL returns the solution for any problem.
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
//    Input, int PROBLEM, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, P00_SOL[M], the solution, if known.
//
{
  double *x;

  if ( problem == 1 )
  {
    x = p01_sol ( m, know );
  }
  else if ( problem == 2 )
  {
    x = p02_sol ( m, know );
  }
  else if ( problem == 3 )
  {
    x = p03_sol ( m, know );
  }
  else if ( problem == 4 )
  {
    x = p04_sol ( m, know );
  }
  else if ( problem == 5 )
  {
    x = p05_sol ( m, know );
  }
  else if ( problem == 6 )
  {
    x = p06_sol ( m, know );
  }
  else if ( problem == 7 )
  {
    x = p07_sol ( m, know );
  }
  else if ( problem == 8 )
  {
    x = p08_sol ( m, know );
  }
  else if ( problem == 9 )
  {
    x = p09_sol ( m, know );
  }
  else if ( problem == 10 )
  {
    x = p10_sol ( m, know );
  }
  else if ( problem == 11 )
  {
    x = p11_sol ( m, know );
  }
  else if ( problem == 12 )
  {
    x = p12_sol ( m, know );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_SOL - Fatal error!\n";
    cerr << "  Illegal value of PROBLEM = " << problem << "\n";
    exit ( 1 );
  }
  return x;
}
//****************************************************************************80

string p00_title ( int problem )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns a title for any problem.
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
//    Input, int PROBLEM, the number of the problem.
//
//    Output, string P00_TITLE, a title for the problem.
//
{
  string title;

  if ( problem == 1 )
  {
    title = p01_title ( );
  }
  else if ( problem == 2 )
  {
    title = p02_title ( );
  }
  else if ( problem == 3 )
  {
    title = p03_title ( );
  }
  else if ( problem == 4 )
  {
    title = p04_title ( );
  }
  else if ( problem == 5 )
  {
    title = p05_title ( );
  }
  else if ( problem == 6 )
  {
    title = p06_title ( );
  }
  else if ( problem == 7 )
  {
    title = p07_title ( );
  }
  else if ( problem == 8 )
  {
    title = p08_title ( );
  }
  else if ( problem == 9 )
  {
    title = p09_title ( );
  }
  else if ( problem == 10 )
  {
    title = p10_title ( );
  }
  else if ( problem == 11 )
  {
    title = p11_title ( );
  }
  else if ( problem == 12 )
  {
    title = p12_title ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_TITLE - Fatal error!\n";
    cerr << "  Illegal value of PROBLEM = " << problem << "\n";
    exit ( 1 );
  }
  return title;
}
//****************************************************************************80

void p01_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_AB evaluates the limits of the optimization region for problem 01.
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
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -5.0;
    b[i] = +5.0;
  }
  return;
}
//****************************************************************************80

double *p01_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the objective function for problem 01.
//
//  Discussion:
//
//    The function is continuous, convex, and unimodal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hugues Bersini, Marco Dorigo, Stefan Langerman, Gregory Seront, 
//    Luca Gambardella,
//    Results of the first international contest on evolutionary optimisation,
//    In Proceedings of 1996 IEEE International Conference on Evolutionary 
//    Computation,
//    IEEE Press, pages 611-615, 1996.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P01_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] + pow ( x[i+j*m] - 1.0, 2 );
    }
  }
  return f;
}
//****************************************************************************80

double *p01_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P01_SOL returns the solution for problem 01.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P01_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 1.0;
    }
  }
  else 
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns a title for problem 01.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "The sphere model.";

  return title;
}
//****************************************************************************80

void p02_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_AB evaluates the limits of the optimization region for problem 02.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -5.12;
    b[i] = +5.12;
  }
  return;
}
//****************************************************************************80

double *p02_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the objective function for problem 02.
//
//  Discussion:
//
//    This function is also known as the weighted sphere model.
//
//    The function is continuous, convex, and unimodal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P02_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];
  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] + ( double ) ( i + 1 ) * pow ( x[i+j*m], 2 );
    }
  }
  return f;
}
//****************************************************************************80

double *p02_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P02_SOL returns the solution for problem 02.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P02_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns a title for problem 02.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "The axis-parallel hyper-ellipsoid function.";

  return title;
}
//****************************************************************************80

void p03_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_AB evaluates the limits of the optimization region for problem 03.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -65.536;
    b[i] = +65.536;
  }
  return;
}
//****************************************************************************80

double *p03_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the objective function for problem 03.
//
//  Discussion:
//
//    This function is also known as the weighted sphere model.
//
//    The function is continuous, convex, and unimodal.
//
//     There is a typographical error in Molga and Smutnicki, so that the
//     formula for this function is given incorrectly.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P03_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;
  double x_sum;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    x_sum = 0.0;
    for ( i = 0; i < m; i++ )
    {
      x_sum = x_sum + x[i+j*m];
      f[j] = f[j] + x_sum * x_sum;
    }
  }
  return f;
}
//****************************************************************************80

double *p03_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P03_SOL returns the solution for problem 03.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P03_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns a title for problem 03.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "The rotated hyper-ellipsoid function.";

  return title;
}
//****************************************************************************80

void p04_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_AB evaluates the limits of the optimization region for problem 04.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -2.048;
    b[i] = +2.048;
  }
  return;
}
//****************************************************************************80

double *p04_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the objective function for problem 04.
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
//  Reference:
//
//    Howard Rosenbrock,
//    An Automatic Method for Finding the Greatest or Least Value of a Function,
//    Computer Journal,
//    Volume 3, 1960, pages 175-184.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P04_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] + pow ( 1.0 - x[i+j*m], 2 );
    }
    for ( i = 0; i < m - 1; i++ )
    {
      f[j] = f[j] + pow ( x[i+1+j*m] - x[i+j*m], 2 );
    }
  }
  return f;
}
//****************************************************************************80

double *p04_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P04_SOL returns the solution for problem 04.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P04_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 1.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns a title for problem 04.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Rosenbrock's valley.";

  return title;
}
//****************************************************************************80

void p05_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_AB evaluates the limits of the optimization region for problem 05.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -5.12;
    b[i] = +5.12;
  }
  return;
}
//****************************************************************************80

double *p05_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the objective function for problem 05.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P05_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;
  double pi = 3.141592653589793;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = ( double ) ( 10 * m );
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] + pow ( x[i+j*m], 2 ) - 10.0 * cos ( 2.0 * pi * x[i+j*m] );
    }
  }
  return f;
}
//****************************************************************************80

double *p05_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P05_SOL returns the solution for problem 05.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P05_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns a title for problem 05.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Rastrigin's function.";

  return title;
}
//****************************************************************************80

void p06_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_AB evaluates the limits of the optimization region for problem 06.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -500.0;
    b[i] = +500.0;
  }
  return;
}
//****************************************************************************80

double *p06_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the objective function for problem 06.
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
//  Reference:
//
//    Hans-Paul Schwefel,
//    Numerical optimization of computer models,
//    Wiley, 1981,
//    ISBN13: 978-0471099888,
//    LC: QA402.5.S3813.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P06_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] - x[i+j*m] * sin ( sqrt ( r8_abs ( x[i+j*m] ) ) );
    }
  }
  return f;
}
//****************************************************************************80

double *p06_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P06_SOL returns the solution for problem 06.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P06_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 420.9687;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns a title for problem 06.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Schwefel's function.";

  return title;
}
//****************************************************************************80

void p07_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_AB evaluates the limits of the optimization region for problem 07.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -600.0;
    b[i] = +600.0;
  }
  return;
}
//****************************************************************************80

double *p07_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_F evaluates the objective function for problem 07.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P07_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;
  double p;
  double s;
  double *y;

  y = r8vec_indicator_new ( m );

  for ( i = 0; i < m; i++ )
  {
    y[i] = sqrt ( y[i] );
  }

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    p = 1.0;
    for ( i = 0; i < m; i++ )
    {
      s = s + pow ( x[i+j*m], 2 );
      p = p * cos ( x[i+j*m] / y[i] );
    }
    f[j] = s / 4000.0 - p + 1.0;
  }
  delete [] y;

  return f;
}
//****************************************************************************80

double *p07_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P07_SOL returns the solution for problem 07.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P07_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns a title for problem 07.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Griewank's function.";

  return title;
}
//****************************************************************************80

void p08_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_AB evaluates the limits of the optimization region for problem 08.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -1.0;
    b[i] = +1.0;
  }
  return;
}
//****************************************************************************80

double *p08_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_F evaluates the objective function for problem 08.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P08_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      f[j] = f[j] + pow ( r8_abs ( x[i+j*m] ), i + 2 );
    }
  }
  return f;
}
//****************************************************************************80

double *p08_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P08_SOL returns the solution for problem 08.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P08_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns a title for problem 08.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "The power sum function.";

  return title;
}
//****************************************************************************80

void p09_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_AB evaluates the limits of the optimization region for problem 09.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -32.768;
    b[i] = +32.768;
  }
  return;
}
//****************************************************************************80

double *p09_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_F evaluates the objective function for problem 09.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P09_F[N], the function evaluated at the arguments.
//
{
  double a = 20.0;
  double b = 0.2;
  double c = 0.2;
  double *f;
  int i;
  int j;
  double pi = 3.141592653589793;
  double s1;
  double s2;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    s1 = 0.0;
    s2 = 0.0;
    for ( i = 0; i < m; i++ )
    {
      s1 = s1 + pow ( x[i+j*m], 2 );
      s2 = s2 + cos ( c * pi * x[i+j*m] );
    }
    f[j] = - a * exp ( - b * sqrt ( s1 / ( double ) ( m ) ) )
      - exp ( s2 / ( double ) ( m ) ) + a + exp ( 1.0 );
  }
  return f;
}
//****************************************************************************80

double *p09_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P09_SOL returns the solution for problem 09.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P09_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_TITLE returns a title for problem 09.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Ackley's function.";

  return title;
}
//****************************************************************************80

void p10_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_AB evaluates the limits of the optimization region for problem 10.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;
  double pi = 3.141592653589793;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = pi;
  }
  return;
}
//****************************************************************************80

double *p10_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_F evaluates the objective function for problem 10.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P10_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;
  int p = 10;
  double pi = 3.141592653589793;
  double s;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    s = 0.0;
    for ( i = 0; i < m; i++ )
    {
      s = s - sin ( x[i+j*m] ) 
        * pow ( sin ( x[i+j*m] * x[i+j*m] * ( double ) ( i + 1 ) / pi ), 2 * p );
    }
    f[j] = s;
  }
  return f;
}
//****************************************************************************80

double *p10_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P10_SOL returns the solution for problem 10.
//
//  Discussion:
//
//    The minimum value is - 0.966 * M.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P10_SOL[M], the solution, if known.
//
{
  double *x;

  know = 0;
  x = NULL;

  return x;
}
//****************************************************************************80

string p10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_TITLE returns a title for problem 10.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Michalewicz's function.";

  return title;
}
//****************************************************************************80

void p11_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_AB evaluates the limits of the optimization region for problem 11.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = -5.12;
    b[i] = +5.12;
  }
  return;
}
//****************************************************************************80

double *p11_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P11_F evaluates the objective function for problem 11.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P11_F[N], the function evaluated at the arguments.
//
{
  double *f;
  int i;
  int j;
  double rsq;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    rsq = 0.0;
    for ( i = 0; i < m; i++ )
    {
      rsq = rsq + pow ( x[i+j*m], 2 );
    }
    f[j] = - ( 1.0 + cos ( 12.0 * sqrt ( rsq ) ) ) / ( 0.5 * rsq + 2.0 );
  }
  return f;
}
//****************************************************************************80

double *p11_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P11_SOL returns the solution for problem 11.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P01_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    know = 1;
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = 0.0;
    }
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p11_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P11_TITLE returns a title for problem 11.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "Drop wave function.";

  return title;
}
//****************************************************************************80

void p12_ab ( int m, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_AB evaluates the limits of the optimization region for problem 12.
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Output, double A[M], B[M], the lower and upper bounds.
//
{
  int i;

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
  return;
}
//****************************************************************************80

double *p12_f ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P12_F evaluates the objective function for problem 12.
//
//  Discussion:
//
//    In dimension I, the function is a piecewise linear function with
//    local minima at 0 and 1.0, and a global minimum at ALPHA(I) = I/(M+1).
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
//  Reference:
//
//    Marcin Molga, Czeslaw Smutnicki,
//    Test functions for optimization needs.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of arguments.
//
//    Input, double X[M*N], the arguments.
//
//    Output, double P12_F[N], the function evaluated at the arguments.
//
{
  double *alpha;
  double beta = 2.0;
  double *f;
  double g;
  int i;
  int j;

  alpha = new double[m];

  for ( i = 0; i < m; i++ )
  {
    alpha[i] = ( double ) ( i + 1 ) / ( double ) ( m + 1 );
  }

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      if ( x[i+j*m] <= 0.0 )
      {
        g = x[i+j*m];
      }
      else if ( x[i+j*m] <= 0.8 * alpha[i] )
      {
        g = 0.8 - x[i+j*m] / alpha[i];
      }
      else if ( x[i+j*m] <= alpha[i] )
      {
        g = 5.0 * x[i+j*m] / alpha[i] - 4.0;
      }
      else if ( x[i+j*m] <= ( 1.0 + 4.0 * alpha[i] ) / 5.0 )
      {
        g = 1.0 + 5.0 * ( x[i+j*m] - alpha[i] ) / ( alpha[i] - 1.0 );
      }
      else if ( x[i+j*m] <= 1.0 )
      {
        g = 0.8 + ( x[i+j*m] - 1.0 ) / ( 1.0 - alpha[i] );
      }
      else
      {
        g = x[i+j*m] - 1.0;
      }
      f[j] = f[j] + g;
    }
    f[j] = f[j] / ( double ) ( m );
    f[j] = - pow ( f[j], beta );
  }
  delete [] alpha;

  return f;
}
//****************************************************************************80

double *p12_sol ( int m, int &know )

//****************************************************************************80
//
//  Purpose:
//
//    P12_SOL returns the solution for problem 12.
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
//    Input, int M, the spatial dimension.
//
//    Input/output, int &KNOW.
//    On input, KNOW is 0, or the index of the previously returned solution.
//    On output, KNOW is 0 if there are no more solutions, or it is the
//    index of the next solution.
//
//    Output, double P12_SOL[M], the solution, if known.
//
{
  int i;
  double *x;

  if ( know == 0 )
  {
    x = new double[m];
    for ( i = 0; i < m; i++ )
    {
      x[i] = ( double ) ( i + 1 ) / ( double ) ( m + 1 );
    }
    know = 1;
  }
  else
  {
    know = 0;
    x = NULL;
  }
  return x;
}
//****************************************************************************80

string p12_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P12_TITLE returns a title for problem 12.
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
//    Output, string TITLE, a title for the problem.
//
{
  string title;

  title = "The deceptive function.";

  return title;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double *r8col_uniform_new ( int m, int n, double a[], double b[], int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8COL_UNIFORM_NEW fills an R8COL with scaled pseudorandom numbers.
//
//  Discussion:
//
//    An R8COL is an array of R8 values, regarded as a set of column vectors.
//
//    The user specifies a minimum and maximum value for each row.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 December 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M], B[M], the upper and lower limits.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8COL_UNIFORM_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
      r[i+j*m] = a[i] 
        + ( b[i] - a[i] ) * ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double *r8vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_INDICATOR_NEW sets an R8VEC to the indicator vector {1,2,3...}.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, double R8VEC_INDICATOR_NEW[N], the indicator array.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i <= n-1; i++ )
  {
    a[i] = ( double ) ( i + 1 );
  }

  return a;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
