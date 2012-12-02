# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_interp_nd.hpp"
# include "r8lib.hpp"

//****************************************************************************80

 double csevl ( double x, double a[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, double CSEVL, the Chebyshev series evaluated at X.
//
{
  double b0;
  double b1;
  double b2;
  int i;
  double twox;
  double value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CSEVL - Fatal error!\n";
    cerr << "  Number of terms <= 0.\n";
    exit ( 1 );
  }

  if ( 1000 < n )
  {
    cerr << "\n";
    cerr << "CSEVL - Fatal error!\n";
    cerr << "  Number of terms greater than 1000.\n";
    exit ( 1 );
 }

  if ( x < -1.1 || 1.1 < x )
  {
    cerr << "\n";
    cerr << "CSEVL - Fatal error!\n";
    cerr << "  X outside (-1,+1).\n";
    exit ( 1 );
  }

  twox = 2.0 * x;
  b1 = 0.0;
  b0 = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[i];
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

int i4vec_sum ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SUM sums the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_SUM = 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector to be summed.
//
//    Output, int I4VEC_SUM, the sum of the entries of A.
//
{
  int i;
  int sum;

  sum = 0;
  for ( i = 0; i < n; i++ )
  {
    sum = sum + a[i];
  }

  return sum;
}
//****************************************************************************80

int inits ( double dos[], int nos, double eta )

//****************************************************************************80
//
//  Purpose:
//
//    INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, double ETA, the desired accuracy.
//
//    Output, int INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
{
  double err;
  int i;
  int value;

  if ( nos < 1 )
  {
    cerr << "\n";
    cerr << "INITS - Fatal error!\n";
    cerr << "  Number of coefficients < 1.\n";
    exit ( 1 );
  }

  err = 0.0;

  for ( i = nos - 1; 0 <= i; i-- )
  {
    err = err + r8_abs ( dos[i] );
    if ( eta < err )
    {
      value = i + 1;
      return value;
    }
  }

  value = i;
  cerr << "\n";
  cerr << "INITS - Warning!\n";
  cerr << "  ETA may be too small.\n";

  return value;
}
//****************************************************************************80

double *p00_c ( int prob, int m, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    P00_CW computes a random C parameter vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double P00_C[M], the parameter vector.
//
{
  double b[6] = { 1.5, 0.0, 1.85, 7.03, 20.4, 4.3 };
  double *c;
  double c_sum;
  int i;

  b[2-1] = ( double ) ( m );

  c = r8vec_uniform_01_new ( m, seed );
  c_sum = r8vec_sum ( m, c );

  for ( i = 0; i < m; i++ )
  {
    c[i] = b[prob-1] * c[i] / c_sum;
  }
  return c;
}
//****************************************************************************80

double *p00_d ( int prob, int m, int id, double c[], double w[], int n, 
  double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_D returns a derivative component of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the evalution points.
//
//    Output, double P00_D[N], the ID-th derivative component.
//
{
  double *d;

  if ( id < 0 || m < id )
  {
    cerr << "\n";
    cerr << "P00_D - Fatal error!\n";
    cerr << "  Illegal spatial coordinate ID = " << id << "\n";
    exit ( 1 );
  }

  if ( prob == 1 )
  {
    d = p01_d ( m, id, c, w, n, x );
  }
  else if ( prob == 2 )
  {
    d = p02_d ( m, id, c, w, n, x );
  }
  else if ( prob == 3 )
  {
    d = p03_d ( m, id, c, w, n, x );
  }
  else if ( prob == 4 )
  {
    d = p04_d ( m, id, c, w, n, x );
  }
  else if ( prob == 5 )
  {
    d = p05_d ( m, id, c, w, n, x );
  }
  else if ( prob == 6 )
  {
    d = p06_d ( m, id, c, w, n, x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_D - Fatal error!\n";
    cerr << "  Illegal function index PROB = " << prob << "\n";
    exit ( 1 );
  }

  return d;
}
//****************************************************************************80

double *p00_f ( int prob, int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_F returns the value of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[M*N], the evalution points.
//
//    Output, double P00_F[N], the function values.
//
{
  double *f;

  if ( prob == 1 )
  {
    f = p01_f ( m, c, w, n, x );
  }
  else if ( prob == 2 )
  {
    f = p02_f ( m, c, w, n, x );
  }
  else if ( prob == 3 )
  {
    f = p03_f ( m, c, w, n, x );
  }
  else if ( prob == 4 )
  {
    f = p04_f ( m, c, w, n, x );
  }
  else if ( prob == 5 )
  {
    f = p05_f ( m, c, w, n, x );
  }
  else if ( prob == 6 )
  {
    f = p06_f ( m, c, w, n, x );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_F - Fatal error!\n";
    cerr << "  Illegal function index PROB = " << prob << "\n";
    exit ( 1 );
  }

  return f;
}

//****************************************************************************80

int p00_prob_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of test functions available.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//   Output, int P00_PROB_NUM, the number of test functions.
//
{
  int prob_num;

  prob_num = 6;

  return prob_num;
}
//****************************************************************************80

double p00_q ( int prob, int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_Q returns the integral of any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Input, int M, the spatial dimension.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P00_Q, the integral.
//
{
  double q;

  if ( prob == 1 )
  {
    q = p01_q ( m, c, w );
  }
  else if ( prob == 2 )
  {
    q = p02_q ( m, c, w );
  }
  else if ( prob == 3 )
  {
    q = p03_q ( m, c, w );
  }
  else if ( prob == 4 )
  {
    q = p04_q ( m, c, w );
  }
  else if ( prob == 5 )
  {
    q = p05_q ( m, c, w );
  }
  else if ( prob == 6 )
  {
    q = p06_q ( m, c, w );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_Q - Fatal error!\n";
    cerr << "  Illegal function index PROB = " << prob << "\n";
    exit ( 1 );
  }

  return q;
}
//****************************************************************************80

string p00_title ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title for any function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the index of the function.
//
//    Output, string P00_TITLE, the function title.
//
{
  string title;

  if ( prob == 1 )
  {
    title = p01_title ( );
  }
  else if ( prob == 2 )
  {
    title = p02_title ( );
  }
  else if ( prob == 3 )
  {
    title = p03_title ( );
  }
  else if ( prob == 4 )
  {
    title = p04_title ( );
  }
  else if ( prob == 5 )
  {
    title = p05_title ( );
  }
  else if ( prob == 6 )
  {
    title = p06_title ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_TITLE - Fatal error!\n";
    cerr << "  Illegal function index PROB = " << prob << "\n";
    exit ( 1 );
  }
  return title;
}
//****************************************************************************80

double *p00_w ( int prob, int m, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    P00_W computes a random W parameter vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem number.
//
//    Input, int M, the spatial dimension.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, double P00_W[M], the parameter vector.
//
{
  double *w;

  w = r8vec_uniform_01_new ( m, seed );

  return w;
}
//****************************************************************************80

double *p01_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_D evaluates any derivative component for problem p01.
//
//  Discussion:
//
//    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
//
//    Default values are:
//
//    c(1:m) = 1/m
//    w(1) = 0.3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P01_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;
  double pi = 3.141592653589793;

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 2.0 * pi * w[0];
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = d[j] + c[i] * x[i+j*m];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    d[j] = - c[id] * sin ( d[j] );
  }
  return d;
}
//****************************************************************************80

double *p01_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_F evaluates the function for problem p01.
//
//  Discussion:
//
//    f(x) = cos ( 2 * pi * w(1) + sum ( c(1:m) * x(1:m) ) )
//
//    Default values are:
//
//    c(1:m) = 1/m
//    w(1) = 0.3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P01_F[N], the function values.
//
{
  double *f;
  int i;
  int j;
  double pi = 3.141592653589793;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 2.0 * pi * w[0];
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = f[j] + c[i] * x[i+j*m];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    f[j] = cos ( f[j] );
  }
  return f;
}
//****************************************************************************80

double p01_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_Q evaluates the integral for problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P01_Q, the integral.
//
{
  double c_prod;
  double c_sum;
  int i;
  double pi = 3.141592653589793;
  double q;

  c_sum = r8vec_sum ( m, c );

  c_prod = 1.0;
  for ( i = 0; i < m; i++ )
  {
    c_prod = c_prod * sin ( 0.5 * c[i] ) / c[i];
  }

  q = pow ( 2.0, m ) * cos ( 2.0 * pi * w[0] + 0.5 * c_sum ) * c_prod;

  return q;
}
//****************************************************************************80

string p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the name of problem p01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P01_TITLE, the title of the problem.
//
{
  string title;

  title = "Oscillatory";

  return title;
}
//****************************************************************************80

double *p02_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_D evaluates an derivative component for problem p02.
//
//  Discussion:
//
//    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P02_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = d[j] * ( pow ( c[i], - 2 ) + pow ( x[i+j*m] - w[i], 2 ) );
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = - 2.0 / d[j] * ( x[id+j*m] - w[id] ) / 
      ( pow ( c[id], - 2 ) + pow ( x[id+j*m] - w[id], 2 ) );
  }

  return d;
}
//****************************************************************************80

double *p02_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_F evaluates the function for problem p02.
//
//  Discussion:
//
//    f(x) = 1 / product ( c(1:m)^(-2) + ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P02_F[N], the function values.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = f[j] * ( pow ( c[i], - 2 ) + pow ( x[i+j*m] - w[i], 2 ) );
    }
  }

  for ( j = 0; j < n; j++ )
  {
    f[j] = 1.0 / f[j];
  }

  return f;
}
//****************************************************************************80

double p02_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_Q evaluates the integral for problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P02_Q, the integral.
//
{
  int i;
  double q;

  q = 1.0;

  for ( i = 0; i < m; i++ )
  {
    q = q * 
      ( atan ( ( 1.0 - w[i] ) * c[i] ) 
      + atan (         w[i]   * c[i] ) 
      ) * c[i]; 
  }

  return q;
}
//****************************************************************************80

string p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns the title of problem p02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P02_TITLE, the title of the problem.
//
{
  string title;

  title = "Product Peak";

  return title;
}
//****************************************************************************80

double *p03_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_D evaluates any derivative component for problem p03.
//
//  Discussion:
//
//    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
//
//    Default values are:
//
//    c(1:m) = 1/m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P03_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 1.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = d[j] + c[i] * x[i+j*m];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = - c[id] * ( double ) ( m + 1 ) / pow ( d[j], m + 2 );
  }
  return d;
}
//****************************************************************************80

double *p03_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_F evaluates the function for problem p03.
//
//  Discussion:
//
//    f(x) = 1 / ( 1 + sum ( c(1:m) * x(1:m) ) ) ^ ( m + 1 )
//
//    Default values are:
//
//    c(1:m) = 1/m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P03_F[N], the function values.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 1.0;
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = f[j] + c[i] * x[i+j*m];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    f[j] = 1.0 / pow ( f[j], m + 1 );
  }
  return f;
}
//****************************************************************************80

double p03_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_Q evaluates the integral for problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P03_Q, the integral.
//
{
  int *ivec;
  double q;
  int rank;
  int s;
//
//  Here, we need to generate all possible DIM_NUM tuples with
//  values of 0 or 1.
//
  ivec = new int[m];

  q = 0.0;
  rank = 0;

  while ( 1 )
  {
    tuple_next ( 0, 1, m, rank, ivec );

    if ( rank == 0 )
    {
      break;
    }

    s = i4vec_sum ( m, ivec );

    q = q + r8_mop ( s ) / ( 1.0 + r8vec_i4vec_dot_product ( m, c, ivec ) );
  }

  q = q / ( r8_factorial ( m ) * r8vec_product ( m, c ) );

  delete [] ivec;

  return q;
}
//****************************************************************************80

string p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns the title of problem p03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P03_TITLE, the title of the problem.
//
{
  string title;

  title = "Corner Peak";

  return title;
}
//****************************************************************************80

double *p04_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_D evaluates any derivative component for problem p04.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1 / m
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P04_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 0.0;
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = d[j] + pow ( c[i] * ( x[i+j*m] - w[i] ), 2 );
    }
  }
  for ( j = 0; j < n; j++ )
  {
    d[j] = exp ( - d[j] ) * pow ( c[id], 2 ) * ( - 2.0 ) * ( x[id+j*m] - w[id] );
  }
  return d;
}
//****************************************************************************80

double *p04_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_F evaluates the function for problem p04.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m)^2 * ( x(1:m) - w(1:m) )^2 )
//
//    Default values are:
//
//    c(1:m) = 1 / m
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P04_F[N], the function values.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = f[j] + pow ( c[i] * ( x[i+j*m] - w[i] ), 2 );
    }
  }
  for ( j = 0; j < n; j++ )
  {
    f[j] = exp ( - f[j] );
  }
  return f;
}
//****************************************************************************80

double p04_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_Q evaluates the integral for problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P04_Q, the integral.
//
{
  int i;
  double pi = 3.141592653589793;
  double q;

  q = 1.0;

  for ( i = 0; i < m; i++ )
  {
    q = q * sqrt ( pi ) 
      * ( r8_error ( c[i] * ( 1.0 - w[i] ) ) 
        + r8_error ( c[i] *         w[i] ) ) 
      / ( 2.0 * c[i] );
  }

  return q;
}
//****************************************************************************80

string p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns the title of problem p04.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P04_TITLE, the title of the problem.
//
{
  string title;

  title = "Gaussian";

  return title;
}
//****************************************************************************80

double *p05_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_D evaluates any derivative component for problem p05.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
//
//    Default values are:
//
//    c(1:m) = 2.0
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P05_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 0.0;
  }

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = d[j] + c[i] * r8_abs ( x[i+j*m] - w[i] );
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = exp ( - d[j] );
  }

  for ( j = 0; j < n; j++ )
  {
    if ( x[id+j*m] - w[id] <= 0.0 )
    {
      d[j] = d[j] * c[id];
    }
    else
    {
      d[j] = - d[j] * c[id];
    }
  }

  return d;
}
//****************************************************************************80

double *p05_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_F evaluates the function for problem p05.
//
//  Discussion:
//
//    f(x) = exp ( - sum ( c(1:m) * abs ( x(1:m) - w(1:m) ) ) )
//
//    Default values are:
//
//    c(1:m) = 2.0
//    w(1:m) = 0.5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P05_F[N], the function values.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  for ( j = 0; j < n; j++ )
  {
    f[j] = 0.0;
  }
  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = f[j] + c[i] * r8_abs ( x[i+j*m] - w[i] );
    }
  }
  for ( j = 0; j < n; j++ )
  {
    f[j] = exp ( - f[j] );
  }
  return f;
}
//****************************************************************************80

double p05_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_Q evaluates the integral for problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P05_Q, the integral.
//
{
  int i;
  double q;

  q = 1.0;

  for ( i = 0; i < m; i++ )
  {
//
//  W < 0 < 1
//
//  | X - W | = X - W from 0 to 1.
//
    if ( w[i] < 0.0 )
    {
      q = q * 
        ( exp ( - c[i] * (     - w[i] ) ) 
        - exp ( - c[i] * ( 1.0 - w[i] ) ) ) / c[i];
    }
//
//  0 < W < 1
//
//  | X - W | = W - X from 0 to Z, 
//            = X - W from      Z to 1.
//
    else if ( w[i] < 1.0 )
    {
      q = q * ( 2.0 
          - exp ( - c[i] * (       w[i] ) ) 
          - exp ( - c[i] * ( 1.0 - w[i] ) ) ) / c[i];
    }
//
//  0 < 1 < W
//
//  | X - W | = W - X from 0 to 1.
//
    else
    {
      q = q * 
        ( exp ( - c[i] * ( w[i] - 1.0 ) ) 
        - exp ( - c[i] * ( w[i]       ) ) ) / c[i];

    }
  }

  return q;
}
//****************************************************************************80

string p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns the title of problem p05.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P05_TITLE, the title of the problem.
//
{
  string title;

  title = "Continuous";

  return title;
}
//****************************************************************************80

double *p06_d ( int m, int id, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_D evaluates any derivative component for problem p06.
//
//  Discussion:
//
//    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
//           0                          otherwise
//
//    Default values are:
//
//    c(1:m) = 0.5^(1/m)
//    w(1:2) = 0.5^(1/m)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, int ID, the spatial coordinate to differentiate.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P06_D[N], the ID-th derivative component.
//
{
  double *d;
  int i;
  int j;

  d = new double[n];

  if ( m == 1 )
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = c[0] * exp ( c[0] * x[0+j*m] );
    }

    for ( j = 0; j < n; j++ )
    {
      if ( w[0] < x[0+j*m] )
      {
        d[j] = 0.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      d[j] = 0.0;
    }
    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        d[j] = d[j] + c[i] * x[i+j*m];
      }
    }
    for ( j = 0; j < n; j++ )
    {
      d[j] = c[id] * exp ( d[j] );
    }
    for ( j = 0; j < n; j++ )
    {
      if ( w[0] < x[0+j*m] || w[1] < x[1+j*m] )
      {
        d[j] = 0.0;
      }
    }
  }

  return d;
}
//****************************************************************************80

double *p06_f ( int m, double c[], double w[], int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_F evaluates the function for problem p06.
//
//  Discussion:
//
//    f(x) = exp ( c(1:m) * x(1:m) ) if x(1) <= w(1) and x(2) <= w(2).
//           0                          otherwise
//
//    Default values are:
//
//    c(1:m) = 0.5^(1/m)
//    w(1:2) = 0.5^(1/m)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the evaluation points.
//
//    Output, double P06_F[N], the function values.
//
{
  double *f;
  int i;
  int j;

  f = new double[n];

  if ( m == 1 )
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = exp ( c[0] * x[0+j*m] );
    }
    for ( j = 0; j < n; j++ )
    {
      if ( w[0] < x[0+j*m] )
      {
        f[j] = 0.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      f[j] = 0.0;
    }
    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        f[j] = f[j] + c[i] * x[i+j*m];
      }
    }
    for ( j = 0; j < n; j++ )
    {
      f[j] = exp ( f[j] );
    }
    for ( j = 0; j < n; j++ )
    {
      if ( w[0] < x[0+j*m] || w[1] < x[1+j*m] )
      {
        f[j] = 0.0;
      }
    }
  }
  return f;
}
//****************************************************************************80

double p06_q ( int m, double c[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_Q evaluates the integral for problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Alan Genz,
//    A Package for Testing Multiple Integration Subroutines,
//    in Numerical Integration: Recent Developments, Software
//    and Applications,
//    edited by Patrick Keast and Graeme Fairweather,
//    Reidel, 1987, pages 337-340,
//    ISBN: 9027725144,
//    LC: QA299.3.N38.
//
//  Parameters:
//
//    Input, int M, the dimension of the argument.
//
//    Input, double C[M], W[M], the problem parameters.
//
//    Output, double P06_Q, the integral.
//
{
  int i;
  double q;
//
//  To simplify the calculation, force W(3:M) to be at least 1.0.
//
  for ( i = 2; i < m; i++ )
  {
    w[i] = 1.0;
  }

  q = 1.0;

  for ( i = 0; i < m; i++ )
  {
    if ( w[i] <= 0.0 )
    {
      q = q * 0.0;
    }
    else if ( w[i] <= 1.0 )
    {
      if ( c[i] == 0.0 )
      {
        q = q * w[i];
      }
      else
      {
        q = q * ( exp ( c[i] * w[i] ) - 1.0 ) / c[i];
      }
    }
    else
    {
      if ( c[i] != 0.0 )
      {
        q = q * ( exp ( c[i] * w[i] ) - 1.0 ) / c[i];
      }
    }
  }

  return q;
}
//****************************************************************************80

string p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns the title of problem p06.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 August 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string P06_TITLE, the title of the problem.
//
{
  string title;

  title = "Discontinuous";

  return title;
}
//****************************************************************************80

double r8_error ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERROR evaluates the error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERROR, the error function of X.
//
{
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  static int nterf = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xbig = 0.0;
  double y;

  if ( nterf == 0 )
  {
    nterf = inits ( erfcs, 21, 0.1 * r8_mach ( 3 ) );
    xbig = sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= sqeps )
  {
    value = 2.0 * x / sqrtpi;
  }
  else if ( y <= 1.0 )
  {
    value = x * ( 1.0 + csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
  }
  else if ( y <= xbig )
  {
    value = 1.0 - r8_errorc ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = 1.0;
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_errorc ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ERRORC evaluates the co-error function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_ERRORC, the co-error function of X.
//
{
  static double erc2cs[49] = {
    -0.6960134660230950112739150826197E-01,
    -0.4110133936262089348982212084666E-01,
    +0.3914495866689626881561143705244E-02,
    -0.4906395650548979161280935450774E-03,
    +0.7157479001377036380760894141825E-04,
    -0.1153071634131232833808232847912E-04,
    +0.1994670590201997635052314867709E-05,
    -0.3642666471599222873936118430711E-06,
    +0.6944372610005012589931277214633E-07,
    -0.1371220902104366019534605141210E-07,
    +0.2788389661007137131963860348087E-08,
    -0.5814164724331161551864791050316E-09,
    +0.1238920491752753181180168817950E-09,
    -0.2690639145306743432390424937889E-10,
    +0.5942614350847910982444709683840E-11,
    -0.1332386735758119579287754420570E-11,
    +0.3028046806177132017173697243304E-12,
    -0.6966648814941032588795867588954E-13,
    +0.1620854541053922969812893227628E-13,
    -0.3809934465250491999876913057729E-14,
    +0.9040487815978831149368971012975E-15,
    -0.2164006195089607347809812047003E-15,
    +0.5222102233995854984607980244172E-16,
    -0.1269729602364555336372415527780E-16,
    +0.3109145504276197583836227412951E-17,
    -0.7663762920320385524009566714811E-18,
    +0.1900819251362745202536929733290E-18,
    -0.4742207279069039545225655999965E-19,
    +0.1189649200076528382880683078451E-19,
    -0.3000035590325780256845271313066E-20,
    +0.7602993453043246173019385277098E-21,
    -0.1935909447606872881569811049130E-21,
    +0.4951399124773337881000042386773E-22,
    -0.1271807481336371879608621989888E-22,
    +0.3280049600469513043315841652053E-23,
    -0.8492320176822896568924792422399E-24,
    +0.2206917892807560223519879987199E-24,
    -0.5755617245696528498312819507199E-25,
    +0.1506191533639234250354144051199E-25,
    -0.3954502959018796953104285695999E-26,
    +0.1041529704151500979984645051733E-26,
    -0.2751487795278765079450178901333E-27,
    +0.7290058205497557408997703680000E-28,
    -0.1936939645915947804077501098666E-28,
    +0.5160357112051487298370054826666E-29,
    -0.1378419322193094099389644800000E-29,
    +0.3691326793107069042251093333333E-30,
    -0.9909389590624365420653226666666E-31,
    +0.2666491705195388413323946666666E-31 };
  static double erfccs[59] = {
    +0.715179310202924774503697709496E-01,
    -0.265324343376067157558893386681E-01,
    +0.171115397792085588332699194606E-02,
    -0.163751663458517884163746404749E-03,
    +0.198712935005520364995974806758E-04,
    -0.284371241276655508750175183152E-05,
    +0.460616130896313036969379968464E-06,
    -0.822775302587920842057766536366E-07,
    +0.159214187277090112989358340826E-07,
    -0.329507136225284321486631665072E-08,
    +0.722343976040055546581261153890E-09,
    -0.166485581339872959344695966886E-09,
    +0.401039258823766482077671768814E-10,
    -0.100481621442573113272170176283E-10,
    +0.260827591330033380859341009439E-11,
    -0.699111056040402486557697812476E-12,
    +0.192949233326170708624205749803E-12,
    -0.547013118875433106490125085271E-13,
    +0.158966330976269744839084032762E-13,
    -0.472689398019755483920369584290E-14,
    +0.143587337678498478672873997840E-14,
    -0.444951056181735839417250062829E-15,
    +0.140481088476823343737305537466E-15,
    -0.451381838776421089625963281623E-16,
    +0.147452154104513307787018713262E-16,
    -0.489262140694577615436841552532E-17,
    +0.164761214141064673895301522827E-17,
    -0.562681717632940809299928521323E-18,
    +0.194744338223207851429197867821E-18,
    -0.682630564294842072956664144723E-19,
    +0.242198888729864924018301125438E-19,
    -0.869341413350307042563800861857E-20,
    +0.315518034622808557122363401262E-20,
    -0.115737232404960874261239486742E-20,
    +0.428894716160565394623737097442E-21,
    -0.160503074205761685005737770964E-21,
    +0.606329875745380264495069923027E-22,
    -0.231140425169795849098840801367E-22,
    +0.888877854066188552554702955697E-23,
    -0.344726057665137652230718495566E-23,
    +0.134786546020696506827582774181E-23,
    -0.531179407112502173645873201807E-24,
    +0.210934105861978316828954734537E-24,
    -0.843836558792378911598133256738E-25,
    +0.339998252494520890627359576337E-25,
    -0.137945238807324209002238377110E-25,
    +0.563449031183325261513392634811E-26,
    -0.231649043447706544823427752700E-26,
    +0.958446284460181015263158381226E-27,
    -0.399072288033010972624224850193E-27,
    +0.167212922594447736017228709669E-27,
    -0.704599152276601385638803782587E-28,
    +0.297976840286420635412357989444E-28,
    -0.126252246646061929722422632994E-28,
    +0.539543870454248793985299653154E-29,
    -0.238099288253145918675346190062E-29,
    +0.109905283010276157359726683750E-29,
    -0.486771374164496572732518677435E-30,
    +0.152587726411035756763200828211E-30 };
  static double erfcs[21] = {
    -0.49046121234691808039984544033376E-01,
    -0.14226120510371364237824741899631,
    +0.10035582187599795575754676712933E-01,
    -0.57687646997674847650827025509167E-03,
    +0.27419931252196061034422160791471E-04,
    -0.11043175507344507604135381295905E-05,
    +0.38488755420345036949961311498174E-07,
    -0.11808582533875466969631751801581E-08,
    +0.32334215826050909646402930953354E-10,
    -0.79910159470045487581607374708595E-12,
    +0.17990725113961455611967245486634E-13,
    -0.37186354878186926382316828209493E-15,
    +0.71035990037142529711689908394666E-17,
    -0.12612455119155225832495424853333E-18,
    +0.20916406941769294369170500266666E-20,
    -0.32539731029314072982364160000000E-22,
    +0.47668672097976748332373333333333E-24,
    -0.65980120782851343155199999999999E-26,
    +0.86550114699637626197333333333333E-28,
    -0.10788925177498064213333333333333E-29,
    +0.12811883993017002666666666666666E-31 };
  double eta;
  static int nterc2 = 0;
  static int nterf = 0;
  static int nterfc = 0;
  static double sqeps = 0.0;
  static double sqrtpi = 1.77245385090551602729816748334115;
  double value;
  static double xmax = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nterf == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nterf = inits ( erfcs, 21, eta );
    nterfc = inits ( erfccs, 59, eta );
    nterc2 = inits ( erc2cs, 49, eta );

    xsml = - sqrt ( - log ( sqrtpi * r8_mach ( 3 ) ) );
    xmax = sqrt (- log ( sqrtpi * r8_mach ( 1 ) ) );
    xmax = xmax - 0.5 * log ( xmax ) / xmax - 0.01;
    sqeps = sqrt ( 2.0 * r8_mach ( 3 ) );
  }

  if ( x <= xsml )
  {
    value = 2.0;
    return value;
  }

  if ( xmax < x )
  {
    cerr << "\n";
    cerr << "R8_ERRORC - Warning!\n";
    cerr << "  X so big that ERFC underflows.\n";
    value = 0.0;
    return value;
  }

  y = r8_abs ( x );

  if ( y < sqeps )
  {
    value = 1.0 - 2.0 * x / sqrtpi;
    return value;
  }
  else if ( y <= 1.0 )
  {
    value = 1.0 - x * ( 1.0 
      + csevl ( 2.0 * x * x - 1.0, erfcs, nterf ) );
    return value;
  }

  y = y * y;

  if ( y <= 4.0 )
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + csevl ( ( 8.0 / y - 5.0 ) / 3.0, erc2cs, nterc2 ) );
  }
  else 
  {
    value = exp ( - y ) / r8_abs ( x ) * ( 0.5 
      + csevl ( 8.0 / y - 1.0, erfccs, nterfc ) );
  }

  if ( x < 0.0 )
  {
    value = 2.0 - value;
  }

  return value;
}
//****************************************************************************80

double r8_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R8_MACH(3) = B^(-T), the smallest relative spacing.
//      R8_MACH(4) = B^(1-T), the largest relative spacing.
//      R8_MACH(5) = log10(B).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528:
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, double R8_MACH, the value of the chosen parameter.
//
{
  double value;

  if ( i == 1 )
  {
    value = 4.450147717014403E-308;
  }
  else if ( i == 2 )
  {
    value = 8.988465674311579E+307;
  }
  else if ( i == 3 )
  {
    value = 1.110223024625157E-016;
  }
  else if ( i == 4 )
  {
    value = 2.220446049250313E-016;
  }
  else if ( i == 5 )
  {
    value = 0.301029995663981E+000;
  }
  else
  {
    cerr << "\n";
    cerr << "R8_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

void tuple_next ( int m1, int m2, int n, int &rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int &RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  int i;
  int j;

  if ( m2 < m1 )
  {
    rank = 0;
    return;
  }

  if ( rank <= 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = m1;
    }
    rank = 1;
  }
  else
  {
    rank = rank + 1;
    i = n - 1;

    for ( ; ; )
    {

      if ( x[i] < m2 )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = m1;

      if ( i == 0 )
      {
        rank = 0;
        for ( j = 0; j < n; j++ )
        {
          x[j] = m1;
        }
        break;
      }
      i = i - 1;
    }
  }
  return;
}
