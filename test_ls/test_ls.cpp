# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "test_ls.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *p00_a ( int prob, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P00_A returns the matrix A for any least squares problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P00_A[M*N], the matrix.
//
{
  double *a;

  if ( prob == 1 )
  {
    a = p01_a ( m, n );
  }
  else if ( prob == 2 )
  {
    a = p02_a ( m, n );
  }
  else if ( prob == 3 )
  {
    a = p03_a ( m, n );
  }
  else if ( prob == 4 )
  {
    a = p04_a ( m, n );
  }
  else if ( prob == 5 )
  {
    a = p05_a ( m, n );
  }
  else if ( prob == 6 )
  {
    a = p06_a ( m, n );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_A - Fatal error!\n";
    cerr << "  Illegal value of PROB = " << prob << "\n";
    exit ( 1 );
  }

  return a;
}
//****************************************************************************80

double *p00_b ( int prob, int m )

//****************************************************************************80
//
//  Purpose:
//
//    P00_B returns the right hand side B for any least squares problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int M, the number of equations.
//
//    Output, double P00_B[M], the right hand side.
//
{
  double *b;

  if ( prob == 1 )
  {
    b = p01_b ( m );
  }
  else if ( prob == 2 )
  {
    b = p02_b ( m );
  }
  else if ( prob == 3 )
  {
    b = p03_b ( m );
  }
  else if ( prob == 4 )
  {
    b = p04_b ( m );
  }
  else if ( prob == 5 )
  {
    b = p05_b ( m );
  }
  else if ( prob == 6 )
  {
    b = p06_b ( m );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_B - Fatal error!\n";
    cerr << "  Illegal value of PROB = " << prob << "\n";
    exit ( 1 );
  }

  return b;
}
//****************************************************************************80

int p00_m ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_M returns the number of equations M for any least squares problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Output, int P00_M, the number of equations.
//
{
  int m;

  if ( prob == 1 )
  {
    m = p01_m ( );
  }
  else if ( prob == 2 )
  {
    m = p02_m ( );
  }
  else if ( prob == 3 )
  {
    m = p03_m ( );
  }
  else if ( prob == 4 )
  {
    m = p04_m ( );
  }
  else if ( prob == 5 )
  {
    m = p05_m ( );
  }
  else if ( prob == 6 )
  {
    m = p06_m ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_M - Fatal error!\n";
    cerr << "  Illegal value of PROB = " << prob << "\n";
    exit ( 1 );
  }

  return m;
}
//****************************************************************************80

int p00_n ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_N returns the number of variables N for any least squares problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Output, int P00_N, the number of variables.
//
{
  int n;

  if ( prob == 1 )
  {
    n = p01_n ( );
  }
  else if ( prob == 2 )
  {
    n = p02_n ( );
  }
  else if ( prob == 3 )
  {
    n = p03_n ( );
  }
  else if ( prob == 4 )
  {
    n = p04_n ( );
  }
  else if ( prob == 5 )
  {
    n = p05_n ( );
  }
  else if ( prob == 6 )
  {
    n = p06_n ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_N - Fatal error!\n";
    cerr << "  Illegal value of PROB = " << prob << "\n";
    exit ( 1 );
  }

  return n;
}
//****************************************************************************80

int p00_prob_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of least squares problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double P00_PROB_NUM, the number of problems.
//
{
  int prob_num;

  prob_num = 6;

  return prob_num;
}
//****************************************************************************80

double *p00_x ( int prob, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P00_X returns the least squares solution X for any least squares problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the problem index.
//
//    Input, int N, the number of variables.
//
//    Output, double P00_X[N], the least squares solution.
//
{
  double *x;

  if ( prob == 1 )
  {
    x = p01_x ( n );
  }
  else if ( prob == 2 )
  {
    x = p02_x ( n );
  }
  else if ( prob == 3 )
  {
    x = p03_x ( n );
  }
  else if ( prob == 4 )
  {
    x = p04_x ( n );
  }
  else if ( prob == 5 )
  {
    x = p05_x ( n );
  }
  else if ( prob == 6 )
  {
    x = p06_x ( n );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_X - Fatal error!\n";
    cerr << "  Illegal value of PROB = " << prob << "\n";
    exit ( 1 );
  }

  return x;
}
//****************************************************************************80

double *p01_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P01_A returns the matrix A for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P01_A[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    a[i+0*m] = 1.0;
    for ( j = 1; j < n; j++ )
    {
      a[i+j*m] = a[i+(j-1)*m] * ( double ) ( i + 1 );
    }
  }
  return a;
}
//****************************************************************************80

double *p01_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P01_B returns the right hand side B for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P01_B[M], the right hand side.
//
{
  double *b;
  double b_save[5] = { 1.0, 2.3, 4.6, 3.1, 1.2 };

  b = r8vec_copy_new ( m, b_save );

  return b;
}
//****************************************************************************80

int p01_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_M returns the number of equations M for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P01_M, the number of equations.
//
{
  int m;

  m = 5;

  return m;
}
//****************************************************************************80

int p01_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_N returns the number of variables N for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P01_N, the number of variables.
//
{
  int n;

  n = 3;

  return n;
}
//****************************************************************************80

double *p01_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P01_X returns the least squares solution X for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P01_X[N], the least squares solution.
//
{
  double *x;
  double x_save[3] = { -3.0200000, 4.4914286, -0.72857143 };

  x = r8vec_copy_new ( n, x_save );

  return x;
}
//****************************************************************************80

double *p02_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P02_A returns the matrix A for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    Numerical Computing with MATLAB,
//    SIAM, 2004,
//    ISBN13: 978-0-898716-60-3,
//    LC: QA297.M625,
//    ebook: http://www.mathworks.com/moler/chapters.html
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P02_A[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    a[i+(n-1)*m] = 1.0;
    for ( j = n - 2; 0 <= j; j-- )
    {
      a[i+j*m] = a[i+(j+1)*m] * ( double ) ( i ) / 5.0;
    }
  }
  return a;
}
//****************************************************************************80

double *p02_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P02_B returns the right hand side B for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P02_B[M], the right hand side.
//
{
  double *b;
  double b_save[6] = { 150.697, 179.323, 203.212, 226.505, 249.633, 281.422 };

  b = r8vec_copy_new ( m, b_save );

  return b;
}
//****************************************************************************80

int p02_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_M returns the number of equations M for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P02_M, the number of equations.
//
{
  int m;

  m = 6;

  return m;
}
//****************************************************************************80

int p02_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_N returns the number of variables N for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P02_N, the number of variables.
//
{
  int n;

  n = 3;

  return n;
}
//****************************************************************************80

double *p02_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P02_X returns the least squares solution X for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P02_X[N], the least squares solution.
//
{
  double *x;
  double x_save[3] = { 5.7013, 121.1341, 152.4745 };

  x = r8vec_copy_new ( n, x_save );

  return x;
}
//****************************************************************************80

double *p03_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P03_A returns the matrix A for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Cleve Moler,
//    Numerical Computing with MATLAB,
//    SIAM, 2004,
//    ISBN13: 978-0-898716-60-3,
//    LC: QA297.M625,
//    ebook: http://www.mathworks.com/moler/chapters.html
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P03_A[M*N], the matrix.
//
{
  double *a;
  double a_save[5*3] = {
    1.0, 4.0, 7.0, 10.0, 13.0, 
    2.0, 5.0, 8.0, 11.0, 14.0, 
    3.0, 6.0, 9.0, 12.0, 15.0 };

  a = r8mat_copy_new ( m, n, a_save );

  return a;
}
//****************************************************************************80

double *p03_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P03_B returns the right hand side B for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P03_B[M], the right hand side.
//
{
  double *b;
  double b_save[5] = { 16.0, 17.0, 18.0, 19.0, 20.0 };

  b = r8vec_copy_new ( m, b_save );

  return b;
}
//****************************************************************************80

int p03_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_M returns the number of equations M for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P03_M, the number of equations.
//
{
  int m;

  m = 5;

  return m;
}
//****************************************************************************80

int p03_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_N returns the number of variables N for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P03_N, the number of variables.
//
{
  int n;

  n = 3;

  return n;
}
//****************************************************************************80

double *p03_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P03_X returns the least squares solution X for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P03_X[N], the least squares solution.
//
{
  double *x;
  double x_save[3] = { -7.5555556, 0.1111111, 7.7777778 };

  x = r8vec_copy_new ( n, x_save );

  return x;
}
//****************************************************************************80

double *p04_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P04_A returns the matrix A for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P04_A[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( i4_power ( j + 1, i ) );
    }
  }

  return a;
}
//****************************************************************************80

double *p04_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P04_B returns the right hand side B for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P04_B[M], the right hand side.
//
{
  double *b;
  double b_save[3] = { 15.0, 55.0, 225.0 };

  b = r8vec_copy_new ( m, b_save );

  return b;
}
//****************************************************************************80

int p04_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_M returns the number of equations M for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P04_M, the number of equations.
//
{
  int m;

  m = 3;

  return m;
}
//****************************************************************************80

int p04_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_N returns the number of variables N for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P04_N, the number of variables.
//
{
  int n;

  n = 5;

  return n;
}
//****************************************************************************80

double *p04_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P04_X returns the least squares solution X for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P04_X[N], the least squares solution.
//
{
  double *x;
  double x_save[5] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

  x = r8vec_copy_new ( n, x_save );

  return x;
}
//****************************************************************************80

double *p05_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P05_A returns the matrix A for problem 5.
//
//  Discussion:
//
//    A is the Hilbert matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P05_A[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = 1.0 / ( double ) ( i + j + 1 );
    }
  }

  return a;
}
//****************************************************************************80

double *p05_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P05_B returns the right hand side B for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P05_B[M], the right hand side.
//
{
  double *b;
  int i;

  b = new double[m];

  b[0] = 1.0;
  for ( i = 1; i < m; i++ )
  {
    b[i] = 0.0;
  }

  return b;
}
//****************************************************************************80

int p05_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_M returns the number of equations M for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P05_M, the number of equations.
//
{
  int m;

  m = 10;

  return m;
}
//****************************************************************************80

int p05_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_N returns the number of variables N for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P05_N, the number of variables.
//
{
  int n;

  n = 10;

  return n;
}
//****************************************************************************80

double *p05_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P05_X returns the least squares solution X for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P05_X[N], the least squares solution.
//
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = r8_mop ( i + 2 ) * ( double ) ( i + 1 ) 
      * r8_choose ( n + i, n - 1 ) * r8_choose ( n, n - i - 1 );
  }

  return x;
}
//****************************************************************************80

double *p06_a ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    P06_A returns the matrix A for problem 6.
//
//  Discussion:
//
//    A is a symmetric, orthogonal matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Input, int N, the number of variables.
//
//    Output, double P06_A[M*N], the matrix.
//
{
  double *a;
  double angle;
  int i;
  int j;
  static double pi = 3.141592653589793;

  a = new double[m*n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      angle = ( double ) ( ( i + 1 ) * ( j + 1 ) ) * pi / ( double ) ( n + 1 );
      a[i+j*m] = sin ( angle ) * sqrt ( 2.0 / ( double ) ( n + 1 ) );
    }
  }
  return a;
}
//****************************************************************************80

double *p06_b ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    P06_B returns the right hand side B for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of equations.
//
//    Output, double P06_B[M], the right hand side.
//
{
  double *b;
  int i;

  b = new double[m];

  b[0] = 1.0;
  for ( i = 1; i < m; i++ )
  {
    b[i] = 0.0;
  }
  return b;
}
//****************************************************************************80

int p06_m ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_M returns the number of equations M for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P06_M, the number of equations.
//
{
  int m;

  m = 10;

  return m;
}
//****************************************************************************80

int p06_n ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_N returns the number of variables N for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P06_N, the number of variables.
//
{
  int n;

  n = 10;

  return n;
}
//****************************************************************************80

double *p06_x ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    P06_X returns the least squares solution X for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of variables.
//
//    Output, double P06_X[N], the least squares solution.
//
{
  double angle;
  int i;
  static double pi = 3.141592653589793;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = ( double ) ( i + 1 ) * pi / ( double ) ( n + 1 );
    x[i] = sin ( angle ) * sqrt ( 2.0 / ( double ) ( n + 1 ) );
  }
  return x;
}
