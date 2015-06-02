# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "filon.hpp"

//****************************************************************************80

double filon_fun_cos ( int n, double *f ( int n, double x[] ), double a, 
  double b, double t )

//****************************************************************************80
//
//  Purpose:
//
//    FILON_FUN_COS uses Filon's method on integrals with a cosine factor.
//
//  Discussion:
//
//    The integral to be approximated has the form:
//
//      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
//
//    where T is user specified.
//
//    The function is interpolated over each subinterval by
//    a parabolic arc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Chase, Lloyd Fosdick,
//    An Algorithm for Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 453-457.
//
//    Stephen Chase, Lloyd Fosdick,
//    Algorithm 353:
//    Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 457-458.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of data points.
//    N must be odd, and greater than 1.
//
//    Input, double *F ( int n, double x[] ), the function which evaluates the 
//    integrand.
//
//    Input, double A, B, the limits of integration.
//
//    Input, double T, the multiplier of the X argument of the cosine.
//
//    Output, double FILON_FUN_COS, the approximate value of the integral.
//
{
  double alpha;
  double beta;
  double c2n;
  double c2nm1;
  double cost;
  double *ftab;
  double gamma;
  double h;
  int i;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    cerr << "\n";
    cerr << "FILON_FUN_COS - Fatal error!\n";
    cerr << "  N < 2\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    cerr << "\n";
    cerr << "FILON_FUN_COS - Fatal error!\n";
    cerr << "  N must be odd.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.0 * sint * sint ) 
      / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
//
//  Tabulate the function.
//
  ftab = f ( n, x );

  c2n = 0.5 * ftab[0] * cos ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    c2n = c2n + ftab[i] * cos ( t * x[i] );
  }
  c2n = c2n + 0.5 * ftab[n-1] * cos ( t * x[n-1] );

  c2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    c2nm1 = c2nm1 + ftab[i] * cos ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[n-1] * sin ( t * x[n-1] )  
              - ftab[0]   * sin ( t * x[0] ) ) 
    + beta * c2n 
    + gamma * c2nm1 );

  delete [] ftab;
  delete [] x;

  return value;
}
//****************************************************************************80

double filon_tab_cos ( int n, double ftab[], double a, double b, double t )

//****************************************************************************80
//
//  Purpose:
//
//    FILON_TAB_COS uses Filon's method on integrals with a cosine factor.
//
//  Discussion:
//
//    The integral to be approximated has the form:
//
//      Integral ( A <= X <= B ) F(X) * COS(T*X) dX
//
//    where T is user specified.
//
//    The function is interpolated over each subinterval by
//    a parabolic arc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Chase, Lloyd Fosdick,
//    An Algorithm for Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 453-457.
//
//    Stephen Chase, Lloyd Fosdick,
//    Algorithm 353:
//    Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 457-458.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of data points.
//    N must be odd, and greater than 1.
//
//    Input, double FTAB[N], contains the value of the function
//    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
//
//    Input, double A, B, the limits of integration.
//
//    Input, double T, the multiplier of the X argument of the cosine.
//
//    Output, double FILON_TAB_COS, the approximate value of the integral.
//
{
  double alpha;
  double beta;
  double c2n;
  double c2nm1;
  double cost;
  double gamma;
  double h;
  int i;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    cerr << "\n";
    cerr << "FILON_TAB_COS - Fatal error!\n";
    cerr << "  N < 2\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    cerr << "\n";
    cerr << "FILON_TAB_COS - Fatal error!\n";
    cerr << "  N must be odd.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost - 2.0 * sint * sint ) 
      / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost 
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }

  c2n = + 0.5 * ftab[0] * cos ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    c2n = c2n + ftab[i] * cos ( t * x[i] );
  }
  c2n = c2n + 0.5 * ftab[n-1] * cos ( t * x[n-1] );

  c2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    c2nm1 = c2nm1 + ftab[i] * cos ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[n-1] * sin ( t * x[n-1] )  
              - ftab[0]   * sin ( t * x[0] ) ) 
    + beta * c2n 
    + gamma * c2nm1 );

  delete [] x;

  return value;
}
//****************************************************************************80

double filon_fun_sin ( int n, double *f ( int n, double x[] ), double a, 
  double b, double t )

//****************************************************************************80
//
//  Purpose:
//
//    FILON_FUN_SIN uses Filon's method on integrals with a sine factor.
//
//  Discussion:
//
//    The integral to be approximated has the form
//
//      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
//
//    where T is user specified.
//
//    The function is interpolated over each subinterval by
//    a parabolic arc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Chase, Lloyd Fosdick,
//    An Algorithm for Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 453-457.
//
//    Stephen Chase, Lloyd Fosdick,
//    Algorithm 353:
//    Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 457-458.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of data points, 
//    including the endpoints.  N must be odd, and greater than 1.
//
//    Input, external F, the subroutine which evaluates the integrand,
//    of the form subroutine F ( N, X, FX ).
//
//    Input, double A, B, the limits of integration.
//
//    Input, double T, multiplier of the X argument of the sine.
//
//    Output, double FILON_FUN_SIN, the approximate value of the integral.
//
{
  double alpha;
  double beta;
  double cost;
  double *ftab;
  double gamma;
  double h;
  int i;
  double s2n;
  double s2nm1;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    cerr << "\n";
    cerr << "FILON_FUN_SIN - Fatal error!\n";
    cerr << "  N < 2\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    cerr << "\n";
    cerr << "FILON_FUN_SIN - Fatal error!\n";
    cerr << "  N must be odd.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost 
      - 2.0 * sint * sint ) / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
//
//  Tabulate the function.
//
  ftab = f ( n, x );

  s2n = + 0.5 * ftab[0] * sin ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    s2n = s2n + ftab[i] * sin ( t * x[i] );
  }
  s2n = s2n + 0.5 * ftab[n-1] * sin ( t * x[n-1] );

  s2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    s2nm1 = s2nm1 + ftab[i] * sin ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[0]   * cos ( t * x[0] ) 
              - ftab[n-1] * cos ( t * x[n-1] ) )
    + beta * s2n 
    + gamma * s2nm1 );
 
  delete [] ftab;
  delete [] x;

  return value;
}
//****************************************************************************80

double filon_tab_sin ( int n, double ftab[], double a, double b, double t )

//****************************************************************************80
//
//  Purpose:
//
//    FILON_TAB_SIN uses Filon's method on integrals with a sine factor.
//
//  Discussion:
//
//    The integral to be approximated has the form
//
//      Integral ( A <= X <= B ) F(X) * SIN(T*X) dX
//
//    where T is user specified.
//
//    The function is interpolated over each subinterval by
//    a parabolic arc.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Chase, Lloyd Fosdick,
//    An Algorithm for Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 453-457.
//
//    Stephen Chase, Lloyd Fosdick,
//    Algorithm 353:
//    Filon Quadrature,
//    Communications of the Association for Computing Machinery,
//    Volume 12, Number 8, August 1969, pages 457-458.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, int N, the number of data points, 
//    including the endpoints.  N must be odd, and greater than 1.
//
//    Input, double FTAB[N], contains the value of the function
//    at A, A+H, A+2*H, ... , B-H, B, where H = (B-A)/(N-1).
//
//    Input, double A, B, the limits of integration.
//
//    Input, double T, multiplier of the X argument of the sine.
//
//    Output, double FILON_TAB_SIN, the approximate value of the integral.
//
{
  double alpha;
  double beta;
  double cost;
  double gamma;
  double h;
  int i;
  double s2n;
  double s2nm1;
  double sint;
  double theta;
  double value;
  double *x;

  if ( a == b )
  {
    value = 0.0;
    return value;
  }
 
  if ( n <= 1 )
  {
    cerr << "\n";
    cerr << "FILON_TAB_SIN - Fatal error!\n";
    cerr << "  N < 2\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
 
  if ( ( n % 2 ) != 1 )
  {
    cerr << "\n";
    cerr << "FILON_TAB_SIN - Fatal error!\n";
    cerr << "  N must be odd.\n";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }
//
//  Set the X values.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * a   
           + ( double ) (     i     ) * b ) 
           / ( double ) ( n     - 1 );
  }

  h = ( b - a ) / ( double ) ( n - 1 );
  theta = t * h;
  sint = sin ( theta );
  cost = cos ( theta );

  if ( 6.0 * fabs ( theta ) <= 1.0 )
  {
    alpha = 2.0 * pow ( theta, 3 ) /   45.0 
          - 2.0 * pow ( theta, 5 ) /  315.0 
          + 2.0 * pow ( theta, 7 ) / 4725.0;
  
    beta =  2.0                    /     3.0 
          + 2.0 * pow ( theta, 2 ) /    15.0 
          - 4.0 * pow ( theta, 4 ) /   105.0 
          + 2.0 * pow ( theta, 6 ) /   567.0 
          - 4.0 * pow ( theta, 8 ) / 22275.0;

    gamma = 4.0                    /      3.0 
          - 2.0 * pow ( theta, 2 ) /     15.0 
          +       pow ( theta, 4 ) /    210.0 
          -       pow ( theta, 6 ) /  11340.0;
  }
  else
  {
    alpha = ( pow ( theta, 2 ) + theta * sint * cost 
      - 2.0 * sint * sint ) / pow ( theta, 3 );

    beta = ( 2.0 * theta + 2.0 * theta * cost * cost 
      - 4.0 * sint * cost ) / pow ( theta, 3 );

    gamma = 4.0 * ( sint - theta * cost ) / pow ( theta, 3 );
  }
  
  s2n = + 0.5 * ftab[0] * sin ( t * x[0] );
  for ( i = 2; i < n - 1; i = i + 2 )
  {
    s2n = s2n + ftab[i] * sin ( t * x[i] );
  }
  s2n = s2n + 0.5 * ftab[n-1] * sin ( t * x[n-1] );

  s2nm1 = 0.0;
  for ( i = 1; i <= n - 2; i = i + 2 )
  {
    s2nm1 = s2nm1 + ftab[i] * sin ( t * x[i] );
  }

  value = h * ( 
      alpha * ( ftab[0]   * cos ( t * x[0] ) 
              - ftab[n-1] * cos ( t * x[n-1] ) ) 
    + beta * s2n 
    + gamma * s2nm1 );
 
  delete [] x;

  return value;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
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
