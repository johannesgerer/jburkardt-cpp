# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
string i4_to_string ( int i4, string format );
void legendre_compute_glr ( int n, double x[], double w[] );
void legendre_compute_glr0 ( int n, double *p, double *pp );
void legendre_compute_glr1 ( int n, double *roots, double *ders );
void legendre_compute_glr2 ( double p, int n, double *roots, double *ders );
void legendre_handle ( int n, double a, double b );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
double rk2_leg ( double t, double tn, double x, int n );
void timestamp ( );
double ts_mult ( double *u, double h, int n );
double wtime ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LEGENDRE_RULE.
//
//  Discussion:
//
//    This program computes a standard Gauss-Legendre quadrature rule
//    and writes it to a file.
//
//  Usage:
//
//    legendre_rule_fast n a b
//
//    where
//
//    * n is the number of points in the rule;
//    * a is the left endpoint;
//    * b is the right endpoint.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n;

  timestamp ( );
  cout << "\n";
  cout << "LEGENDRE_RULE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute a Gauss-Legendre rule for approximating\n";
  cout << "\n";
  cout << "    Integral ( a <= x <= b ) f(x) dx\n";
  cout << "\n";
  cout << "  of order N.\n";
  cout << "\n";
  cout << "  The computed rule is written to 3 files:\n";
  cout << "\n";
  cout << "  * leg_oN_w.txt - the weight file\n";
  cout << "  * leg_oN_x.txt - the abscissa file.\n";
  cout << "  * leg_oN_r.txt - the region file.\n";
//
//  Get N.
//
  if ( 1 < argc )
  {
    n = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter N:\n";
    cin >> n;
  }

  cout << "\n";
  cout << "  N = " << n << "\n";
//
//  Get A:
//
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter A.\n";
    cin >> a;
  }

  cout << "\n";
  cout << "  A = " << a << "\n";
//
//  Get B:
//
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter B.\n";
    cin >> b;
  }

  cout << "\n";
  cout << "  B = " << b << "\n";
//
//  Construct the rule and output it.
//
  legendre_handle ( n, a, b );

  cout << "\n";
  cout << "LEGENDRE_RULE_FAST:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

string i4_to_string ( int i4, string format )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_STRING converts an I4 to a C++ string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I4, an integer.
//
//    Input, string FORMAT, the format string.
//
//    Output, string I4_TO_STRING, the string.
//
{
  char i4_char[80];
  string i4_string;

  sprintf ( i4_char, format.c_str ( ), i4 );

  i4_string = string ( i4_char );

  return i4_string;
}
//****************************************************************************80

void legendre_compute_glr ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_GLR: Legendre quadrature by the Glaser-Liu-Rokhlin method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double p;
  double pp;
  double w_sum;
//
//  Get the value and derivative of the N-th Legendre polynomial at 0.
//
  legendre_compute_glr0 ( n, &p, &pp );
//
//  If N is odd, then zero is a root.
//  
  if ( n % 2 == 1 )
  {
    x[(n-1)/2] = p;
    w[(n-1)/2] = pp;
  }
//
//  If N is even, we have to call a function to find the first root.
//
  else
  {
    legendre_compute_glr2 ( p, n, &x[n/2], &w[n/2] );
  }
//
//  Get the complete set of roots and derivatives.
//
  legendre_compute_glr1 ( n, x, w );
//
//  Compute the W.
//
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( 1.0 - x[i] ) / ( 1.0 + x[i] ) / w[i] / w[i];
  }
  w_sum = 0.0;
  for ( i = 0; i < n; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }
  return;
}
//****************************************************************************80

void legendre_compute_glr0 ( int n, double *p, double *pp )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_GLR0 gets a starting value for the fast algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, int N, the order of the Legendre polynomial.
//
//    Output, double *P, *PP, the value of the N-th Legendre polynomial
//    and its derivative at 0.
//
{
  double dk;
  int k;
  double pm1;
  double pm2;
  double ppm1;
  double ppm2;

  pm2 = 0.0;
  pm1 = 1.0;
  ppm2 = 0.0;
  ppm1 = 0.0;

  for ( k = 0; k < n; k++)
  {
    dk = ( double ) k;
    *p = - dk * pm2 / ( dk + 1.0 );
    *pp = ( ( 2.0 * dk + 1.0 ) * pm1 - dk * ppm2 ) / ( dk + 1.0 );
    pm2 = pm1;
    pm1 = *p;
    ppm2 = ppm1;
    ppm1 = *pp;
  }
  return;
}
//****************************************************************************80

void legendre_compute_glr1 ( int n, double *x, double *w )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_GLR1 gets the complete set of Legendre points and weights.
//
//  Discussion:
//
//    This routine requires that a starting estimate be provided for one
//    root and its derivative.  This information will be stored in entry
//    (N+1)/2 if N is odd, or N/2 if N is even, of X and W.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, int N, the order of the Legendre polynomial.
//
//    Input/output, double X[N].  On input, a starting value
//    has been set in one entry.  On output, the roots of the Legendre 
//    polynomial.
//
//    Input/output, double W[N].  On input, a starting value
//    has been set in one entry.  On output, the derivatives of the Legendre 
//    polynomial at the zeros.
//
//  Local Parameters:
//
//    Local, int M, the number of terms in the Taylor expansion.
//
{
  double dk;
  double dn;
  double h;
  int j;
  int k;
  int l;
  int m = 30;
  int n2;
  static double pi = 3.141592653589793;
  int s;
  double *u;
  double *up;
  double xp;

  if ( n % 2 == 1 )
  {
    n2 = ( n - 1 ) / 2 - 1;
    s = 1;
  }
  else
  {
    n2 = n / 2 - 1;
    s = 0;
  }

  u = new double[m+2];
  up = new double[m+1];

  dn = ( double ) n;

  for ( j = n2 + 1; j < n - 1; j++ )
  {
    xp = x[j];

    h = rk2_leg ( pi/2.0, -pi/2.0, xp, n ) - xp;

    u[0] = 0.0;
    u[1] = 0.0;
    u[2] = w[j];

    up[0] = 0.0;
    up[1] = u[2];

    for ( k = 0; k <= m - 2; k++ )
    {
      dk = ( double ) k;

      u[k+3] = 
      ( 
        2.0 * xp * ( dk + 1.0 ) * u[k+2]
        + ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1] / ( dk + 1.0 ) 
      ) / ( 1.0 - xp ) / ( 1.0 + xp ) / ( dk + 2.0 );

      up[k+2] = ( dk + 2.0 ) * u[k+3];
    }

    for ( l = 0; l < 5; l++ )
    { 
      h = h - ts_mult ( u, h, m ) / ts_mult ( up, h, m-1 );
    }

    x[j+1] = xp + h;
    w[j+1] = ts_mult ( up, h, m - 1 );    
  }

  for ( k = 0; k <= n2 + s; k++ )
  {
    x[k] = - x[n-1-k];
    w[k] = w[n-1-k];
  }
  return;
}
//****************************************************************************80

void legendre_compute_glr2 ( double pn0, int n, double *x1, double *d1 )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_COMPUTE_GLR2 finds the first real root.
//
//  Discussion:
//
//    This function is only called if N is even.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, double PN0, the value of the N-th Legendre polynomial
//    at 0.
//
//    Input, int N, the order of the Legendre polynomial.
//
//    Output, double *X1, the first real root.
//
//    Output, double *D1, the derivative at X1.
//
//  Local Parameters:
//
//    Local, int M, the number of terms in the Taylor expansion.
//
{
  double dk;
  double dn;
  int k;
  int l;
  int m = 30;
  static double pi = 3.141592653589793;
  double t;
  double *u;
  double *up;

  t = 0.0;
  *x1 = rk2_leg ( t, -pi/2.0, 0.0, n );

  u = new double[m+2];
  up = new double[m+1];

  dn = ( double ) n;
//
//  U[0] and UP[0] are never used.
//  U[M+1] is set, but not used, and UP[M] is set and not used.
//  What gives?
//
  u[0] = 0.0;
  u[1] = pn0;

  up[0] = 0.0;

  for ( k = 0; k <= m - 2; k = k + 2 )
  {
    dk = ( double ) k;

    u[k+2] = 0.0;
    u[k+3] = ( dk * ( dk + 1.0 ) - dn * ( dn + 1.0 ) ) * u[k+1]
      / (dk + 1.0) / (dk + 2.0 ); 

    up[k+1] = 0.0;
    up[k+2] = ( dk + 2.0 ) * u[k+3];
  }
  
  for ( l = 0; l < 5; l++ )
  {
    *x1 = *x1 - ts_mult ( u, *x1, m ) / ts_mult ( up, *x1, m-1 );
  }
  *d1 = ts_mult ( up, *x1, m-1 );

  return;
}
//****************************************************************************80

void legendre_handle ( int n, double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_HANDLE computes the requested Gauss-Legendre rule and outputs it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Input, double A, B, the left and right endpoints of the integration
//    interval.
// 
{
  string output_r;
  string output_w;
  string output_x;
  double *r;
  double t;
  double *w;
  double *x;

  r = new double[2];
  w = new double[n];
  x = new double[n];

  r[0] = a;
  r[1] = b;
//
//  Compute the rule.
//
  t = wtime ( );
  legendre_compute_glr ( n, x, w );
  t = wtime ( ) - t;
  cout << "\n";
  cout << "  Elapsed time during computation was " << t << " seconds.\n";
//
//  Rescale the rule.
//
  rescale ( a, b, n, x, w );
//
//  Write the data to files.
//
  output_w = "leg_o" + i4_to_string ( n, "%d" ) + "_w.txt";
  output_x = "leg_o" + i4_to_string ( n, "%d" ) + "_x.txt";
  output_r = "leg_o" + i4_to_string ( n, "%d" ) + "_r.txt";

  cout << "\n";
  cout << "  Weight file will be   \"" << output_w << "\".\n";
  cout << "  Abscissa file will be \"" << output_x << "\".\n";
  cout << "  Region file will be   \"" << output_r << "\".\n";
            
  r8mat_write ( output_w, 1, n, w );
  r8mat_write ( output_x, 1, n, x );
  r8mat_write ( output_r, 1, 2, r );
  
  delete [] r;
  delete [] w;
  delete [] x;

  return;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void rescale ( double a, double b, int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESCALE rescales a Legendre quadrature rule from [-1,+1] to [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2009
//
//  Author:
//
//    Original MATLAB version by Nick Hale.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Andreas Glaser, Xiangtao Liu, Vladimir Rokhlin, 
//    A fast algorithm for the calculation of the roots of special functions, 
//    SIAM Journal on Scientific Computing,
//    Volume 29, Number 4, pages 1420-1438, 2007.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the new interval.
//
//    Input, int N, the order.
//
//    Input/output, double X[N], on input, the abscissas for [-1,+1].
//    On output, the abscissas for [A,B].
//
//    Input/output, double W[N], on input, the weights for [-1,+1].
//    On output, the weights for [A,B].
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( a + b ) + ( b - a ) * x[i] ) / 2.0;
  }
  for ( i = 0; i < n; i++ )
  {
    w[i] = ( b - a ) * w[i] / 2.0;
  }
  return;
}
//****************************************************************************80   

double rk2_leg ( double t1, double t2, double x, int n )

//****************************************************************************80
//
//  Purpose:
//
//    RK2_LEG advances the value of X(T) using a Runge-Kutta method.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 October 2009
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double T1, T2, the range of the integration interval.
//
//    Input, double X, the value of X at T1.
//
//    Input, int N, the number of steps to take.
//
//    Output, double RK2_LEG, the value of X at T2.
//
{
  double f;
  double h;
  int j;
  double k1;
  double k2;
  int m = 10;
  double snn1;
  double t;

  h = ( t2 - t1 ) / ( double ) m;
  snn1 = sqrt ( ( double ) ( n * ( n + 1 ) ) );
  t = t1;

  for ( j = 0; j < m; j++ )
  {
    f = ( 1.0 - x ) * ( 1.0 + x );
    k1 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + k1;

    t = t + h;

    f = ( 1.0 - x ) * ( 1.0 + x );
    k2 = - h * f / ( snn1 * sqrt ( f ) - 0.5 * x * sin ( 2.0 * t ) );
    x = x + 0.5 * ( k2 - k1 );
  }
  return x;
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
//****************************************************************************80

double ts_mult ( double *u, double h, int n )

//****************************************************************************80
//
//  Purpose:
//
//    TS_MULT evaluates a polynomial.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2013
//
//  Author:
//
//    Original C++ version by Nick Hale.
//    This C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double U[N+1], the polynomial coefficients.
//    U[0] is ignored.
//
//    Input, double H, the polynomial argument.
//
//    Input, int N, the number of terms to compute.
//
//    Output, double TS_MULT, the value of the polynomial.
//
{
  double hk;
  int k;
  double ts;
  
  ts = 0.0;
  hk = 1.0;
  for ( k = 1; k<= n; k++ )
  {
    ts = ts + u[k] * hk;
    hk = hk * h;
  }
  return ts;
}
//****************************************************************************80

double wtime ( )

//****************************************************************************80
//
//  Purpose:
//
//    WTIME estimates the elapsed wall clock time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double WTIME, the current elapsed wall clock time.
//
{
  double now;

  now = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC; 

  return now;
}  

