# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
double *ccn_compute_points_new ( int n );
int i4_min ( int i1, int i2 );
double *nc_compute_new ( int n, double x_min, double x_max, double x[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void rescale ( double a, double b, int n, double x[], double w[] );
void rule_write ( int order, string filename, double x[], double w[],
  double r[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CCN_RULE.
//
//  Discussion:
//
//    This program computes a nested Clenshaw Curtis quadrature rule
//    and writes it to a file.
//
//    The user specifies:
//    * N, the number of points in the rule;
//    * A, the left endpoint;
//    * B, the right endpoint;
//    * FILENAME, which defines the output filenames.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  string filename;
  int n;
  double *r;
  double *w;
  double *x;
  double x_max;
  double x_min;

  timestamp ( );
  cout << "\n";
  cout << "CCN_RULE\n";
  cout << "  C++ version\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute one of a family of nested Clenshaw Curtis rules\n";
  cout << "  for approximating\n";
  cout << "    Integral ( -1 <= x <= +1 ) f(x) dx\n";
  cout << "  of order N.\n";
  cout << "\n";
  cout << "  The user specifies N, A, B and FILENAME.\n";
  cout << "\n";
  cout << "  N is the number of points.\n";
  cout << "  A is the left endpoint.\n";
  cout << "  B is the right endpoint.\n";
  cout << "  FILENAME is used to generate 3 files:\n";
  cout << "    filename_w.txt - the weight file\n";
  cout << "    filename_x.txt - the abscissa file.\n";
  cout << "    filename_r.txt - the region file.\n";
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
    cout << "  Enter the value of N (1 or greater)\n";
    cin >> n;
  }
//
//  Get A.
//
  if ( 2 < argc )
  {
    a = atof ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the left endpoint A:\n";
    cin >> a;
  }
//
//  Get B.
//
  if ( 3 < argc )
  {
    b = atof ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the right endpoint B:\n";
    cin >> b;
  }
//
//  Get FILENAME:
//
  if ( 4 < argc )
  {
    filename = argv[4];
  }
  else
  {
    cout << "\n";
    cout << "  Enter FILENAME, the \"root name\" of the quadrature files.\n";
    cin >> filename;
  }
//
//  Input summary.
//
  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  FILENAME = \"" << filename << "\".\n";
//
//  Construct the rule.
//
  r = new double[2];

  r[0] = a;
  r[1] = b;

  x = ccn_compute_points_new ( n );

  x_min = -1.0;
  x_max = +1.0;
  w = nc_compute_new ( n, x_min, x_max, x );
//
//  Rescale the rule.
//
  rescale ( a, b, n, x, w );
//
//  Output the rule.
//
  rule_write ( n, filename, x, w, r );
//
//  Free memory.
//
  delete [] r;
  delete [] w;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "CCN_RULE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double *ccn_compute_points_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CCN_COMPUTE_POINTS: compute Clenshaw Curtis Nested points.
//
//  Discussion:
//
//    We want to compute the following sequence:
//
//    1/2,
//    0, 1
//    1/4, 3/4
//    1/8, 3/8, 5/8, 7/8,
//    1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
//
//    But we would prefer that the numbers in each row be regrouped in pairs
//    that are symmetric about 1/2, with the number above 1/2 coming first.
//    Thus, the last row might become:
//    (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
//
//    Once we have our sequence, we apply the Chebyshev transformation
//    which maps [0,1] to [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements to compute.
//
//    Output, double CCN_COMPUTE_POINTS_NEW[N], the elements of the sequence.
//
{
  int d;
  int i;
  int k;
  int m;
  double pi = 3.141592653589793;
  int td;
  int tu;
  double *x;

  x = new double[n];
//
//  Handle first three entries specially.
//
  if ( 1 <= n )
  {
    x[0] = 0.5;
  }

  if ( 2 <= n )
  {
    x[1] = 1.0;
  }

  if ( 3 <= n )
  {
    x[2] = 0.0;
  }

  m = 3;
  d = 2;

  while ( m < n )
  {
    tu = d + 1;
    td = d - 1;

    k = i4_min ( d, n - m );

    for ( i = 1; i <= k; i++ )
    {
      if ( ( i % 2 ) == 1 )
      {
        x[m+i-1] = tu / 2.0 / ( double ) ( k );
        tu = tu + 2;
      }
      else
      {
        x[m+i-1] = td / 2.0 / ( double ) ( k );
        td = td - 2;
      }
    }
    m = m + k;
    d = d * 2;
  }
//
//  Apply the Chebyshev transformation.
//
  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( x[i] * pi );
  }
  x[0] = 0.0;

  if ( 2 <= n )
  {
    x[1] = -1.0;
  }

  if ( 3 <= n )
  {
    x[2] = +1.0;
  }

  return x;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

double *nc_compute_new ( int n, double x_min, double x_max, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    NC_COMPUTE_NEW computes a Newton-Cotes quadrature rule.
//
//  Discussion:
//
//    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
//    estimates
//
//      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
//
//    using N abscissas X and weights W:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) ).
//
//    For the CLOSED rule, the abscissas include the end points.
//    For the OPEN rule, the abscissas do not include the end points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 November 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order.
//
//    Input, double X_MIN, X_MAX, the endpoints of the interval.
//
//    Input, double X[N], the abscissas.
//
//    Output, double NC_COMPUTE_NEW[N], the weights.
//
{
  double *d;
  int i;
  int j;
  int k;
  double *w;
  double yvala;
  double yvalb;

  d = new double[n];
  w = new double[n];

  for ( i = 0; i < n; i++ )
  {
//
//  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
//  and zero at the other nodes.
//
    for ( j = 0; j < n; j++ )
    {
      d[j] = 0.0;
    }
    d[i] = 1.0;

    for ( j = 2; j <= n; j++ )
    {
      for ( k = j; k <= n; k++ )
      {
        d[n+j-k-1] = ( d[n+j-k-1-1] - d[n+j-k-1] ) / ( x[n+1-k-1] - x[n+j-k-1] );
      }
    }

    for ( j = 1; j <= n - 1; j++ )
    {
      for ( k = 1; k <= n - j; k++ )
      {
        d[n-k-1] = d[n-k-1] - x[n-k-j] * d[n-k];
      }
    }
//
//  Evaluate the antiderivative of the polynomial at the left and
//  right endpoints.
//
    yvala = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvala = yvala * x_min + d[j] / ( double ) ( j + 1 );
    }
    yvala = yvala * x_min;

    yvalb = d[n-1] / ( double ) ( n );
    for ( j = n - 2; 0 <= j; j-- )
    {
      yvalb = yvalb * x_max + d[j] / ( double ) ( j + 1 );
    }
    yvalb = yvalb * x_max;

    w[i] = yvalb - yvala;
  }

  delete [] d;

  return w;
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
//    John Burkardt.
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

void rule_write ( int order, string filename, double x[], double w[],
  double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_WRITE writes a quadrature rule to three files.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 February 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, double A, the left endpoint.
//
//    Input, double B, the right endpoint.
//
//    Input, string FILENAME, specifies the output filenames.
//    "filename_w.txt", "filename_x.txt", "filename_r.txt"
//    defining weights, abscissas, and region.
//
{
  string filename_r;
  string filename_w;
  string filename_x;
  int i;
  int kind;

  filename_w = filename + "_w.txt";
  filename_x = filename + "_x.txt";
  filename_r = filename + "_r.txt";

  cout << "\n";
  cout << "  Creating quadrature files.\n";
  cout << "\n";
  cout << "  Root file name is     \"" << filename   << "\".\n";
  cout << "\n";
  cout << "  Weight file will be   \"" << filename_w << "\".\n";
  cout << "  Abscissa file will be \"" << filename_x << "\".\n";
  cout << "  Region file will be   \"" << filename_r << "\".\n";

  r8mat_write ( filename_w, 1, order, w );
  r8mat_write ( filename_x, 1, order, x );
  r8mat_write ( filename_r, 1, 2,     r );

  return;
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
