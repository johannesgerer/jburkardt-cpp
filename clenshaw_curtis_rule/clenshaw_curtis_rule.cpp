# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void clenshaw_curtis_compute ( int order, double xtab[], double weight[] );
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
//    MAIN is the main program for CLENSHAW_CURTIS_RULE.
//
//  Discussion:
//
//    This program computes a standard Clenshaw Curtis quadrature rule
//    and writes it to a file.
//
//    The user specifies:
//    * the ORDER (number of points) in the rule;
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
//    21 February 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  string filename;
  int order;
  double *r;
  double *w;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "CLENSHAW_CURTIS_RULE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute a Clenshaw Curtis rule for approximating\n";
  cout << "\n";
  cout << "    Integral ( -1 <= x <= +1 ) f(x) dx\n";
  cout << "\n";
  cout << "  of order ORDER.\n";
  cout << "\n";
  cout << "  The user specifies ORDER, A, B and FILENAME.\n";
  cout << "\n";
  cout << "  ORDER is the number of points.\n";
  cout << "\n";
  cout << "  A is the left endpoint.\n";
  cout << "\n";
  cout << "  B is the right endpoint.\n";
  cout << "\n";
  cout << "  FILENAME is used to generate 3 files:\n";
  cout << "\n";
  cout << "    filename_w.txt - the weight file\n";
  cout << "    filename_x.txt - the abscissa file.\n";
  cout << "    filename_r.txt - the region file.\n";
//
//  Get ORDER.
//
  if ( 1 < argc )
  {
    order = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of ORDER (1 or greater)\n";
    cin >> order;
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
    cout << "  Enter FILENAME, the \"root name\" of the quadrature files).\n";
    cin >> filename;
  }
//
//  Input summary.
//
  cout << "\n";
  cout << "  ORDER = " << order << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
  cout << "  FILENAME = \"" << filename << "\".\n";
//
//  Construct the rule.
//
  r = new double[2];
  w = new double[order];
  x = new double[order];

  r[0] = a;
  r[1] = b;

  clenshaw_curtis_compute ( order, x, w );
//
//  Rescale the rule.
//
  rescale ( a, b, order, x, w );
//
//  Output the rule.
//
  rule_write ( order, filename, x, w, r );
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
  cout << "CLENSHAW_CURTIS_RULE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void clenshaw_curtis_compute ( int order, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function is w(x) = 1.0.
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    The quadrature rule:
//
//      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//    1 <= ORDER.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  double b;
  int i;
  int j;
  double pi = 3.141592653589793;
  double theta;

  if ( order < 1 )
  {
    std::cerr << "\n";
    std::cerr << "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n";
    std::cerr << "  Illegal value of ORDER = " << order << "\n";
    std::exit ( 1 );
  }
  else if ( order == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      x[i] =  std::cos ( ( double ) ( order - 1 - i ) * pi
                       / ( double ) ( order - 1     ) );
    }
    x[0] = -1.0;
    if ( ( order % 2 ) == 1 )
    {
      x[(order-1)/2] = 0.0;
    }
    x[order-1] = +1.0;

    for ( i = 0; i < order; i++ )
    {
      theta = ( double ) ( i ) * pi / ( double ) ( order - 1 );

      w[i] = 1.0;

      for ( j = 1; j <= ( order - 1 ) / 2; j++ )
      {
        if ( 2 * j == ( order - 1 ) )
        {
          b = 1.0;
        }
        else
        {
          b = 2.0;
        }

        w[i] = w[i] - b *  std::cos ( 2.0 * ( double ) ( j ) * theta )
          / ( double ) ( 4 * j * j - 1 );
      }
    }

    w[0] = w[0] / ( double ) ( order - 1 );
    for ( i = 1; i < order - 1; i++ )
    {
      w[i] = 2.0 * w[i] / ( double ) ( order - 1 );
    }
    w[order-1] = w[order-1] / ( double ) ( order - 1 );
  }

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
//    RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
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
