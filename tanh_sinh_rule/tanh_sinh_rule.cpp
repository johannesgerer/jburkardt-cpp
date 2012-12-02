# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void tanh_sinh_compute ( int order, double h, double xtab[], double weight[] );
void tanh_sinh_handle ( int order, string prefix );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TANH_SIH_RULE.
//
//  Discussion:
//
//    This program computes a tanh-sinh quadrature rule
//    and writes it to a file.
//
//    The user specifies:
//    * the ORDER (number of points) in the rule
//    * the OUTPUT option.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  int order;
  string prefix;

  timestamp ( );
  cout << "\n";
  cout << "TANH_SINH_RULE\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute a tanh-sinh rule for approximating\n";
  cout << "\n";
  cout << "    Integral ( -1 <= x <= +1 ) f(x) dx\n";
  cout << "\n";
  cout << "  of order ORDER.\n";
  cout << "\n";
  cout << "  The user specifies ORDER and PREFIX.\n";
  cout << "\n";
  cout << "  PREFIX is used to name the 3 quadrature files:\n";
  cout << "\n";
  cout << "    prefix_w.txt - the weight file\n";
  cout << "    prefix_x.txt - the abscissa file.\n";
  cout << "    prefix_r.txt - the region file.\n";
//
//  Get the order.
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

  cout << "\n";
  cout << "  The requested order of the rule is = " << order << "\n";
//
//  Get the quadrature file root name:
//
  if ( 2 < argc )
  {
    prefix = argv[2];
  }
  else
  {
    cout << "\n";
    cout << "  Enter PREFIX, the \"root name\" of the quadrature files.\n";
    cin >> prefix;
  }

  cout << "\n";
  cout << "  PREFIX is \"" << prefix << "\".\n";
//
//  Construct the rule and output it.
//
  tanh_sinh_handle ( order, prefix );

  cout << "\n";
  cout << "TANH_SINH_RULE:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void tanh_sinh_compute ( int order, double h, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_SINH_COMPUTE computes a tanh-sinh quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the quadrature order.
//
//    Input, double H, the spacing.
//
//    Output, double X[ORDER], the abscissas.
//
//    Output, double W[ORDER], the weights.
//
{
  double ct;
  double ct2;
  int i;
  double pi = 3.141592653589793;
  double st;
  double t;
  double w_sum;

  for ( i = 0; i < order; i++ )
  {
    t = ( double ) ( 2 * i - order + 1 ) * h / 2.0;

    ct = cosh ( t );
    st = sinh ( t );
    ct2 = cosh ( 0.5 * pi * st );

    x[i] = tanh ( 0.5 * pi * st );

    w[i] = 0.5 * pi * h * ct / ct2 / ct2;
  }
//
//  Normalize the weights so that they sum to 2.0.
//
  w_sum = 0.0;
  for ( i = 0; i < order; i++ )
  {
    w_sum = w_sum + w[i];
  }
  for ( i = 0; i < order; i++ )
  {
    w[i] = 2.0 * w[i] / w_sum;
  }

  return;
}
//****************************************************************************80

void tanh_sinh_handle ( int order, string prefix )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_SINH_HANDLE computes the requested tanh-sinh rule and outputs it.
//
//  Discussion:
//
//    The prefix is used to define the output file names:
//
//      prefix + "_r.txt",
//      prefix + "_w.txt",
//      prefix + "_x.txt".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ORDER, the order of the rule.
//
//    Input, string PREFIX, the file prefix.
// 
{
  double h;
  int i;
  string output_r;
  string output_w;
  string output_x;
  double *r;
  double *w;
  double *x;

  r = new double[2];
  w = new double[order];
  x = new double[order];

  r[0] = - 1.0;
  r[1] = + 1.0;
//
//  This choice of H is only one of many.
//  For our choice, the family ORDER = 1, 3, 7, 15, 31, 63, ... is nested.
//
// h = 16.0 / ( double ) ( order + 1 );
// h =  8.0 / ( double ) ( order + 1 );
   h =  4.0 / ( double ) ( order + 1 );

  tanh_sinh_compute ( order, h, x, w );

  output_w = prefix + "_w.txt";
  output_x = prefix + "_x.txt";
  output_r = prefix + "_r.txt";

  cout << "\n";
  cout << "  Creating quadrature files.\n";
  cout << "\n";
  cout << "  Order = " << order << "\n";
  cout << "  Parameter H = " << h << "\n";
  cout << "\n";
  cout << "  Root file name is     \"" << prefix   << "\".\n";
  cout << "\n";
  cout << "  Weight file will be   \"" << output_w << "\".\n";
  cout << "  Abscissa file will be \"" << output_x << "\".\n";
  cout << "  Region file will be   \"" << output_r << "\".\n";
            
  r8mat_write ( output_w, 1, order, w );
  r8mat_write ( output_x, 1, order, x );
  r8mat_write ( output_r, 1, 2,     r );
  
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
