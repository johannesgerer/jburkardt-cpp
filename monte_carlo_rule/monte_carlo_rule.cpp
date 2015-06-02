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
double *r8mat_uniform_01_new ( int m, int n, int &seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MONTE_CARLO_RULE.
//
//  Discussion:
//
//    MONTE_CARLO_RULE generates N points in the M-dimensional unit hypercube,
//    and writes out files so that the data can be regarded as a quadrature rule.
//
//  Usage:
//
//    monte_carlo_rule m n seed
//
//    where
//
//    * M, the spatial dimension,
//    * N, the number of points to generate,
//    * SEED, the seed, a positive integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  string filename_r;
  string filename_w;
  string filename_x;
  int i;
  int m;
  int n;
  double *r;
  int s;
  int seed;
  double *w;
  double *x;

  timestamp ( );
  cout << "\n";
  cout << "MONTE_CARLO_RULE\n";
  cout << "  C++ version\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Compute the abscissas and weights of a quadrature rule\n";
  cout << "  that is simply a Monte Carlo sampling.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * M, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * SEED, a positive integer.\n";
  cout << "\n";
  cout << "  Output from the program includes\n";
  cout << "  a set of 3 files that define the quadrature rule.\n";
  cout << "\n";
  cout << "    (1) \"mc_m?_n?_s?_r.txt\", the ranges;\n";
  cout << "    (2) \"mc_m?_n?_s?_w.txt\", the weights;\n";
  cout << "    (3) \"mc_m?_n?_s?_x.txt\", the abscissas.\n";
//
//  Get the spatial dimension M.
//
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the spatial dimension M (1 or greater)\n";
    cin >> m;
  }
//
//  Get the number of points N.
//
  if ( 2 < argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the number of points N (1 or greater):\n";
    cin >> n;
  }
//
//  Get the seed S.
//
  if ( 3 < argc )
  {
    s = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the seed S (1 or greater):\n";
    cin >> s;
  }
//
//  Input summary.
//
  cout << "\n";
  cout << "  M = " << m << "\n";
  cout << "  N = " << n << "\n";
  cout << "  S = " << s << "\n";
//
//  Construct the rule.
//
  r = new double[m*2];

  for ( i = 0; i < m; i++ )
  {
    r[i+0*m] = 0.0;
    r[i+1*m] = 1.0;
  }

  w = new double[n];

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0 / ( double ) n;
  }

  seed = s;
  x = r8mat_uniform_01_new ( m, n, seed );
//
//  Output the rule.
//
  filename_r = "mc_d" + i4_to_string ( m, "%d" ) 
    + "_n" + i4_to_string ( n, "%d" ) 
    + "_s" + i4_to_string ( s, "%d" ) + "_r.txt";

  filename_w = "mc_d" + i4_to_string ( m, "%d" ) 
    + "_n" + i4_to_string ( n, "%d" ) 
    + "_s" + i4_to_string ( s, "%d" ) + "_w.txt";

  filename_x = "mc_d" + i4_to_string ( m, "%d" ) 
    + "_n" + i4_to_string ( n, "%d" ) 
    + "_s" + i4_to_string ( s, "%d" ) + "_x.txt";

  cout << "\n";
  cout << "  Region file will be   \"" << filename_r << "\".\n";
  cout << "  Weight file will be   \"" << filename_w << "\".\n";
  cout << "  Abscissa file will be \"" << filename_x << "\".\n";

  r8mat_write ( filename_r, m, 2, r );
  r8mat_write ( filename_w, 1, n, w );
  r8mat_write ( filename_x, m, n, x );
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
  cout << "MONTE_CARLO_RULE:\n";
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

double *r8mat_uniform_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
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
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
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
