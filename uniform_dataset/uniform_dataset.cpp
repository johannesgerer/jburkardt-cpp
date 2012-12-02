# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for UNIFORM_DATASET.
//
//  Discussion:
//
//    UNIFORM_DATASET generates a uniform pseudorandom dataset 
//    and writes it to a file.
//
//  Usage:
//
//    uniform_dataset m n seed
//
//    where
//
//    * M, the spatial dimension,
//    * N, the number of points to generate,
//    * SEED, the seed, a positive integer.
//
//    creates an M by N uniform random dataset and writes it to the
//    file "uniform_M_N.txt".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  ostringstream m_ostring;
  int n;
  ostringstream n_ostring;
  string output_filename;
  double *r;
  int seed;

  timestamp ( );

  cout << "\n";
  cout << "UNIFORM_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a uniform pseudorandom dataset.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * M, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * SEED, a positive integer.\n";
  cout << "\n";
  cout << "  The program generates the data, writes it to the file\n";
  cout << "\n";
  cout << "    uniform_M_N.txt\n";
  cout << "\n";
  cout << "  where ""M"" and ""N"" are the numeric values specified\n";
  cout << "  by the user.\n";
//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of M\n";
    cin >> m;
  }

  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
//
//  Get the number of points.
//
  if ( 2 < argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the number of points N\n";
    cin >> n;
  }

  cout << "  Number of points N = " << n << "\n";
//
//  Get the seed.
//
  if ( 3 < argc )
  {
    seed = atoi ( argv[3] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter the value of SEED\n";
    cin >> seed;
  }

  cout << "  The seed is = " << seed << "\n";
//
//  Compute the data.
//
  r = r8mat_uniform_01_new ( m, n, &seed );
//
//  Write it to a file.
//
  m_ostring << m;
  n_ostring << n;

  output_filename = "uniform_" + m_ostring.str ( ) + "_" 
    + n_ostring.str ( ) + ".txt";

  r8mat_write ( output_filename, m, n, r );

  cout << "\n";
  cout << "  The data was written to the file \"" 
    << output_filename << "\".\n";
//
//  Terminate.
//
  delete [] r;

  cout << "\n";
  cout << "UNIFORM_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a new unit pseudorandom R8MAT.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
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
