# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( );
void filename_inc ( string *filename );
int *life_init ( double prob, int m, int n, int &seed );
void life_update ( int m, int n, int grid[] );
void life_write ( string output_filename, int m, int n, int grid[] );
double r8_uniform_01 ( int &seed );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LIFE_SERIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Martin Gardner,
//    Mathematical Games:
//    The Fantastic Combinations of John Conway's new solitaire game "Life",
//    Scientific American,
//    Volume 223, Number 4, October 1970, pages 120-123.
//
{
  string filename = "life_000.txt";
  int it;
  int it_max;
  int m;
  int n;
  int *grid;
  double prob;
  int seed;

  timestamp ( );
  cout << "\n";
  cout << "LIFE_SERIAL\n";
  cout << "  C++ version\n";
  cout << "  Carry out a few steps of John Conway's\n";
  cout << "  Game of Life.\n";
  cout << "\n";

  it_max = 10;
  m = 10;
  n = 10;
  prob = 0.20;
  seed = 123456789;

  for ( it = 0; it <= it_max; it++ )
  {
    if ( it == 0 )
    {
      grid = life_init ( prob, m, n, seed );
    }
    else
    {
      life_update ( m, n, grid );
    }
    life_write ( filename, m, n, grid );
    cout << "  " << filename << "\n";
    filename_inc ( &filename );
  }
//
//  Free memory.
//
  delete [] grid;
//
//  Terminate.
//
  cout << "\n";
  cout << "LIFE_SERIAL\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

int *life_init ( double prob, int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    LIFE_INIT initializes the life grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double PROB, the probability that a grid cell
//    should be alive.
//
//    Input, int M, N, the number of rows and columns
//    of interior grid cells.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, int LIFE_INIT[(1+M+1)*(1+N+1)], the initial grid.
//
{
  int *grid;
  int i;
  int j;
  double r;

  grid = new int[ ( m + 2 ) * ( n + 2 )];
  for ( j = 0; j <= n + 1; j++ )
  {
    for ( i = 0; i <= m + 1; i++ )
    {
      grid[i+j*(m+2)] = 0;
    }
  }

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      r = r8_uniform_01 ( seed );
      if ( r <= prob )
      {
        grid[i+j*(m+2)] = 1;
      }
    }
  }

  return grid;
}
//****************************************************************************80

void life_update ( int m, int n, int grid[] )

//****************************************************************************80
//
//  Purpose:
//
//    LIFE_UPDATE updates a Life grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns
//    of interior grid cells.
//
//    Input/output, int GRID[(1+M+1)*(1+N+1)], the data.
//
{
  int i;
  int j;
  int *s;

  s = new int[m*n];

  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      s[i-1+(j-1)*m] = 
          grid[i-1+(j-1)*(m+2)] + grid[i-1+j*(m+2)] + grid[i-1+(j+1)*(m+2)]
        + grid[i  +(j-1)*(m+2)]                     + grid[i  +(j+1)*(m+2)]
        + grid[i+1+(j-1)*(m+2)] + grid[i+1+j*(m+2)] + grid[i+1+(j+1)*(m+2)];
    }
  }
//
//  Any dead cell with 3 live neighbors becomes alive.
//  Any living cell with less than 2 or more than 3 neighbors dies.
//
  for ( j = 1; j <= n; j++ )
  {
    for ( i = 1; i <= m; i++ )
    {
      if ( grid[i+j*(m+2)] == 0 )
      {
        if ( s[i-1+(j-1)*m] == 3 )
        {
          grid[i+j*(m+2)] = 1;
        }
      }
      else if ( grid[i+j*(m+2)] == 1 )
      {
        if ( s[i-1+(j-1)*m] < 2 || 3 < s[i-1+(j-1)*m] )
        {
          grid[i+j*(m+2)] = 0;
        }
      }
    }
  }

  delete [] s;

  return;
}
//****************************************************************************80

void life_write ( string output_filename, int m, int n, int grid[] )

//****************************************************************************80
//
//  Purpose:
//
//    LIFE_WRITE writes a grid to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output file name.
//
//    Input, int M, N, the number of rows and columns
//    of interior grid cells.
//
//    Input, int GRID[(1+M+1)*(1+N+1)], the data.
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
    cerr << "LIFE_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j <= n + 1; j++ )
  {
    for ( i = 0; i <= m + 1; i++ )
    {
      output << " " << grid[i+j*(m+2)];
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

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
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
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
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
