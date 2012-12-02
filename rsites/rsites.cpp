# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int get_seed ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RSITES.
//
//  Discussion:
//
//    RSITES returns a set of N random numbers in M dimensions.
//
//  Modified:
//
//    22 October 2003
//
//  Author:
//
//    Ken Clarkson, AT&T Research (original C source)
//
//  Parameters:
//
//    Commandline argument, int N, the number of points to generate.
//
//    Commandline argument, int M, the spatial dimension.
{
  int i;
  int j;
  int m;
  int n;
  unsigned int seed;

  if ( argc < 2 )
  {
    cerr << "\n";
    cerr << "RSITES:\n";
    cerr << "  C++ version\n";
    cerr << "\n";
    cerr << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cerr << "\n";
    cerr << "  Please enter N, the number of points.\n";

    cin >> n;
  }
  else
  {
    sscanf ( argv[1], "%d", &n );
  }

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "RSITES - Fatal error!\n";
    cerr << "  N must be at least 1.\n";
    exit ( 1 );
  }

  if ( argc < 3 )
  {
    cerr << "\n";
    cerr << "RSITES:\n";
    cerr << "  Please enter M, the spatial dimension.\n";

    cin >> m;
  }
  else
  {
    sscanf ( argv[2], "%d", &m );
  }

  if ( m < 1 )
  {
    cout << "\n";
    cout << "RSITES - Fatal error!\n";
    cout << "  M must be at least 1.\n";
    exit ( 1 );
  }

  if ( argc < 4 )
  {
    seed = get_seed ( );
  }
  else
  {
    sscanf ( argv[3], "%d", &seed );
  }

  srand ( seed );

  cout << "#  output file of N points in M dimensions,\n";
  cout << "#  created by RSITES.\n";
  cout << "#\n";
  cout << "#  Point coordinates are between 0 and 999999.\n";
  cout << "#\n";
  cout << "#  Number of points, N = " << n << "\n";
  cout << "#  Spatial dimension, M = " << m << "\n";
  cout << "#  Value of SEED for SRAND = " << seed << "\n";
  cout << "#\n";

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < m; j++ )
    {
      cout << setw(6) << ( rand ( ) % 1000000 ) << "  ";
    }
    cout << "\n";
  }
  return 0;
}
//****************************************************************************80

int get_seed ( void )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I_MAX 2147483647
  time_t clock;
  int i;
  int ihour;
  int imin;
  int isec;
  int seed;
  struct tm *lt;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  clock = time ( &tloc );
  lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
  ihour = lt->tm_hour;

  if ( 12 < ihour )
  {
    ihour = ihour - 12;
  }
//
//  Move Hours to 0, 1, ..., 11
//
  ihour = ihour - 1;

  imin = lt->tm_min;

  isec = lt->tm_sec;

  seed = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
  seed = seed + 1;
//
//  Remap ISEED from [1,43200] to [1,IMAX].
//
  seed = ( int )
    ( ( ( float ) seed )
    * ( ( float ) I_MAX ) / ( 60.0E+00 * 60.0E+00 * 12.0E+00 ) );
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;
# undef I_MAX
}
