# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>
# include <cstring>
# include <sstream>

using namespace std;

int main ( int argc, char *argv[] );
int get_seed ( );
int i4_uniform_ab ( int ilo, int ihi, int &seed );
double *latin_random_new ( int dim_num, int point_num, int &seed );
int *perm_uniform ( int n, int &seed );
double *r8mat_uniform_01_new ( int m, int n, int &seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LATIN_RANDOM_DATASET.
//
//  Discussion:
//
//    LATIN_RANDOM_DATASET generates a Latin Random Square dataset 
//    and writes it to a file.
//
//  Usage:
//
//    latin_random_dataset m n seed
//
//    where
//
//    * M, the spatial dimension,
//    * N, the number of points to generate,
//    * SEED, the seed, a positive integer.
//
//    creates an M by N dataset and writes it to the
//    file "latin_random_M_N.txt".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 November 2014
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
  cout << "LATIN_RANDOM_DATASET\n";
  cout << "  C++ version\n";
  cout << "  Generate a Latin Random Square dataset.\n";
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

  if ( seed == 0 )
  {
    seed = get_seed ( );
  }
//
//  Compute the data.
//
  r = latin_random_new ( m, n, seed );
//
//  Write it to a file.
//
  m_ostring << m;
  n_ostring << n;

  output_filename = "latin_random_" + m_ostring.str ( ) + "_" 
    + n_ostring.str ( ) + ".txt";

  r8mat_write ( output_filename, m, n, r );

  cout << "\n";
  cout << "  The data was written to the file \"" 
    << output_filename << "\".\n";
//
//  Free memory.
//
  delete [] r;
//
//  Terminate.
//
  cout << "\n";
  cout << "LATIN_RANDOM_DATASET:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
//  Remap SEED from [1,43200] to [1,IMAX].
//
  seed = ( int ) 
    ( ( ( double ) seed )
    * ( ( double ) I_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
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
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  const int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

double *latin_random_new ( int dim_num, int point_num, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    LATIN_RANDOM_NEW returns points in a Latin Random square.
//
//  Discussion:
//
//    In each spatial dimension, there will be exactly one
//    point whose coordinate value lies between consecutive
//    values in the list:
//
//      ( 0, 1, 2, ..., point_num ) / point_num
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input/output, int &SEED, a seed for UNIFORM.
//
//    Output, double LATIN_RANDOM_NEW[DIM_NUM,POINT_NUM], the points.
//
{
  int i;
  int j;
  int *perm;
  double r;
  double *x;

  x = r8mat_uniform_01_new ( dim_num, point_num, seed );
//
//  For spatial dimension I, 
//    pick a random permutation of 1 to POINT_NUM,
//    force the corresponding I-th components of X to lie in the
//    interval ( PERM[J]-1, PERM[J] ) / POINT_NUM.
//
  for ( i = 0; i < dim_num; i++ )
  {
    perm = perm_uniform_new ( point_num, seed );

    for ( j = 0; j < point_num; j++ )
    {
      x[i+j*dim_num] = ( ( ( double ) perm[j] ) + x[i+j*dim_num] ) 
                       / ( ( double ) point_num );
    }
    delete [] perm;
  }
  return x;
}
//****************************************************************************80

int *perm_uniform_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_UNIFORM_NEW selects a random permutation of N objects.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of objects to be permuted.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int PERM_UNIFORM_NEW[N], a permutation of
//    (0, 1, ..., N-1).
//
{
  int i;
  int j;
  int k;
  int *p;

  p = new int[n];

  for ( i = 0; i < n; i++ )
  {
    p[i] = i;
  }

  for ( i = 0; i < n - 1; i++ )
  {
    j = i4_uniform_ab ( i, n - 1, seed );
    k    = p[i];
    p[i] = p[j];
    p[j] = k;
  }

  return p;
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
  const int i4_huge = 2147483647;
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
        seed = seed + i4_huge;
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
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
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
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
      output << "  " << table[i+j*m];
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
