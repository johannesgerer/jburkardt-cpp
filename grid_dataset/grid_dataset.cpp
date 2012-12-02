# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <sstream>
# include <string>

using namespace std;

# include "grid.H"

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for GRID_DATASET.
//
//  Discussion:
//
//    GRID_DATASET generates a grid dataset and writes it to a file.
//
//    Interesting features of this problem are the determination
//    of the side of a grid that will generate "about" N points,
//    the method of dropping the extra points at random, and the
//    ability to center the grid inside the unit hypercube in a
//    number of ways.
//
//  Usage:
//
//    grid_dataset ( m, n, seed, center )
//
//    where
//
//    * M, the spatial dimension,
//    * N, the number of points to generate,
//    * SEED, the seed, a positive integer.
//    * CENTER, the grid centering option.
//      1: 0/(  N-1) ... (  N-1)/(  N-1)
//      2: 1/(  N+1) ...    N   /(  N+1)
//      3: 0/   N    ... (  N-1)/   N
//      4: 1/   N    ...    N   /   N
//      5: 1/(2*N)   ... (2*N-1)/(2*N  )
//
//    The program generates the data, writes it to the file
//
//      grid_M_N_CENTER.txt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 December 2009
//
//  Author:
//
//    John Burkardt
//
{
  int center;
  ostringstream center_ostring;
  char file_out_name[255];
  int m;
  ostringstream m_ostring;
  int n;
  ostringstream n_ostring;
  string output_filename;
  double *r;
  int seed;
  char *string;

  timestamp ( );

  cout << "\n";
  cout << "GRID_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a grid dataset.\n";

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
    cout << "  Recomputed SEED = " << seed << "\n";
  }
//
//  Get CENTER.
//
  if ( 4 < argc )
  {
    center = atoi ( argv[4] );
  }
  else
  {
    cout << "\n";
    cout << "  Enter CENTER, the grid centering option.\n";
    cout << "  Normal values are between 1 and 5:\n";
    cout << "  1: 0/(  N-1) ... (  N-1)/(  N-1)\n";
    cout << "  2: 1/(  N+1) ...    N   /(  N+1)\n";
    cout << "  3: 0/   N    ... (  N-1)/   N\n";
    cout << "  4: 1/   N    ...    N   /   N\n";
    cout << "  5: 1/(2*N)   ... (2*N-1)/(2*N  )\n";

    cin >> center;
  }
 
  cout << "  CENTER = " << center << "\n";
//
//  Compute the data.
//
  r = new double[m*n];

  grid_generate ( m, n, center, &seed, r );
//
//  Write the data to a file.
//
  m_ostring << m;
  n_ostring << n;
  center_ostring << center;

  output_filename = "grid_" + m_ostring.str ( ) + "_" 
    + n_ostring.str ( ) + "_" + center_ostring.str ( ) + ".txt";

  r8mat_write ( output_filename, m, n, r );

  cout << "\n";
  cout << "  The grid data was written to the file \""
    << output_filename << "\"\n";
//
//  Terminate.
//
  delete [] r;

  cout << "\n";
  cout << "GRID_DATASET\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

unsigned long get_seed ( )

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
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL

  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long ) 
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0E+00 * 60.0E+00 * 12.0E+00 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef UNSIGNED_LONG_MAX
}
//****************************************************************************80

void grid_generate ( int dim_num, int n, int center, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_GENERATE generates a grid dataset.
//
//  Discussion:
//
//    N points are needed in a DIM_NUM dimensional space.
//
//    The points are to lie on a uniform grid of side N_SIDE.
//
//    Unless the N = N_SIDE**DIM_NUM for some N_SIDE, we can't use all the
//    points on a grid.  What we do is find the smallest N_SIDE
//    that's big enough, and randomly omit some points.
//
//    If N_SIDE is 4, then the choices in 1D are:
//
//    A: 0,   1/3, 2/3, 1
//    B: 1/5, 2/5, 3/5, 4/5
//    C: 0,   1/4, 2/4, 3/4
//    D: 1/4, 2/4, 3/4, 1
//    E: 1/8, 3/8, 5/8, 7/8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points to generate.
//
//    Input, int CENTER, specifies the 1D grid centering:
//    1: first point is 0.0, last point is 1.0;
//    2: first point is 1/(N+1), last point is N/(N+1);
//    3: first point is 0, last point is (N-1)/N;
//    4: first point is 1/N, last point is 1;
//    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
//
//    Input/output, int *SEED, the random number seed.
//
//    Output, double R[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  int n_grid;
  int n_side;
  int rank;
  int *rank_list;
  int *tuple;

  rank_list = new int[n];
  tuple = new int[dim_num];
//
//  Find the dimension of the smallest grid with N points.
//
  n_side = grid_side ( dim_num, n );
//
//  We need to select N points out of N_SIDE**DIM_NUM set.
//
  n_grid = ( int ) pow ( ( double ) n_side, ( double ) dim_num );
//
//  Generate a random subset of N items from a set of size N_GRID.
//
  ksub_random2 ( n_grid, n, seed, rank_list );
//
//  Must make one dummy call to TUPLE_NEXT_FAST with RANK = 0.
//
  rank = 0;
  tuple_next_fast ( n_side, dim_num, rank, tuple );
//
//  Now generate the appropriate indices, and "center" them.
//
  for ( j = 0; j < n; j++ )
  {
    rank = rank_list[j] - 1;

    tuple_next_fast ( n_side, dim_num, rank, tuple );

    if ( center == 1 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] - 1 ) / ( double ) ( n_side - 1 );
      }
    }
    else if ( center == 2 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] ) / ( double ) ( n_side + 1 );
      }
    }
    else if ( center == 3 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] - 1 ) / ( double ) ( n_side );
      }
    }
    else if ( center == 4 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] ) / ( double ) ( n_side );
      }
    }
    else if ( center == 5 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( 2 * tuple[i] - 1 ) 
                       / ( double ) ( 2 * n_side );
      }
    }
  }

  delete [] rank_list;
  delete [] tuple;

  return;
}
//****************************************************************************80

int grid_side ( int dim_num, int n )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_SIDE finds the smallest DIM_NUM-dimensional grid containing at least N points.
//
//  Discussion:
//
//    Each coordinate of the grid will have N_SIDE distinct values.
//    Thus the total number of points in the grid is N_SIDE**DIM_NUM.
//    This routine seeks the smallest N_SIDE such that N <= N_SIDE**DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points to generate.
//
//    Output, int GRID_SIDE, the length of one side of the smallest 
//    grid in DIM_NUM dimensions that contains at least N points.
//
{
  double exponent;
  int n_grid;
  int n_side;

  if ( n <= 0 )
  {
    n_side = 0;
    return n_side;
  }

  if ( dim_num <= 0 )
  {
    n_side = -1;
    return n_side;
  }

  exponent = 1.0E+00 / ( double ) dim_num;

  n_side = ( int ) pow ( n, exponent );

  if ( pow ( ( double ) n_side, ( double ) dim_num ) < n )
  {
    n_side = n_side + 1;
  }

  return n_side;
}
//****************************************************************************80

void ksub_random2 ( int n, int k, int *seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    FORTRAN90 version by John Burkardt.
//
//  Reference:
//
//    A Nijenhuis and H Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, int A(K).  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int available;
  int candidate;
  int have;
  int need;
  double r;

  if ( k < 0 || n < k )
  {
    cout << "\n";
    cout << "KSUB_RANDOM2 - Fatal error!\n";
    cout << "  N = " << n << "\n";
    cout << "  K = " << k << "\n";
    cout << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( k == 0 )
  {
    return;
  }

  need = k;
  have = 0;
  available = n;
  candidate = 0;

  for ( ; ; )
  {
    candidate = candidate + 1;

    r = r8_uniform_01 ( seed );

    if ( r * ( double ) available <= ( double ) need )
    {
      need = need - 1;
      a[have] = candidate;
      have = have + 1;

      if ( need <= 0 )
      {
        break;
      }

    }

    available = available - 1;

  }

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit double precision pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
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
//    11 August 2004
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
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
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
//****************************************************************************80

void tuple_next_fast ( int m, int n, int rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      X
//    ----  ---    ---
//    0     * *    1 1
//    1     1 1    1 2
//    2     1 2    1 3
//    3     1 3    2 1
//    4     2 1    2 2
//    5     2 2    2 3
//    6     2 3    3 1
//    7     3 1    3 2
//    8     3 2    3 3
//    9     3 3    1 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the maximum entry.
//
//    Input, int N, the number of components.
//
//    Input, int RANK, indicates the rank of the tuples.
//    On the very first call only, it is necessary that
//    the user set RANK = 0.  
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  static int *base = NULL;
  int i;

  if ( rank == 0 )
  {
    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; i >= 0; i-- )
    {
      base[i] = base[i+1] * m;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( rank / base[i] ) % m ) + 1;
  }

  return;
}
