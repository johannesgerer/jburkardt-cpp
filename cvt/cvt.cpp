# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>

using namespace std;

# include "cvt.hpp"

//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 ) 
  {
    c = c - 32;
  }   

  return c;
}
//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 ) 
  {
    c1 = c1 - 32;
  } 
  if ( 97 <= c2 && c2 <= 122 ) 
  {
    c2 = c2 - 32;
  }     

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

void cvt ( int dim_num, int n, int batch, int init, int sample, int sample_num, 
  int it_max, int it_fixed, int *seed, double r[],int *it_num, double *it_diff, 
  double *energy )

//****************************************************************************80
//
//  Purpose:
//
//    CVT computes a Centroidal Voronoi Tessellation.
//
//  Discussion:
//
//    This routine initializes the data, and carries out the
//    CVT iteration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int INIT, specifies how the points are to be initialized.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine;
//     4, points are already initialized on input.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Input, int IT_FIXED, the maximum number of iterations to take
//    with a fixed set of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the approximate CVT points.
//    If INIT = 4 on input, then it is assumed that these values have been
//    initialized.  On output, the CVT iteration has been applied to improve
//    the value of the points.
//
//    Output, int *IT_NUM, the number of iterations taken.  Generally,
//    this will be equal to IT_MAX, unless the iteration tolerance was
//    satisfied early.
//
//    Output, double *IT_DIFF, the L2 norm of the difference
//    between the iterates.
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  bool DEBUG = true;
  int i;
  bool initialize;
  int seed_base;
  int seed_init;

  if ( batch < 1 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value BATCH < 1.\n";
    exit ( 1 );
  }

  if ( seed <= 0 )
  {
    cout << "\n";
    cout << "CVT - Fatal error!\n";
    cout << "  The input value SEED <= 0.\n";
    exit ( 1 );
  }

  if ( DEBUG )
  {
    cout << "\n";
    cout << "  Step       SEED          L2-Change        Energy\n";
    cout << "\n";
  }

  *it_num = 0;
  *it_diff = 0.0;
  *energy = 0.0;
  seed_init = *seed;
//
//  Initialize the data, unless the user has already done that.
//
  if ( init != 4 )
  {
    initialize = true;
    cvt_sample ( dim_num, n, n, init, initialize, seed, r );
  }
  if ( DEBUG )
  {
    cout                          << "  "
         << setw(4)  << *it_num   << "  "
         << setw(12) << seed_init << "\n";
  }
//
//  If the initialization and sampling steps use the same random number
//  scheme, then the sampling scheme does not have to be initialized.
//
  if ( init == sample )
  {
    initialize = false;
  }
  else
  {
    initialize = true;
  }
//
//  Carry out the iteration.
//
  while ( *it_num < it_max )
  {
//
//  If it's time to update the seed, save its current value
//  as the starting value for all iterations in this cycle.
//  If it's not time to update the seed, restore it to its initial
//  value for this cycle.
//
    if ( ( (*it_num) % it_fixed ) == 0 )
    {
      seed_base = *seed;
    }
    else
    {
      *seed = seed_base;
    }

    *it_num = *it_num + 1;
    seed_init = *seed;

    cvt_iterate ( dim_num, n, batch, sample, initialize, sample_num, seed,
      r, it_diff, energy );

    initialize = false;

    if ( DEBUG )
    {
      cout                          << "  "
           << setw(4)  << *it_num   << "  "
           << setw(12) << seed_init << "  "
           << setw(14) << *it_diff  << "  "
           << setw(14) << *energy    << "\n";
    }
  }
  return;
}
//****************************************************************************80

double cvt_energy ( int dim_num, int n, int batch, int sample, bool initialize,
  int sample_num, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_ENERGY computes the CVT energy of a dataset.
//
//  Discussion:
//
//    For a given number of generators, a CVT is a minimizer (or at least
//    a local minimizer) of the CVT energy.  During a CVT iteration,
//    it should generally be the case that the CVT energy decreases from
//    step to step, and that perturbations or adjustments of an
//    approximate CVT will almost always have higher CVT energy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of generators.
//
//    Input, int BATCH, the maximum number of sample points to generate
//    at one time.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, bool INITIALIZE, is TRUE if the pseudorandom process 
//    should be reinitialized.
//
//    Input, int SAMPLE_NUM, the number of sample points to use.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Input, double R[DIM_NUM*N], the coordinates of the points.
//
//    Output, double CVT_ENERGY, the estimated CVT energy.
//
{
  double energy;
  int get;
  int have;
  int i;
  int j;
  int *nearest;
  double *s;

  nearest = new int[batch];
  s = new double [dim_num*batch];

  have = 0;
  energy = 0.0;

  while ( have < sample_num )
  {
    get = i4_min ( sample_num - have, batch );

    cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s );

    have = have + get;

    find_closest ( dim_num, n, get, s, r, nearest );

    for ( j = 0; j < get; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        energy = energy + ( s[i+j*dim_num] - r[i+nearest[j]*dim_num] )
                        * ( s[i+j*dim_num] - r[i+nearest[j]*dim_num] );
      }
    }

  }

  energy = energy / ( double ) ( sample_num );

  delete [] nearest;
  delete [] s;

  return energy;
}
//****************************************************************************80

void cvt_iterate ( int dim_num, int n, int batch, int sample, bool initialize, 
  int sample_num, int *seed, double r[], double *it_diff, double *energy )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_ITERATE takes one step of the CVT iteration.
//
//  Discussion:
//
//    The routine is given a set of points, called "generators", which
//    define a tessellation of the region into Voronoi cells.  Each point
//    defines a cell.  Each cell, in turn, has a centroid, but it is
//    unlikely that the centroid and the generator coincide.
//
//    Each time this CVT iteration is carried out, an attempt is made
//    to modify the generators in such a way that they are closer and
//    closer to being the centroids of the Voronoi cells they generate.
//
//    A large number of sample points are generated, and the nearest generator
//    is determined.  A count is kept of how many points were nearest to each
//    generator.  Once the sampling is completed, the location of all the
//    generators is adjusted.  This step should decrease the discrepancy
//    between the generators and the centroids.
//
//    The centroidal Voronoi tessellation minimizes the "energy",
//    defined to be the integral, over the region, of the square of
//    the distance between each point in the region and its nearest generator.
//    The sampling technique supplies a discrete estimate of this
//    energy.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int BATCH, sets the maximum number of sample points
//    generated at one time.  It is inefficient to generate the sample
//    points 1 at a time, but memory intensive to generate them all
//    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
//    BATCH must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, bool INITIALIZE, is TRUE if the SEED must be reset to SEED_INIT
//    before computation.  Also, the pseudorandom process may need to be
//    reinitialized.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input/output, int *SEED, the random number seed.
//
//    Input/output, double R[DIM_NUM*N], the Voronoi
//    cell generators.  On output, these have been modified
//
//    Output, double *IT_DIFF, the L2 norm of the difference
//    between the iterates.
//
//    Output, double *ENERGY,  the discrete "energy", divided
//    by the number of sample points.
//
{
  int *count;
  int get;
  int have;
  int i;
  int j;
  int j2;
  int *nearest;
  double *r2;
  double *s;
  bool success;
  double term;
//
//  Take each generator as the first sample point for its region.
//  This can slightly slow the convergence, but it simplifies the
//  algorithm by guaranteeing that no region is completely missed
//  by the sampling.
//
  *energy = 0.0;
  r2 = new double[dim_num*n];
  count = new int[n];
  nearest = new int[sample_num];
  s = new double[dim_num*sample_num];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r[i+j*dim_num];
    }
  }
  for ( j = 0; j < n; j++ )
  {
    count[j] = 1;
  }
//
//  Generate the sampling points S.
//
  have = 0;

  while ( have < sample_num )
  {
    get = i4_min ( sample_num - have, batch );

    cvt_sample ( dim_num, sample_num, get, sample, initialize, seed, s );

    initialize = false;
    have = have + get;
//
//  Find the index N of the nearest cell generator to each sample point S.
//
    find_closest ( dim_num, n, get, s, r, nearest );
//
//  Add S to the centroid associated with generator N.
//
    for ( j = 0; j < get; j++ )
    {
      j2 = nearest[j];
      for ( i = 0; i < dim_num; i++ )
      {
        r2[i+j2*dim_num] = r2[i+j2*dim_num] + s[i+j*dim_num];
      }
      for ( i = 0; i < dim_num; i++ )
      {
        *energy = *energy + pow ( r[i+j2*dim_num] - s[i+j*dim_num], 2 );
      }
      count[j2] = count[j2] + 1;
    }
  }
//
//  Estimate the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r2[i+j*dim_num] = r2[i+j*dim_num] / ( double ) ( count[j] );
    }
  }
//
//  Determine the sum of the distances between generators and centroids.
//
  *it_diff = 0.0;

  for ( j = 0; j < n; j++ )
  {
    term = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      term = term + ( r2[i+j*dim_num] - r[i+j*dim_num] ) 
                  * ( r2[i+j*dim_num] - r[i+j*dim_num] );
    }
    *it_diff = *it_diff + sqrt ( term );
  }
//
//  Replace the generators by the centroids.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      r[i+j*dim_num] = r2[i+j*dim_num];
    }
  }
//
//  Normalize the discrete energy estimate.
//
  *energy = *energy / sample_num;

  delete [] count;
  delete [] nearest;
  delete [] r2;
  delete [] s;

  return;
}
//****************************************************************************80

void cvt_sample ( int dim_num, int n, int n_now, int sample, bool initialize, 
  int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_SAMPLE returns sample points.
//
//  Discussion:
//
//    N sample points are to be taken from the unit box of dimension DIM_NUM.
//
//    These sample points are usually created by a pseudorandom process
//    for which the points are essentially indexed by a quantity called
//    SEED.  To get N sample points, we generate values with indices
//    SEED through SEED+N-1.
//
//    It may not be practical to generate all the sample points in a 
//    single call.  For that reason, the routine allows the user to
//    request a total of N points, but to require that only N_NOW be
//    generated now (on this call).  
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input, int N_NOW, the number of sample points to be generated
//    on this call.  N_NOW must be at least 1.
//
//    Input, int SAMPLE, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input, bool INITIALIZE, is TRUE if the pseudorandom process should be
//    reinitialized.
//
//    Input/output, int *SEED, the random number seed.
//
//    Output, double R[DIM_NUM*N_NOW], the sample points.
//
{
  double exponent;
  static int *halton_base = NULL;
  static int *halton_leap = NULL;
  static int *halton_seed = NULL;
  int halton_step;
  int i;
  int j;
  int k;
  static int ngrid;
  static int rank;
  int rank_max;
  static int *tuple = NULL;

  if ( n_now < 1 )
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  N_NOW < 1.\n";
    exit ( 1 );
  }

  if ( sample == -1 )
  {
    if ( initialize )
    {
      random_initialize ( *seed );
    }

    for ( j = 0; j < n_now; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) random ( ) / ( double ) RAND_MAX;
      }
    }
    *seed = ( *seed ) + n_now * dim_num;
  }
  else if ( sample == 0 )
  {
    r8mat_uniform_01 ( dim_num, n_now, seed, r );
  }
  else if ( sample == 1 )
  {
    halton_seed = new int[dim_num];
    halton_leap = new int[dim_num];
    halton_base = new int[dim_num];
      
    halton_step = *seed;

    for ( i = 0; i < dim_num; i++ )
    {
      halton_seed[i] = 0;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_leap[i] = 1;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      halton_base[i] = prime ( i + 1 );
    }

    i4_to_halton_sequence ( dim_num, n_now, halton_step, halton_seed,
      halton_leap, halton_base, r );

    delete [] halton_seed;
    delete [] halton_leap;
    delete [] halton_base;

    *seed = *seed + n_now;
  }
  else if ( sample == 2 )
  {
    exponent = 1.0 / ( double ) ( dim_num );
    ngrid = ( int ) pow ( ( double ) n, exponent );
    rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );
    tuple = new int[dim_num];

    if ( rank_max < n )
    {
      ngrid = ngrid + 1;
      rank_max = ( int ) pow ( ( double ) ngrid, ( double ) dim_num );
    }

    if ( initialize )
    {
      rank = -1;
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
    }

    rank = ( *seed ) % rank_max;

    for ( j = 0; j < n_now; j++ )
    {
      tuple_next_fast ( ngrid, dim_num, rank, tuple );
      rank = rank + 1;
      rank = rank % rank_max;
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = double ( 2 * tuple[i] - 1 ) / double ( 2 * ngrid );
      }
    }
    delete [] tuple;
    *seed = *seed + n_now;
  }
  else if ( sample == 3 )
  {
    user ( dim_num, n_now, seed, r );
  }
  else
  {
    cout << "\n";
    cout << "CVT_SAMPLE - Fatal error!\n";
    cout << "  The value of SAMPLE = " << sample << " is illegal.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void data_read ( char *file_in_name, int dim_num, int n, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    DATA_READ reads generator coordinate data from a file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_IN_NAME, the name of the input file.
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R[DIM_NUM*N], the point coordinates.
//
{
  bool error;
  ifstream file_in;
  int i;
  int j;
  char line[255];
  double *x;

  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << file_in_name << "\"\n";
    exit ( 1 );
  }

  x = new double[dim_num];

  j = 0;

  while ( j < n )
  {
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, dim_num, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      r[i+j*dim_num] = x[i];
    }
    j = j + 1;

  }

  file_in.close ( );

  delete [] x;

  cout << "\n";
  cout << "DATA_READ:\n";
  cout << "  Read coordinate data from file.\n";

  return;
}
//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'  
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

void find_closest ( int dim_num, int n, int sample_num, double s[], double r[],
  int nearest[] )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_CLOSEST finds the nearest R point to each S point.
//
//  Discussion:
//
//    This routine finds the closest Voronoi cell generator by checking every
//    one.  For problems with many cells, this process can take the bulk
//    of the CPU time.  Other approaches, which group the cell generators into
//    bins, can run faster by a large factor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of cell generators.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, double S[DIM_NUM*SAMPLE_NUM], the points to be checked.
//
//    Input, double R[DIM_NUM*N], the cell generators.
//
//    Output, int NEAREST[SAMPLE_NUM], the (0-based) index of the nearest 
//    cell generator.
//
{
  double dist_sq_min;
  double dist_sq;
  int i;
  int jr;
  int js;

  for ( js = 0; js < sample_num; js++ )
  {
    dist_sq_min = r8_huge ( );
    nearest[js] = -1;

    for ( jr = 0; jr < n; jr++ )
    {
      dist_sq = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        dist_sq = dist_sq + ( s[i+js*dim_num] - r[i+jr*dim_num] ) 
                          * ( s[i+js*dim_num] - r[i+jr*dim_num] );
      }

      if ( jr == 0 || dist_sq < dist_sq_min )
      {
        dist_sq_min = dist_sq;
        nearest[js] = jr;
      }
    }
  }

  return;
}
//****************************************************************************80

int get_seed ( void )

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
//  Remap ISEED from [1,43200] to [1,IMAX].
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

bool halham_leap_check ( int dim_num, int leap[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEAP[DIM_NUM], the successive jumps in the sequence.
//    Each entry must be greater than 0.
//
//    Output, bool HALHAM_LEAP_CHECK, is true if LEAP is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( leap[i] < 1 ) 
    {
      cout << "\n";
      cout << "HALHAM_LEAP_CHECK - Fatal error!\n";
      cout << "  Leap entries must be greater than 0.\n";
      cout << "  leap[" << i << "] = " << leap[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//****************************************************************************80

bool halham_n_check ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of the subsequence.
//    N must be positive.
//
//    Output, bool HALHAM_N_CHECK, is true if N is legal.
//
{
  bool value;

  if ( n < 1 ) 
  {
    cout << "\n";
    cout << "HALHAM_N_CHECK - Fatal error!\n";
    cout << "  N < 0.";
    cout << "  N = " << n << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//****************************************************************************80

bool halham_dim_num_check ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//    DIM_NUM must be positive.
//
//    Output, bool HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
//
{
  bool value;

  if ( dim_num < 1 ) 
  {
    cout << "\n";
    cout << "HALHAM_DIM_NUM_CHECK - Fatal error!\n";
    cout << "  DIM_NUM < 0.";
    cout << "  DIM_NUM = " << dim_num << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//****************************************************************************80

bool halham_seed_check ( int dim_num, int seed[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int SEED[DIM_NUM], the sequence index
//    corresponding to STEP = 0.  Each entry must be 0 or greater.
//
//    Output, bool HALHAM_SEED_CHECK, is true if SEED is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( seed[i] < 0 ) 
    {
      cout << "\n";
      cout << "HALHAM_SEED_CHECK - Fatal error!\n";
      cout << "  SEED entries must be nonnegative.\n";
      cout << "  seed[" << i << "] = " << seed[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//****************************************************************************80

bool halham_step_check ( int step )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the index of the subsequence element.
//    STEP must be 1 or greater.
//
//    Output, bool HALHAM_STEP_CHECK, is true if STEP is legal.
//
{
  int i;
  bool value;

  if ( step < 0 ) 
  {
    cout << "\n";
    cout << "HALHAM_STEP_CHECK - Fatal error!\n";
    cout << "  STEP < 0.";
    cout << "  STEP = " << step << "\n";
    value = false;
  }
  else
  {
    value = true;
  }

  return value;
}
//****************************************************************************80

bool halton_base_check ( int dim_num, int base[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_BASE_CHECK checks BASE for a Halton sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases.
//
//    Output, bool HALTON_BASE_CHECK, is true if BASE is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( base[i] <= 1 ) 
    {
      cout << "\n";
      cout << "HALTON_BASE_CHECK - Fatal error!\n";
      cout << "  Bases must be greater than 1.\n";
      cout << "  base[" << i << "] = " << base[i] << "\n";
      value = false;
      break;
    }
  }

  return value;
}
//****************************************************************************80

int i4_log_10 ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_LOG_10 returns the whole part of the logarithm base 10 of an I4.
//
//  Discussion:
//
//    It should be the case that 10^I4_LOG_10(I) <= |I| < 10^(I4_LOG_10(I)+1).
//    (except for I = 0).
//
//    The number of decimal digits in I is I4_LOG_10(I) + 1.
//
//  Example:
//
//        I    I4_LOG_10(I)
//
//        0     0
//        1     0
//        2     0
//
//        9     0
//       10     1
//       11     1
//
//       99     1
//      100     2
//      101     2
//
//      999     2
//     1000     3
//     1001     3
//
//     9999     3
//    10000     4
//    10001     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer.
//
//    Output, int I4_LOG_10, the whole part of the logarithm of abs ( I ).
//
{
  int ten_pow;
  int value;

  i = abs ( i );

  ten_pow = 10;
  value = 0;

  while ( ten_pow <= i )
  {
    ten_pow = ten_pow * 10;
    value = value + 1;
  }

  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1 and I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 ) 
  {
    return i1;
  }
  else 
  {
    return i2;
  }

}
//**********************************************************************

void i4_to_halton_sequence ( int dim_num, int n, int step, int seed[], int leap[],
  int base[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Halton sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
//
//    of the original Halton sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is DIM_NUM * N, where DIM_NUM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first DIM_NUM entries of R, followed by the DIM_NUM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*DIM_NUM through (DIM_NUM-1) + (J-1)*DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements of the sequence.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
//    to STEP = 0. 
//
//    Input, int LEAP[DIM_NUM], the succesive jumps in the Halton sequence.
//
//    Input, int BASE[DIM_NUM], the Halton bases.
//
//    Output, double R[DIM_NUM*N], the next N elements of the
//    leaped Halton subsequence, beginning with element STEP.
//
{
  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }
  
  if ( !halham_n_check ( n ) )
  {
    exit ( 1 );
  }
  
  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }
  
  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < dim_num; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      seed2[j] = seed[i] + ( step + j ) * leap[i];
    }

    for ( j = 0; j < n; j++ )
    {
      r[i+j*dim_num] = 0.0;
    }

    for ( j = 0; j < n; j++ )
    {
      base_inv = 1.0 / ( ( double ) base[i] );

      while ( seed2[j] != 0 )
      {
        digit = seed2[j] % base[i];
        r[i+j*dim_num] = r[i+j*dim_num] + ( ( double ) digit ) * base_inv;
        base_inv = base_inv / ( ( double ) base[i] );
        seed2[j] = seed2[j] / base[i];
      }
    }
  }

  delete [] seed2;

  return;
}
//****************************************************************************80

char *i4_to_s ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_S converts an I4 to a string.
//
//  Example:
//
//    INTVAL  S
//
//         1  1
//        -1  -1
//         0  0
//      1952  1952
//    123456  123456
//   1234567  1234567
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer to be converted.
//
//    Output, char *I4_TO_S, the representation of the integer.
//
{
  int digit;
  int j;
  int length;
  int ten_power;
  char *s;

  length = i4_log_10 ( i );

  ten_power = ( int ) pow ( ( double ) 10, ( double ) length );

  if ( i < 0 )
  {
    length = length + 1;
  }
//
//  Add one position for the trailing null.
//
  length = length + 1;

  s = new char[length];

  if ( i == 0 )
  {
    s[0] = '0';
    s[1] = '\0';
    return s;
  }
//
//  Now take care of the sign.
//
  j = 0;
  if ( i < 0 )
  {
    s[j] = '-';
    j = j + 1;
    i = abs ( i );
  }
//
//  Find the leading digit of I, strip it off, and stick it into the string.
//
  while ( 0 < ten_power )
  {
    digit = i / ten_power;
    s[j] = digit_to_ch ( digit );
    j = j + 1;
    i = i - digit * ten_power;
    ten_power = ten_power / 10;
  }
//
//  Tack on the trailing NULL.
//
  s[j] = '\0';
  j = j + 1;

  return s;
}
//****************************************************************************80

int prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, pages 870-873.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cout << "\n";
    cout << "PRIME - Fatal error!\n";
    cout << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal real number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
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
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, char *TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
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
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, char *TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
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
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  On output, SEED has 
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return;
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

unsigned long random_initialize ( int seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator, 
//    it will behave as though it were seeded with value 1.  
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRAND (the appropriate routine 
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRAND, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  unsigned long ul_seed;

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed ( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RAND with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  ul_seed = ( unsigned long ) seed;
  srand ( ul_seed );

  return ul_seed;
# undef DEBUG
}
//****************************************************************************80

bool s_eqi ( char *s1, char *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S1, char *S2, pointers to two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  nchar = i4_min ( nchar1, nchar2 );

//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < nchar1 ) 
  {
    for ( i = nchar; i < nchar1; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < nchar2 ) 
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;

}
//****************************************************************************80

int s_len_trim ( char* s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n ) 
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

double s_to_r8 ( char *s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( char *s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;
  double x;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

void timestamp ( void )

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
//    02 October 2003
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
//    Rank          X
//    ----          ----
//   -1            -1 -1
//
//    0             1  1
//    1             1  2
//    2             1  3
//    3             2  1
//    4             2  2
//    5             2  3
//    6             3  1
//    7             3  2
//    8             3  3
//    9             1  1
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
//  Parameters:
//
//    Input, int M, the maximum entry in each component.
//    M must be greater than 0.
//
//    Input, int N, the number of components.
//    N must be greater than 0.
//
//    Input, integer RANK, indicates the rank of the tuples.
//    Typically, 0 <= RANK < N**M; values larger than this are legal
//    and meaningful, and are equivalent to the corresponding value
//    MOD N**M.  If RANK < 0, this indicates that this is the first
//    call for the given values of (M,N).  Initialization is done,
//    and X is set to a dummy value.
//
//    Output, int X[N], the next tuple, or a dummy value if initialization
//    is being done.
//
{
  static int *base = NULL;
  int i;
//
  if ( rank < 0 )
  {
    if ( m <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value M <= 0 is not legal.\n";
      cout << "  M = " << m << "\n";
      exit ( 1 );
    }
    if ( n <= 0 )
    {
      cout << "\n";
      cout << "TUPLE_NEXT_FAST - Fatal error!\n";
      cout << "  The value N <= 0 is not legal.\n";
      cout << "  N = " << n << "\n";
      exit ( 1 );
    }

    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
    for ( i = 0; i < n; i++ )
    {
      x[i] = -1;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( rank / base[i] ) % m ) + 1;
    }
  }
  return;
}
//****************************************************************************80

void user ( int dim_num, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    USER samples points in a user-specified region with given density.
//
//  Discussion:
//
//    This routine can be used to 
//
//    * specify an interesting initial configuration for the data,
//      by specifing that USER be used for initialization (INIT = 3);
//
//    * specify the shape of the computational region, by specifying
//      that sample points are to be generated by this routine, 
//      (SAMPLE = 3) and then returning sample points uniformly at random.
//
//    * specify the distribution or density function, by specifying
//      that sample points are to be generated by this routine, 
//      (SAMPLE = 3 ) and then returning sample points according to a 
//      given probability density function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//
//    Input, integer N, the number of sample points desired.
//
//    Input/output, int *SEED, the "seed" value.  On output, SEED has 
//    been updated.
//
//    Output, double R[DIM_NUM*N], the sample values.
//
{
# define PI 3.141592653589793

  double angle;
  int j;
  double radius;

  for ( j = 0; j < n; j++ )
  {
    angle = 2.0 * PI * ( double ) random ( ) / ( double ) RAND_MAX;
    radius = sqrt ( ( double ) random ( ) / ( double ) RAND_MAX );
    r[0+j*2] = radius * cos ( angle );
    r[1+j*2] = radius * sin ( angle );
  }

  return;
# undef PI
}
