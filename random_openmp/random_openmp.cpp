# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include <omp.h>

using namespace std;

int main ( );
void monte_carlo ( int n, int &seed );
double random_value ( int &seed );
void timestamp ( );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for RANDOM_OPENMP.

  Discussion:

    This program simply explores one issue in the generation of random
    numbers in a parallel program.  If the random number generator uses
    an integer seed to determine the next entry, then it is not easy for
    a parallel program to reproduce the same exact sequence.

    But what is worse is that it might not be clear how the separate
    OpenMP threads should handle the SEED value - as a shared or private
    variable?  It seems clear that each thread should have a private
    seed that is initialized to a distinct value at the beginning of
    the computation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt
*/
{
  int n;
  int seed;

  timestamp ( );

  cout << "\n";
  cout << "RANDOM_OPENMP\n";
  cout << "  C++ version\n";
  cout << "  An OpenMP program using random numbers.\n";
  cout << "  The random numbers depend on a seed.\n";
  cout << "  We need to insure that each OpenMP thread\n";
  cout << "  starts with a different seed.\n";
  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  n = 100;
  seed = 123456789;
  monte_carlo ( n, seed );
/*
  Terminate.
*/
  cout << "\n";
  cout << "RANDOM_OPENMP\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
/******************************************************************************/

void monte_carlo ( int n, int &seed )

/******************************************************************************/
/*
  Purpose:

    MONTE_CARLO carries out a Monte Carlo calculation with random values.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt

  Parameter:

    Input, int N, the number of values to generate.

    Input, int &SEED, a seed for the random number generator.
*/
{
  int i;
  int my_id;
  int *my_id_vec;
  int my_seed;
  int *my_seed_vec;
  double *x;

  x = new double[n];

  my_id_vec = new int[n];
  my_seed_vec = new int[n];

# pragma omp master
{
  cout << "\n";
  cout << "  Thread   Seed  I   X(I)\n";
  cout << "\n";
}

# pragma omp parallel private ( i, my_id, my_seed ) shared ( my_id_vec, my_seed_vec, n, x )
{
  my_id = omp_get_thread_num ( );
  my_seed = seed + my_id;
  cout << "  " << setw(6) << my_id
       << "  " << setw(12) << my_seed << "\n";

# pragma omp for
  for ( i = 0; i < n; i++ )
  {
    my_id_vec[i] = my_id;
    x[i] = random_value ( my_seed );
    my_seed_vec[i] = my_seed;
//  cout << "  " << setw(6) << my_id
//       << "  " << setw(12) << my_seed
//       << "  " << setw(6) << i
//       << "  " << setw(14) << x[i] << "\n";
  }

}
//
//  C++ OpenMP IO from multiple processors comes out chaotically.
//  For this reason only, we'll save the data from the loop and
//  print it in the sequential section!
//
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(6) << my_id_vec[i]
         << "  " << setw(12) << my_seed_vec[i]
         << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }

  delete [] my_id_vec;
  delete [] my_seed_vec;
  delete [] x;

  return;
}
/******************************************************************************/

double random_value ( int &seed )

/******************************************************************************/
/*
  Purpose:

    RANDOM_VALUE generates a random value R.

  Discussion:

    This is not a good random number generator.  It is a SIMPLE one.
    It illustrates a model which works by accepting an integer seed value
    as input, performing some simple operation on the seed, and then
    producing a "random" real value using some simple transformation.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    03 September 2012

  Author:

    John Burkardt

  Parameters:

    Input/output, int &SEED, a seed for the random 
    number generator.

    Output, double RANDOM_VALUE, the random value.
*/
{
  double r;

  seed = ( seed % 65536 );
  seed = ( ( 3125 * seed ) % 65536 );
  r = ( double ) ( seed ) / 65536.0;

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

