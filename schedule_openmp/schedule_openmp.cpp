# include <cstdlib>
# include <iostream>
# include <iomanip>

# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
int prime_default ( int n );
int prime_static ( int n );
int prime_dynamic ( int n );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SCHEDULE_OPENMP.
//
//  Discussion:
//
//    This program demonstrates the difference between default,
//    static and dynamic scheduling for a loop parallelized in OpenMP.
//
//    The purpose of scheduling is to deal with loops in which there is
//    known or suspected imbalance in the work load.  In this example,
//    if the work is divided in the default manner between two threads,
//    the second thread has 3 times the work of the first.  
//
//    Both static and dynamic scheduling, if used, even out the work
//    so that both threads have about the same load.  This could be
//    expected to decrease the run time of the loop by about 1/3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int n_factor;
  int n_hi;
  int n_lo;
  int primes;
  double time1;
  double time2;
  double time3;

  cout << "\n";
  cout << "SCHEDULE_OPENMP\n";
  cout << "  C++/OpenMP version\n";
  cout << "  Count the primes from 1 to N.\n";
  cout << "  This is an unbalanced work load, particular for two threads.\n";
  cout << "  Demonstrate default, static and dynamic scheduling.\n";
  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  n_lo = 1;
  n_hi = 131072;
  n_factor = 2;

  cout << "\n";
  cout << "                           Default        Static       Dynamic\n";
  cout << "         N     Pi(N)          Time          Time          Time\n";
  cout << "\n";

  n = n_lo;

  while ( n <= n_hi )
  {
    time1 = omp_get_wtime ( );
    primes = prime_default ( n );
    time1 = omp_get_wtime ( ) - time1;

    time2 = omp_get_wtime ( );
    primes = prime_static ( n );
    time2 = omp_get_wtime ( ) - time2;

    time3 = omp_get_wtime ( );
    primes = prime_dynamic ( n );
    time3 = omp_get_wtime ( ) - time3;

    cout << "  " << setw(8) << n
         << "  " << setw(8) << primes
         << "  " << setw(12) << time1
         << "  " << setw(12) << time2
         << "  " << setw(12) << time3 << "\n";

    n = n * n_factor;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "SCHEDULE_OPENMP\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

int prime_default ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_DEFAULT counts primes, using default scheduling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number to check.
//
//    Output, int PRIME_DEFAULT, the number of prime numbers up to N.
//
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
//****************************************************************************80

int prime_static ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_STATIC counts primes using static scheduling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number to check.
//
//    Output, int PRIME_STATIC, the number of prime numbers up to N.
//
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total ) schedule ( static, 100 )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
//****************************************************************************80

int prime_dynamic ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_DYNAMIC counts primes using dynamic scheduling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number to check.
//
//    Output, int PRIME_DYNAMIC, the number of prime numbers up to N.
//
{
  int i;
  int j;
  int prime;
  int total = 0;

# pragma omp parallel \
  shared ( n ) \
  private ( i, j, prime )

# pragma omp for reduction ( + : total ) schedule ( dynamic, 100 )
  for ( i = 2; i <= n; i++ )
  {
    prime = 1;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = 0;
        break;
      }
    }
    total = total + prime;
  }

  return total;
}
