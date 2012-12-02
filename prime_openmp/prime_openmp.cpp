# include <cstdlib>
# include <iostream>
# include <iomanip>

# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
void prime_number_sweep ( int n_lo, int n_hi, int n_factor );
int prime_number ( int n );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRIME_OPENMP.
//
//  Discussion:
//
//    This program calls a version of PRIME_NUMBER that includes
//    OpenMP directives for parallel processing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int n_factor;
  int n_hi;
  int n_lo;

  cout << "\n";
  cout << "PRIME_OPENMP\n";
  cout << "  C++/OpenMP version\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  n_lo = 1;
  n_hi = 131072;
  n_factor = 2;

  prime_number_sweep ( n_lo, n_hi, n_factor );

  n_lo = 5;
  n_hi = 500000;
  n_factor = 10;

  prime_number_sweep ( n_lo, n_hi, n_factor );
//
//  Terminate.
//
  cout << "\n";
  cout << "PRIME_OPENMP\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void prime_number_sweep ( int n_lo, int n_hi, int n_factor )

//****************************************************************************80
//
//  Purpose:
//
//   PRIME_NUMBER_SWEEP does repeated calls to PRIME_NUMBER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LO, the first value of N.
//
//    Input, int N_HI, the last value of N.
//
//    Input, int N_FACTOR, the factor by which to increase N after
//    each iteration.
//
{
  int i;
  int n;
  int primes;
  double wtime;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Call PRIME_NUMBER to count the primes from 1 to N.\n";
  cout << "\n";
  cout << "         N        Pi          Time\n";
  cout << "\n";

  n = n_lo;

  while ( n <= n_hi )
  {
    wtime = omp_get_wtime ( );

    primes = prime_number ( n );

    wtime = omp_get_wtime ( ) - wtime;

    cout << "  " << setw(8) << n
         << "  " << setw(8) << primes
         << "  " << setw(14) << wtime << "\n";

    n = n * n_factor;
  }
 
  return;
}
//****************************************************************************80

int prime_number ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_NUMBER returns the number of primes between 1 and N.
//
//  Discussion:
//
//    A naive algorithm is used.
//
//    Mathematica can return the number of primes less than or equal to N
//    by the command PrimePi[N].
//
//                N  PRIME_NUMBER
//
//                1           0
//               10           4
//              100          25
//            1,000         168
//           10,000       1,229
//          100,000       9,592
//        1,000,000      78,498
//       10,000,000     664,579
//      100,000,000   5,761,455
//    1,000,000,000  50,847,534
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the maximum number to check.
//
//    Output, int PRIME_NUMBER, the number of prime numbers up to N.
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
