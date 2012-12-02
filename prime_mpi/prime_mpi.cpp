# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include "mpi.h"

using namespace std;

int main ( int argc, char *argv[] );
int prime_number ( int n, int id, int p );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRIME_MPI.
//
//  Discussion:
//
//    This program calls a version of PRIME_NUMBER that includes
//    MPI calls for parallel processing.
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
{
  int i;
  int id;
  int n;
  int n_factor;
  int n_hi;
  int n_lo;
  int p;
  int primes;
  int primes_part;
  double wtime;

  n_lo = 1;
  n_hi = 262144;
  n_factor = 2;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get the number of processes.
//
  p = MPI::COMM_WORLD.Get_size (  );
//
//  Determine this processes's rank.
//
  id = MPI::COMM_WORLD.Get_rank ( );

  if ( id == 0 )
  {
    timestamp ( );
    cout << "\n";
    cout << "PRIME_MPI\n";
    cout << "  C++/MPI version\n";
    cout << "\n";
    cout << "  An MPI example program to count the number of primes.\n";
    cout << "  The number of processes is " << p << "\n";
    cout << "\n";
    cout << "         N        Pi          Time\n";
    cout << "\n";
  }

  n = n_lo;

  while ( n <= n_hi )
  {
    if ( id == 0 )
    {
      wtime = MPI::Wtime ( );
    }
    MPI::COMM_WORLD.Bcast ( &n, 1, MPI::INT, 0 );

    primes_part = prime_number ( n, id, p );

    MPI::COMM_WORLD.Reduce ( &primes_part, &primes, 1, MPI::INT, MPI::SUM, 
      0 );

    if ( id == 0 )
    {
      wtime = MPI::Wtime ( ) - wtime;

      cout << "  " << setw(8) << n
           << "  " << setw(8) << primes
           << "  " << setw(14) << wtime << "\n";
    }
    n = n * n_factor;
  }
//
//  Terminate MPI.
//
  MPI::Finalize ( );
//
//  Terminate.
//
  if ( id == 0 ) 
  {
    cout << "\n";
    cout << "PRIME_MPI - Master process:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

int prime_number ( int n, int id, int p )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_NUMBER returns the number of primes between 1 and N.
//
//  Discussion:
//
//    In order to divide the work up evenly among P processors, processor
//    ID starts at 2+ID and skips by P.
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
//    Input, int ID, the ID of this process,
//    between 0 and P-1.
//
//    Input, int P, the number of processes.
//
//    Output, int PRIME_NUMBER, the number of prime numbers up to N.
//
{
  int i;
  int j;
  int prime;
  int total;

  total = 0;

  for ( i = 2 + id; i <= n; i = i + p )
  {
    prime = 1;
    for ( j = 2; j < i; j++ )
    {
      if ( ( i % j ) == 0 )
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
