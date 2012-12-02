# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "prime_serial.hpp"

int main ( );
void prime_number_sweep ( int n_lo, int n_hi, int n_factor );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRIME_SERIAL_PRB.
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

  timestamp ( );

  cout << "\n";
  cout << "PRIME_SERIAL_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PRIME_SERIAL library.\n";

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
  cout << "PRIME_SERIAL_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

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
  double ctime;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Call PRIME_NUMBER to count the primes from 1 to N.\n";
  cout << "\n";
  cout << "         N        Pi          Time\n";
  cout << "\n";

  n = n_lo;

  while ( n <= n_hi )
  {
    ctime = cpu_time ( );

    primes = prime_number ( n );

    ctime = cpu_time ( ) - ctime;

    cout << "  " << setw(8) << n
         << "  " << setw(8) << primes
         << "  " << setw(14) << ctime << "\n";

    n = n * n_factor;
  }

  return;
}


