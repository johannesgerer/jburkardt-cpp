# include <cstdlib>
# include <iostream>

using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRIME_SUM.
//
//  Discussion:
//
//    PRIME_SUM adds up the prime numbers from 2 to N.
//
//    Invoke the program with a command like 
//
//      prime_sum n
//
//    The sum is printed out.
//
//  Example:
//
//    prime_sum   10  =      17 = 2 + 3 + 5 + 7
//    prime_sum  100  =    1060
//    prime_sum 1000  =   76127
//    prime_sum 10000 = 5736396
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  bool prime;
  int total;
//
//  If called with no arguments, use N = 100.
//
  if ( argc < 1 )
  {
    n = 100;
  }
  else
  {
    *argv++;
    n = atoi ( *argv );
  }

  total = 0;
//
//  Consider each integer I from 2 to N.
//
  for ( i = 2; i <= n; i++ )
  {
//
//  Determine if I is prime.
//
    prime = true;

    for ( j = 2; j < i; j++ )
    {
      if ( i % j == 0 )
      {
        prime = false;
        break;
      }
    }
//
//  If prime, add to sum.
//
    if ( prime )
    {
      total = total + i;
    }
  }
//
//  Print the sum.
//
  cout << "PRIME_SUM(2:" << n << ") = " << total << "\n";

  return 0;
}
