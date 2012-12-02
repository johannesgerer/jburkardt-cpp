# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "i8lib.H"

int main ( );
void test015 ( );
void test190 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for I8LIB_PRB.
//
//  Discussion:
//
//    I8LIB_PRB calls the I8LIB tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "I8LIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the I8LIB library.\n";

  test015 ( );
  test190 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "I8LIB_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests I8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  long long int cnk;
  long long int k;
  long long int n;

  cout << "\n";
  cout << "TEST015\n";
  cout << "  I8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "     N     K    CNK\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk = i8_choose ( n, k );

      cout << setw(6) << n   << "  "
           << setw(6) << k   << "  "
           << setw(6) << cnk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test190 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST190 tests I8_XOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 June 2010
//
//  Author:
//
//    John Burkardt
//
{
  long long int i;
  long long int i_lo = 0LL;
  long long int i_hi = 100LL;
  long long int j;
  long long int k;
  long long int l;
  long long int seed;
  int test;
  int test_num = 10;

  seed = 123456789LL;

  cout << "\n";
  cout << "TEST190\n";
  cout << "  I8_XOR returns the bitwise exclusive OR of\n";
  cout << "  two I8's.\n";
  cout << "  The operator ^ should generally be used instead.\n";
  cout << "\n";
  cout << "       I       J  I8_XOR(I,J)\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    i = i8_uniform ( i_lo, i_hi, &seed );
    j = i8_uniform ( i_lo, i_hi, &seed );
    k = i8_xor ( i, j );
    l = i ^ j;

    cout << "  " << setw(6) << i 
         << "  " << setw(6) << j 
         << "  " << setw(6) << k 
         << "  " << setw(6) << l << "\n";
  }

  return;
}
