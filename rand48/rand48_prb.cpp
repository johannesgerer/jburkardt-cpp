# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( void );
void test06 ( void );
void test07 ( void );
void test08 ( void );
void test09 ( void );
void timestamp ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RAND48_PRB.
//
//  Discussion:
//
//    RAND48_PRB tests the RAND48 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "RAND48_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the RAND48 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "RAND48_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests MRAND48 + SRAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  long int seed;
  long int value;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  MRAND48 returns signed long integers in [-2^31,+2^31].\n";
  cout << "  SRAND48 is used to initialize the seed (but only 32 bits).\n";

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests MRAND48 + SEED48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  MRAND48 returns signed long integers in [-2^31,+2^31].\n";
  cout << "  SEED48 is used to initialize the seed (all 48 bits).\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = mrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests JRAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  JRAND48 returns signed long integers in [-2^31,+2^31].\n";
  cout << "  The 48 bit seed is an explicit argument, 3 16 bit values.\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = jrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test04 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests LRAND48 + SRAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  long int seed;
  long int value;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  LRAND48 returns unsigned long integers in [0,+2^31].\n";
  cout << "  SRAND48 is used to initialize the seed (32 bits only).\n";

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test05 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests LRAND48 + SEED48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  LRAND48 returns unsigned long integers in [0,+2^31].\n";
  cout << "  SEED48 is used to initialize the seed (all 48 bits).\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );


  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = lrand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test06 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests NRAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  long int value;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  NRAND48 returns nonnegative long integers in [0,+2^31].\n";
  cout << "  The 48 bit seed is an explicit argument of 3 16 bit values.\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";


  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";


  for ( i = 1; i <= 10; i++ )
  {
    value = nrand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test07 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests DRAND48 + SRAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  long int seed;
  double value;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  DRAND48 returns double precision floating point values in [0.0,1.0].\n";
  cout << "  SRAND48 is used to initialize the seed (32 bits only).\n";

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789L;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  srand48 ( seed );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test08 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests DRAND48 + SEED48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  double value;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  DRAND48 returns double precision real values in [0.0,1.0].\n";
  cout << "  The 48 bit seed is an explicit argument of 3 16 bit values.\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  seed48 ( seedvec );

  for ( i = 1; i <= 10; i++ )
  {
    value = drand48 ( );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
}
//****************************************************************************80

void test09 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests ERAND48.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int modulus = 65536;
  long long int seed;
  unsigned short int seedvec[3];
  double value;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  ERAND48 returns double precision real values in [0.0,1.0].\n";
  cout << "  The 48 bit seed is an explicit argument of 3 16 bit values.\n";

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 987654321LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }

  seed = 123456789LL;
  seedvec[0] = seed % modulus;
  seed = seed / modulus;
  seedvec[1] = seed % modulus;
  seed = seed / modulus;
  seedvec[2] = seed % modulus;
  seed = seed / modulus;

  cout << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "  The seed vector is " << seedvec[0] 
       << ", " << seedvec[1] 
       << ", " << seedvec[2] << "\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    value = erand48 ( seedvec );

    cout << "  " <<  setw(6) << i
         << "  " << setw(12) << value << "\n";
  }
  return;
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
