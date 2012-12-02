# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

# include "ziggurat.hpp"

using namespace std;

int main ( void );
double cpu_time ( void );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void test04 ( void );
void test05 ( int sample_num );
void test06 ( int sample_num );
void test07 ( int sample_num );
void test08 ( int sample_num );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ZIGGURAT_PRB.
//
//  Discussion:
//
//    ZIGGURAT_PRB calls sample problems for the ZIGGURAT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  int sample_num = 1000000;

  timestamp ( );

  cout << "\n";
  cout << "ZIGGURAT_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the ZIGGURAT library.\n";
//
//  Make sure that SEED controls the sequence, and can restart it.
//
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Measure the time it takes to generate 10,000 variables.
//
  test05 ( sample_num );
  test06 ( sample_num );
  test07 ( sample_num );
  test08 ( sample_num );
//
//  Terminate.
//
  cout << "\n";
  cout << "ZIGGURAT_PRB\n";
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
//    TEST01 tests SHR3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  unsigned long int seed;
  unsigned long int value;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  SHR3 returns pseudorandom uniformly distributed\n";
  cout << "  unsigned long integers.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    cout << "\n";
    cout << "  " << setw(6)  << 0
         << "  " << setw(12) << ( long int ) seed
         << "  " << setw(12) << seed << "\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      value = shr3 ( &seed );

      cout << "  " << setw(6)  << i
           << "  " << setw(12) << ( long int ) seed
           << "  " << setw(12) << seed
           << "  " << setw(12) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R4_UNI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  unsigned long int seed;
  float value;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R4_UNI returns pseudorandom uniformly distributed\n";
  cout << "  floats between 0 and 1.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    cout << "\n";
    cout << "  " << setw(6)  << 0
         << "  " << setw(12) << ( long int ) seed
         << "  " << setw(12) << seed << "\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_uni ( &seed );

      cout << "  " << setw(6)  << i
           << "  " << setw(12) << ( long int ) seed
           << "  " << setw(12) << seed
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R4_NOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  float fn[128];
  int i;
  int j;
  int kn[128];
  unsigned long int seed;
  float value;
  float wn[128];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  R4_NOR returns pseudorandom normally distributed\n";
  cout << "  real numbers between 0 and 1.\n";

  r4_nor_setup ( kn, fn, wn );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    cout << "\n";
    cout << "  " << setw(6)  << 0
         << "  " << setw(12) << ( long int ) seed
         << "  " << setw(12) << seed << "\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_nor ( &seed, kn, fn, wn );

      cout << "  " << setw(6)  << i
           << "  " << setw(12) << ( long int ) seed
           << "  " << setw(12) << seed
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test04 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests R4_EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  float fe[256];
  int i;
  int j;
  int ke[256];
  unsigned long seed;
  float value;
  float we[256];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  R4_EXP returns pseudorandom exponentially distributed\n";
  cout << "  real numbers between 0 and 1.\n";

  r4_exp_setup ( ke, fe, we );

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      seed = 123456789;
    }
    else
    {
      seed = 987654321;
    }

    cout << "\n";
    cout << "  " << setw(6)  << 0
         << "  " << setw(12) << ( long int ) seed
         << "  " << setw(12) << seed << "\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_exp ( &seed, ke, fe, we );

      cout << "  " << setw(6)  << i
           << "  " << setw(12) << ( long int ) seed
           << "  " << setw(12) << seed
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test05 ( int sample_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 times SHR3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  int sample;
  unsigned long int seed;
  unsigned long int value;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Measure the time it takes SHR3 to generate\n";
  cout << "  " << sample_num << " unsigned long int values.\n";

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = shr3 ( &seed );
  }
  ctime = cpu_time ( ) - ctime;

  cout << "\n";
  cout << "  " << ctime << " seconds.\n";

  return;
}
//****************************************************************************80

void test06 ( int sample_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 times R4_UNI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  int sample;
  unsigned long int seed;
  float value;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Measure the time it takes R4_UNI to generate\n";
  cout << "  " << sample_num << " uniformly random float values.\n";

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_uni ( &seed );
  }
  ctime = cpu_time ( ) - ctime;

  cout << "\n";
  cout << "  " << ctime << " seconds.\n";

  return;
}
//****************************************************************************80

void test07 ( int sample_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 times R8_NOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  float fn[128];
  int kn[128];
  int sample;
  unsigned long seed;
  float value;
  float wn[129];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Measure the time it takes R8_NOR to generate\n";
  cout << "  " << sample_num << " normal random float values.\n";

  r4_nor_setup ( kn, fn, wn );

  seed = 123456789;

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_nor ( &seed, kn, fn, wn );
  }
  ctime = cpu_time ( ) - ctime;

  cout << "\n";
  cout << "  " << ctime << " seconds.\n";

  return;
}
//****************************************************************************80

void test08 ( int sample_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 times R4_EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  float fe[256];
  int ke[256];
  int sample;
  unsigned long int seed;
  float value;
  float we[256];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Measure the time it takes R4_EXP to generate\n";
  cout << "  " << sample_num << " exponential random float values.\n";

  r4_exp_setup ( ke, fe, we );

  seed = 123456789;


  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_exp ( &seed, ke, fe, we );
  }
  ctime = cpu_time ( ) - ctime;

  cout << "\n";
  cout << "  " << ctime << " seconds.\n";

  return;
}
//****************************************************************************80

double cpu_time ( void )

//****************************************************************************80
//
//  Purpose:
//
//    CPU_TIME returns the current reading on the CPU clock.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
