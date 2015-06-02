# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <stdint.h>

# include "ziggurat_inline.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( int sample_num );
void test06 ( int sample_num );
void test07 ( int sample_num );
void test08 ( int sample_num );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ZIGGURAT_INLINE_PRB.
//
//  Discussion:
//
//    ZIGGURAT_ORIGINAL_PRB tests the ZIGGURAT_INLINE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int sample_num = 1000000;

  timestamp ( );
  cout << "\n";
  cout << "ZIGGURAT_INLINE_PRB\n";
  cout << "  C++ version:\n";
  cout << "  Test the ZIGGURAT_INLINE library.\n";
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
//  List 10 successive values of the unsigned int 32 bit generators.
//
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ZIGGURAT_INLINE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests KISS_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  uint32_t jcong_value;
  uint32_t jsr_value;
  uint32_t value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  KISS_VALUE is an inline generator of pseudorandom uniformly\n";
  cout << "  distributed uint32_t's.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      jsr_value = 123456789;
      jcong_value = 234567891;
      w_value = 345678912;
      z_value = 456789123;
    }
    else
    {
      jsr_value = 987654321;
      jcong_value = 198765432;
      w_value = 219876543;
      z_value = 321987654;
    }

    cout << "\n";
    cout << "  Call ZIGSET with these seeds:\n";
    cout << "\n";
    cout << "  jsr_value =   " << jsr_value << "\n";
    cout << "  jcong_value = " << jcong_value << "\n";
    cout << "  w_value =     " << w_value << "\n";
    cout << "  z_value =     " << z_value << "\n";
    cout << "\n";
//
//  Call ZIGSET to set the seed.
//
    zigset ( jsr_value, jcong_value, w_value, z_value );

    for ( i = 1; i <= 10; i++ )
    {
      value = kiss_value ( );
      cout << "  " << setw(6)  << i
           << "  " << setw(12) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests R4_UNI_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  uint32_t jcong_value;
  uint32_t jsr_value;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R4_UNI_VALUE is an inline generator of pseudorandom uniformly\n";
  cout << "  distributed floats between 0 and 1.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      jsr_value = 123456789;
      jcong_value = 234567891;
      w_value = 345678912;
      z_value = 456789123;
    }
    else
    {
      jsr_value = 987654321;
      jcong_value = 198765432;
      w_value = 219876543;
      z_value = 321987654;
    }

    cout << "\n";
    cout << "  Call ZIGSET with these seeds:\n";
    cout << "\n";
    cout << "  jsr_value =   " << jsr_value << "\n";
    cout << "  jcong_value = " << jcong_value << "\n";
    cout << "  w_value =     " << w_value << "\n";
    cout << "  z_value =     " << z_value << "\n";
    cout << "\n";
//
//  Call ZIGSET to set the seed.
//
    zigset ( jsr_value, jcong_value, w_value, z_value );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_uni_value ( );
      cout << "  " << setw(6)  << i
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R4_NOR_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  uint32_t jcong_value;
  uint32_t jsr_value;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  R4_NOR_VALUE is an inline generator of pseudorandom normally\n";
  cout << "  distributed floats.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      jsr_value = 123456789;
      jcong_value = 234567891;
      w_value = 345678912;
      z_value = 456789123;
    }
    else
    {
      jsr_value = 987654321;
      jcong_value = 198765432;
      w_value = 219876543;
      z_value = 321987654;
    }

    cout << "\n";
    cout << "  Call ZIGSET with these seeds:\n";
    cout << "\n";
    cout << "  jsr_value =   " << jsr_value << "\n";
    cout << "  jcong_value = " << jcong_value << "\n";
    cout << "  w_value =     " << w_value << "\n";
    cout << "  z_value =     " << z_value << "\n";
    cout << "\n";
//
//  Call ZIGSET to set the seed.
//
    zigset ( jsr_value, jcong_value, w_value, z_value );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_nor_value ( );
      cout << "  " << setw(6)  << i
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests R4_EXP_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  uint32_t jcong_value;
  uint32_t jsr_value;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  R4_EXP_VALUE is an inline generator of pseudorandom exponentially\n";
  cout << "  distributed floats.\n";

  for ( j = 0; j < 3; j++ )
  {
    if ( ( j % 2 ) == 0 )
    {
      jsr_value = 123456789;
      jcong_value = 234567891;
      w_value = 345678912;
      z_value = 456789123;
    }
    else
    {
      jsr_value = 987654321;
      jcong_value = 198765432;
      w_value = 219876543;
      z_value = 321987654;
    }

    cout << "\n";
    cout << "  Call ZIGSET with these seeds:\n";
    cout << "\n";
    cout << "  jsr_value =   " << jsr_value << "\n";
    cout << "  jcong_value = " << jcong_value << "\n";
    cout << "  w_value =     " << w_value << "\n";
    cout << "  z_value =     " << z_value << "\n";
    cout << "\n";
//
//  Call ZIGSET to set the seed.
//
    zigset ( jsr_value, jcong_value, w_value, z_value );

    for ( i = 1; i <= 10; i++ )
    {
      value = r4_exp_value ( );
      cout << "  " << setw(6)  << i
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
//    TEST05 times KISS_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  uint32_t jcong_value;
  uint32_t jsr_value;
  int sample;
  uint32_t value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Measure the time it takes KISS_VALUE to generate\n";
  cout << "  " << sample_num << " uint32_t values.\n";
//
//  Call ZIGSET to set the seed.
//
  jsr_value = 123456789;
  jcong_value = 234567891;
  w_value = 345678912;
  z_value = 456789123;
  zigset ( jsr_value, jcong_value, w_value, z_value );

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = kiss_value ( );
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
//    TEST06 times R4_UNI_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  uint32_t jcong_value;
  uint32_t jsr_value;
  int sample;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Measure the time it takes R4_UNI_VALUE to generate\n";
  cout << "  " << sample_num << " uniformly random float values.\n";
//
//  Call ZIGSET to set the seed.
//
  jsr_value = 123456789;
  jcong_value = 234567891;
  w_value = 345678912;
  z_value = 456789123;
  zigset ( jsr_value, jcong_value, w_value, z_value );

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_uni_value ( );
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
//    TEST07 times R8_NOR_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  uint32_t jcong_value;
  uint32_t jsr_value;
  int sample;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Measure the time it takes R8_NOR_VALUE to generate\n";
  cout << "  " << sample_num << " normal random float values.\n";
//
//  Call ZIGSET to set the seed.
//
  jsr_value = 123456789;
  jcong_value = 234567891;
  w_value = 345678912;
  z_value = 456789123;
  zigset ( jsr_value, jcong_value, w_value, z_value );

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_nor_value ( );
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
//    TEST08 times R4_EXP_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ctime;
  uint32_t jcong_value;
  uint32_t jsr_value;
  int sample;
  float value;
  uint32_t w_value;
  uint32_t z_value;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Measure the time it takes R4_EXP_VALUE to generate\n";
  cout << "  " << sample_num << " exponential random float values.\n";
//
//  Call ZIGSET to set the seed.
//
  jsr_value = 123456789;
  jcong_value = 234567891;
  w_value = 345678912;
  z_value = 456789123;
  zigset ( jsr_value, jcong_value, w_value, z_value );

  ctime = cpu_time ( );

  for ( sample = 0; sample < sample_num; sample++ )
  {
    value = r4_exp_value ( );
  }
  ctime = cpu_time ( ) - ctime;

  cout << "\n";
  cout << "  " << ctime << " seconds.\n";

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests CONG_SEEDED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  uint32_t jcong_new;
  uint32_t jcong_in;
  uint32_t jcong_old;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  CONG_SEEDED is a generator of pseudorandom uniformly\n";
  cout << "  distributed unsigned 32 bit integers.\n";
  cout << "\n";
  cout << "    Input Seed   Output Seed  Output Value\n";
  cout << "\n";

  jcong_new = 234567891;

  for ( j = 1; j <= 10; j++ )
  {
    jcong_old = jcong_new;
    jcong_in = jcong_new;
    jcong_new = cong_seeded ( jcong_in );
    cout << "  " << setw(12) << jcong_old
         << "  " << setw(12) << jcong_in
         << "  " << setw(12) << jcong_new << "\n";
  }

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests KISS_SEEDED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  uint32_t jcong_in;
  uint32_t jcong_old;
  uint32_t jsr_in;
  uint32_t jsr_old;
  uint32_t w_in;
  uint32_t w_old;
  uint32_t value;
  uint32_t z_in;
  uint32_t z_old;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  KISS_SEEDED is a generator of pseudorandom uniformly\n";
  cout << "  distributed unsigned 32 bit integers.\n";
  cout << "\n";
  cout << "              JCONG           JSR             W             Z         Value\n";
  cout << "\n";

  jcong_in = 234567891;
  jsr_in = 123456789;
  w_in = 345678912;
  z_in = 456789123;

  for ( j = 1; j <= 10; j++ )
  {
    jcong_old = jcong_in;
    jsr_old = jsr_in;
    w_old = w_in;
    z_old = z_in;
    value = kiss_seeded ( jcong_in, jsr_in, w_in, z_in );
    cout << "  In "
         << "  " << setw(12) << jcong_old
         << "  " << setw(12) << jsr_old
         << "  " << setw(12) << w_old
         << "  " << setw(12) << z_old << "\n";
    cout << "  Out"
         << "  " << setw(12) << jcong_in
         << "  " << setw(12) << jsr_in
         << "  " << setw(12) << w_in
         << "  " << setw(12) << z_in
         << "  " << setw(12) << value << "\n";
  }

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests MWC_SEEDED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  uint32_t w_in;
  uint32_t w_old;
  uint32_t value;
  uint32_t z_in;
  uint32_t z_old;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  MWC_SEEDED is a generator of pseudorandom uniformly\n";
  cout << "  distributed unsigned 32 bit integers.\n";
  cout << "\n";
  cout << "       Input W       Input Z      Output W      Output Z  Output Value\n";
  cout << "\n";

  w_in = 345678912;
  z_in = 456789123;

  for ( j = 1; j <= 10; j++ )
  {
    w_old = w_in;
    z_old = z_in;
    value = mwc_seeded ( w_in, z_in );
    cout << "  " << setw(12) << w_old
         << "  " << setw(12) << z_old
         << "  " << setw(12) << w_in
         << "  " << setw(12) << z_in
         << "  " << setw(12) << value << "\n";
  }

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests SHR3_SEEDED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2013
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  uint32_t jsr_new;
  uint32_t jsr_in;
  uint32_t jsr_old;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  SHR3_SEEDED is a generator of pseudorandom uniformly\n";
  cout << "  distributed unsigned 32 bit integers.\n";
  cout << "\n";
  cout << "    Input Seed   Output Seed  Output Value\n";
  cout << "\n";

  jsr_new = 123456789;

  for ( j = 1; j <= 10; j++ )
  {
    jsr_old = jsr_new;
    jsr_in = jsr_new;
    jsr_new = shr3_seeded ( jsr_in );
    cout << "  " << setw(12) << jsr_old
         << "  " << setw(12) << jsr_in
         << "  " << setw(12) << jsr_new << "\n";
  }

  return;
}
