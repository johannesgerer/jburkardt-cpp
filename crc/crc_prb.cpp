# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "crc.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CRC_PRB.
//
//  Discussion:
//
//    CRC_PRB calls the CRC test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CRC_PRB\n";
  cout << "  C++ version";
  cout << "  Test the CRC library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CRC_PRB\n";
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
//    TEST01 tests MAKE_CRC_TABLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST01\n";
  cout << "  MAKE_CRC_TABLE sets up the 256 entry CRC table.\n";
  cout << "\n";
  cout << "  Each entry is an unsigned long integer.\n";

  make_crc_table();

  print_crc_table();

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CRC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  unsigned char c;
  int len = 10;
  unsigned long my_crc;
  unsigned char s[11] = "Dr. Crypto";
//
//  How do I do the CRC for, say, an integer vector?
//
  cout << "\n";
  cout << "TEST02\n";
  cout << "  CRC computes the CRC of a string of bytes.\n";
  cout << "\n";
  cout << "  The original string: \"" << (char*) s << "\"\n";

  my_crc = crc ( s, len );
  cout << "\n";
  cout << "  For the original string,\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  c = s[3];
  s[3] = s[4];
  s[4] = c;

  cout << "\n";
  cout << "  The modified string: \"" << (char*) s << "\"\n";

  my_crc = crc ( s, len );
  cout << "\n";
  cout << "  After swapping C[3] and C[4],\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  c = s[3];
  s[3] = s[4];
  s[4] = c;

  s[5] = s[5] + 1;

  my_crc = crc ( s, len );

  cout << "\n";
  cout << "  After incrementing C[5] by 1:\n";
  cout << "  The modified string: \"" << (char*) s << "\"\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  my_crc = crc ( s, len-1 );

  cout << "\n";
  cout << "  Using just the first LEN-1 values:\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests UPDATE_CRC_S.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  int len = 11;
  unsigned long my_crc;
  unsigned char s[11] = "Dr. Crypto";

  cout << "\n";
  cout << "TEST03\n";
  cout << "  UPDATE_CRC_S computes the running CRC by processing\n";
  cout << "  a part of the string at a time.\n";
  cout << "  This should be the same as computing the full CRC.\n";
  cout << "\n";
  cout << "  The original string: \"" << (char*) s<< "\"\n";

  my_crc = crc ( s, len );
  cout << "\n";
  cout << "  One step CRC:\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  my_crc = 0xffffffffL;

  my_crc = update_crc_s ( my_crc, &s[0], 5 );
  my_crc = update_crc_s ( my_crc, &s[5], 6 );

  my_crc = my_crc ^ 0xffffffffL;

  cout << "\n";
  cout << "  Incremental CRC:\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests UPDATE_CRC_C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2012
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int len = 11;
  unsigned long my_crc;
  unsigned char s[11] = "Dr. Crypto";

  cout << "\n";
  cout << "TEST04\n";
  cout << "  UPDATE_CRC_C computes the running CRC one byte at a time.\n";
  cout << "  This should be the same as computing the full CRC.\n";
  cout << "\n";
  cout << "  The original string: \"" << (char*) s<< "\"\n";

  my_crc = crc ( s, len );
  cout << "\n";
  cout << "  One step CRC:\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  my_crc = 0xffffffffL;
  for ( i = 0; i < len; i++ )
  {
    my_crc = update_crc_c ( my_crc, s[i] );
  }
  my_crc = my_crc ^ 0xffffffffL;

  cout << "\n";
  cout << "  CRC computed one byte at a time:\n";
  cout << "  CRC = " << dec << my_crc << " (dec)\n";
  cout << "      = " << hex << my_crc << " (hex)\n";

  return;
}
