# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "crc.hpp"

//
//  Global, unsigned long crc_table[256], the table of cyclic 
//  redundancy checksums for all 8-bit messages.
//
//  Global, int crc_table_computed, a flag recording whether 
//  crc_table has been computed.
//

unsigned long crc_table[256];
bool crc_table_computed = false;

//****************************************************************************80

unsigned long crc ( unsigned char* buf, int len )

//****************************************************************************80
//
//  Purpose:
//
//    CRC returns the CRC of the bytes in BUF[0...LEN-1].
//
//  Discussion:
//
//    Recall that ^ is the bitwise XOR operator and that
//    0xNNN introduces a hexadecimal constant..
//
//  Modified:
//
//    09 March 2002
//
//  Reference:
//
//    Glenn Randers-Pehrson, et al,
//    PNG (Portable Network Graphics) Specification,
//    Version 1.2, July 1999.
//
//  Parameters:
//
//    Input, unsigned char* BUF, a string whose CRC is to be computed.
//
//    Input, int LEN, the number of characters in the string.
//
//    Output, unsigned long CRC, the CRC of the string.
//
{
  return update_crc_s ( 0xffffffffL, buf, len ) ^ 0xffffffffL;
}
//****************************************************************************80

void make_crc_table ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAKE_CRC_TABLE makes the table for a fast CRC computation.
//
//  Discussion:
//
//    Recall that & is the bitwise AND operator, ^ is the bitwise XOR operator,
//    >> is the bitwise right shift operator, and 0xNNN introduces a 
//    hexadecimal constant.
//
//  Modified:
//
//    09 March 2002
//
//  Reference:
//
//    Glenn Randers-Pehrson, et al,
//    PNG (Portable Network Graphics) Specification,
//    Version 1.2, July 1999.
//
{
  unsigned long c;
  int k;
  int n;

  for ( n = 0; n < 256; n++ )
  {
    c = ( unsigned long ) n;
    for ( k = 0; k < 8; k++ )
    {
      if ( c & 1 )
      {
        c = 0xedb88320L ^ ( c >> 1 );
      } 
      else
      {
        c = c >> 1;
      }
    }
    crc_table[n] = c;
  }
  crc_table_computed = true;

  return;
}
//****************************************************************************80

void print_crc_table ( )

//****************************************************************************80
//
//  Purpose:
//
//    PRINT_CRC_TABLE prints the CRC table.
//
//  Discussion:
//
//    Thanks to Saasha Metsarantala for pointing out that the "dec" operator 
//    is required to force decimal output of an integer, 24 July 2012.
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
  int n;

  if ( !crc_table_computed )
  {
    make_crc_table ( );
  }

  cout << "\n";
  cout << "CRC_TABLE\n";
  cout << "  N  dec(CRC(N))  hex(CRC(N))\n";
  cout << "\n";

  for ( n = 0; n < 256; n++ )
  {
    cout << "  " 
      << dec << setw(3)  << n << "  " 
      << dec << setw(12) << crc_table[n] << "  "
      << hex << setw(8)  << crc_table[n] << "\n";
  }
  return;
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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2003
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
//****************************************************************************80

unsigned long update_crc_c ( unsigned long crc, unsigned char c )

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE_CRC_C updates a running CRC with one more character.
//
//  Discussion:
//
//    The value of CRC should have been initialized to all 1's, and the
//    transmitted value is the 1's complement of the final running CRC.
//
//    Recall that & is the bitwise AND operator, ^ is the bitwise XOR operator,
//    and >> is the bitwise right shift operator.
//
//  Modified:
//
//    29 December 2007
//
//  Reference:
//
//    Glenn Randers-Pehrson, et al,
//    PNG (Portable Network Graphics) Specification,
//    Version 1.2, July 1999.
//
//  Parameters:
//
//    Input, unsigned long CRC, the current value of the cyclic redundancy
//    checksum.
//
//    Input, unsigned char C, the next character to be processed.
//
//    Output, unsigned long UPDATE_CRC_C, the updated CRC.
//
{
  unsigned long crc2 = crc;

  if ( !crc_table_computed )
  {
    make_crc_table ( );
  }

  crc2 = crc_table[ ( crc2 ^ c ) & 0xff ] ^ ( crc2 >> 8 );
  
  return crc2;
}
//****************************************************************************80

unsigned long update_crc_s ( unsigned long crc, unsigned char* buf, 
  int len )

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE_CRC_S updates a running CRC with a new string of characters.
//
//  Discussion:
//
//    The value of CRC should have been initialized to all 1's, and the
//    transmitted value is the 1's complement of the final running CRC.
//
//    Recall that & is the bitwise AND operator, ^ is the bitwise XOR operator,
//    and >> is the bitwise right shift operator.
//
//  Modified:
//
//    29 December 2007
//
//  Reference:
//
//    Glenn Randers-Pehrson, et al,
//    PNG (Portable Network Graphics) Specification,
//    Version 1.2, July 1999.
//
//  Parameters:
//
//    Input, unsigned long CRC, the current value of the cyclic redundancy
//    checksum.
//
//    Input, unsigned char* BUF, the "next" string of characters to be
//    processed.
//
//    Input, int LEN, the number of characters in BUF.
//
//    Output, unsigned long UPDATE_CRC_S, the updated CRC.
//
{
  unsigned long c = crc;
  int n;

  if ( !crc_table_computed )
  {
    make_crc_table ( );
  }

  for ( n = 0; n < len; n++ )
  {
    c = crc_table[ ( c ^ buf[n] ) & 0xff ] ^ ( c >> 8 );
  }
  
  return c;
}
