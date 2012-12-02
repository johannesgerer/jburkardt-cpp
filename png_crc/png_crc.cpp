# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>
# include <ctime>

using namespace std;

int main ( long argc, char *argv[] );

unsigned long crc ( unsigned char *buf, int len );
int handle ( char file_in_name[] );
void make_crc_table ( void );
void print_crc_table ( void );
unsigned long update_crc_c ( unsigned long crc, unsigned char c );
void timestamp ( void );

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

int main ( long argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    PNG_CRC computes the PNG CRC of a file. 
//
//  Usage:
//
//    png_crc files
//
//  Modified:
//
//    21 January 2003
//
//  Reference:
//
//    David Duce, editor,
//    Portable Network Graphics (PNG) Specification,
//    Second Edition.
//    Technical Report REC-PNG-20031110,
//    World Wide Web Consortium, 2003.
//
{
  char file_in_name[80];
  int i;
  bool VERBOSE = true;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "PNG_CRC:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compute the PNG CRC of a file.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "PNG_CRC:\n";
    cout << "  Please enter the name of a file to be indexed.\n";

    cin.getline ( file_in_name, sizeof ( file_in_name ) );

    handle ( file_in_name );
  }
//
//  Otherwise, get the file(s) from the argument list. 
//
  else 
  {

    for ( i = 1 ; i < argc ; ++i ) 
    {
      handle ( argv[i] );
    }

  } 
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "PNG_CRC:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

int handle ( char file_in_name[] )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE computes the PNG CRC of a single file.
//
//  Modified:
//
//    21 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char FILE_IN_NAME[], the name of the file.
//
{
  char c;
  unsigned char c2;
  ifstream file_in;
  unsigned long file_in_crc;
  bool VERBOSE = true;
//
//  Initialize the CRC.
//
  file_in_crc = 0xffffffffL;
//
//  Open the file.
//
  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "PNG_CRC - Fatal error!\n";
    cout << "  Cannot open the input file \"" << file_in_name << "\"\n";
    return 1;
  }

  while ( 1 )
  {
    file_in.read ( &c, 1 );
    c2 = ( unsigned char ) c;

    if ( file_in.eof ( ) )
    {
      break;
    }

    file_in_crc = update_crc_c ( file_in_crc, c2 );
  }

  file_in.close ( );

  cout << file_in_name << " " << hex << file_in_crc << dec << "\n";

  return 0;
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
//  Reference:
//
//    G Randers-Pehrson, et al,
//    PNG (Portable Network Graphics) Specification,
//    Version 1.2, July 1999.
//
//  Modified:
//
//    09 March 2002
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

void timestamp ( void )

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
//  Modified:
//
//    03 October 2003
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
//  Reference:
//
//    G Randers-Pehrson, et al,
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
