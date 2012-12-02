# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

int main ( long argc, char *argv[] );
void handle ( char file_in_name[] );
void timestamp ( void );

//****************************************************************************80

int main ( long argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HEXDUMP.
//
//  Discussion:
//
//    HEXDUMP is a simple hexadecimal dump program.  
//
//    HEXDUMP is similar to the UNIX octal dump program OD with the "-h" switch,
//    but automatically includes a listing of printable ASCII characters.
//  
//    The original version of this program is discussed on page 219 
//    of the reference.
//
//  Usage:
//
//    hexdump FILE_NAME > OUTPUT  for command-line use;
//    hexdump                     for interactive use.
//
//  Sample output:
//
//    HEXDUMP: Hexadecimal dump of file: box.3ds.
//
//     Address               Hexadecimal values                  Printable
//     -------  -----------------------------------------------  ----------------
//
//           0  4d 4d 53 02 00 00 02 00 0a 00 00 00 03 00 00 00  MMS.............
//          16  3d 3d 68 01 00 00 3e 3d 0a 00 00 00 03 00 00 00  ==h...>=........
//          32  00 01 0a 00 00 00 00 00 80 3f 00 40 4e 01 00 00  .........?.@N...
//          48  42 6f 78 30 31 00 00 41 42 01 00 00 10 41 68 00  Box01..AB....Ah.
//          64  00 00 08 00 26 76 0a c2 1f 0c a8 c1 00 00 00 00  ....&v..........
//       ...
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    Original C version by Howard Burdick.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Howard Burdick,
//    Digital Imaging, Theory and Applications,
//    McGraw Hill, 1997,
//    ISBN: 0079130593,
//    LC: TA1637.B87.
//
{
  char file_in_name[80];
  int i;
  bool VERBOSE = true;

  if ( VERBOSE )
  {
    timestamp ( );

    cout << "\n";
    cout << "HEXDUMP:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Produce a hexadecimal dump of a file.\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "HEXDUMP:\n";
    cout << "  Please enter the name of a file to be analyzed.\n";

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

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "HEXDUMP:\n";
    cout << "  Normal end of execution.\n";

    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

void handle ( char file_in_name[] )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE creates a hexdump for a single file.
//
//  Discussion:
//
//    Thanks to John McCabe of Rutgers for pointing out and correcting a
//    bug which occurred on some machines if the byte value 0xFF was
//    encountered.
//
//  Modified:
//
//    30 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char FILE_IN_NAME[], the name of the file to handle.
//
{
  long int addr;
  unsigned char buffer[20];
  long int cnt;
  long int cnt2;
  ifstream file_in;
  long n;
//
//  Open the file.
//
  file_in.open ( file_in_name );

  if ( !file_in ) 
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open \"" << file_in_name << "\"\n";
    return;
  }
     
  cout << "\n";
  cout << "Hexdump of \"" << file_in_name << "\":\n";
  cout << "\n";
  cout << 
    "Address               Hexadecimal values                  Printable\n";
  cout << 
    "-------  -----------------------------------------------  -------------\n";
  cout << "\n";
//
//  Dump the file contents.
//
  addr = 0;

  while ( 1 )
  {
    file_in.read ( ( char * ) buffer, 16 );

    cnt = file_in.gcount();

    if ( cnt <= 0 )
    {
      break;
    }
//
//  Print the address in decimal and hexadecimal.
//
    cout << setw(7) << ( int ) addr << "  ";

    addr = addr + 16;
//
//  Print 16 data items, in pairs, in hexadecimal.
//
    cnt2 = 0;
    for ( n = 0; n < 16; n++ )
    {   
      cnt2 = cnt2 + 1;
      if ( cnt2 <= cnt )
      {
        cout << hex << setw(2) << setfill ( '0' ) << ( int ) buffer[n];
      }
      else
      {
        cout << "  ";
      }
      cout << " ";
    }

    cout << setfill ( ' ' );
//
//  Print the printable characters, or a period if unprintable.
//
    cout << " ";
    cnt2 = 0;
    for ( n = 0; n < 16; n++ )
    {
      cnt2 = cnt2 + 1;
      if ( cnt2 <= cnt )
      {
        if ( buffer[n] < 32 || 126 < buffer[n] )
        {
          cout << '.';
        }
        else
        {
          cout << buffer[n];
        }
      }
    }
    cout << "\n";
    cout << dec;

    if ( file_in.eof ( ) )
    {
      break;
    }

  }
//
// Close the file.
//
  file_in.close ( );

  return;
}
//********************************************************************** 

void timestamp ( void )

//**********************************************************************
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
