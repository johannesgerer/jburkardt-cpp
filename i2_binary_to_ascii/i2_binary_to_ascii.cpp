# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>

using namespace std;

int main ( );
void i2_binary_to_ascii_convert ( string input_filename, string output_filename );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    I2_BINARY_TO_ASCII converts a file of I2 data from binary to ASCII format.
//
//  Discussion:
//
//    An "I2" is an integer represented by two bytes, or 16 bits of data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  string input_filename = "sealevel.bin";
  string output_filename = "sealevel.txt";

  timestamp ( );
  cout << "\n";
  cout << "I2_BINARY_TO_ASCII:\n";
  cout << "  C++ version\n";
  cout << "  Given a file of 2 byte integers stored in binary format,\n";
  cout << "  create an equivalent ASCII text file.\n";

  i2_binary_to_ascii_convert ( input_filename, output_filename );
//
//  Terminate.
//
  cout << "\n";
  cout << "I2_BINARY_TO_ASCII:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i2_binary_to_ascii_convert ( string input_filename, string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    I2_BINARY_TO_ASCII_CONVERT converts I2 data from binary to text form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file containing
//    the binary data.
//
//    Input, string OUTPUT_FILENAME, the name of the output file, to be
//    created and filled with the corresponding ASCII data.
//
{
  bool byte_swap = true;
  char c;
  int i;
  ifstream input;
  int j;
  int k;
  char *mem_char;
  unsigned short int *mem_unsigned_short_int;
  ofstream output;
  ifstream::pos_type size;
//
//  Read the information into memory.
//
  input.open ( input_filename.c_str ( ), ios::in|ios::binary|ios::ate );
  size = input.tellg ( );
  mem_char = new char[size];
  input.seekg ( 0, ios::beg );
  input.read ( mem_char, size );
  input.close ( );
//
//  If requested, swap bytes.
//
  if ( byte_swap )
  {
    k = 0;
    for ( j = 0; j < 121; j++ )
    {
      for ( i = 0; i < 301; i++ )
      {
        c = mem_char[k];
        mem_char[k] = mem_char[k+1];
        mem_char[k+1] = c;
        k = k + 2;
      }
    }
  }
//
//  Write the information to an ASCII file.
//
  mem_unsigned_short_int = ( unsigned short int * ) mem_char;

  output.open ( output_filename.c_str ( ) );

  k = 0;
  for ( j = 0; j < 121; j++ )
  {
    for ( i = 0; i < 301; i++ )
    {
      output << mem_unsigned_short_int[k] << "\n";
      k = k + 1;
    }
  }
  output.close ( );

  delete [] mem_char;

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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

