# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "tga_io.hpp"

int main ( long argc, char *argv[] );
void test01 ( char *file_name );
void test02 ( char *file_input_name, char *file_output_name );
void timestamp ( void );

//****************************************************************************80

int main ( long argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TGA_IO_PRB.
//
//  Discussion:
//
//    TGA_IO_PRB calls the TGA_IO built-in test routines.
//
//  Modified:
//
//    02 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  bool tga_byte_swap;
  char *file_name1 = "earth.tga";
  char *file_name1_copy = "earth_copy.tga";
  char *file_name2 = "shuttle.tga";
  char *file_name2_copy = "shuttle_copy.tga";

  timestamp ( );

  cout << "\n";
  cout << "TGA_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TGA_IO library.\n";
//
//  Set the byte swapping option.
//
  tga_byte_swap = true;
  tga_byte_swap_set ( tga_byte_swap );

  cout << "\n";
  cout << "  The TGA_BYTE_SWAP option being used is " << tga_byte_swap << "\n";

  test01 ( file_name1 );
  test01 ( file_name2 );

  test02 ( file_name1, file_name1_copy );
  
  cout << "\n";
  cout << "TGA_IO_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( char *file_input_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 calls TGA_READ to read a TGA file.
//
//  Modified:
//
//    31 January 2006
//
//  Author:
//
//    John Burkardt
//
{
  bool success;
  tga_header *th;
  tga_color_map *tc;
  tga_data *td;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  TGA_READ can read a TGA file.\n";
  cout << "  Here, we read the file \"" << file_input_name << "\".\n" << flush;

  success = tga_read ( file_input_name, &th, &tc, &td );

  if ( !success ) 
  {
    cout << "\n";
    cout << "  TGA READ failed.\n";
    exit ( 1 );
  }

  tga_header_print ( th );

  delete [] th;
  delete [] tc;
  delete [] td;

  return;
}
//****************************************************************************80

void test02 ( char *file_input_name, char *file_output_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 calls TGA_WRITE to write a TGA file.
//
//  Modified:
//
//    02 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  ofstream file_output_unit;
  bool success;
  tga_header *th;
  tga_color_map *tc;
  tga_data *td;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  TGA_WRITE can write a TGA file.\n";
  cout << "\n";
  cout << "  Call TGA_READ to read \"" << file_input_name << "\".\n" << flush;

  success = tga_read ( file_input_name, &th, &tc, &td );

  if ( !success ) 
  {
    cout << "\n";
    cout << "  TGA READ failed.\n";
    exit ( 1 );
  }

  tga_header_print ( th );

  cout << "\n";
  cout << "  Now write the data to \"" << file_output_name << "\".\n" << flush;
  
  file_output_unit.open ( file_output_name, ios::binary );

  if ( !file_output_unit.is_open ( ) )
  {
    cout << "\n";
    cout << "TEST02 - Fatal error!\n";
    cout << "  Could not open the file \"" << file_output_name << "\".\n";
    return;
  }

  tga_header_write ( file_output_unit, *th );

  tga_color_map_write ( file_output_unit, *th, *tc );

  tga_data_write ( file_output_unit, *th, *tc, *td );

  file_output_unit.close ( );

  delete [] th;
  delete [] tc;
  delete [] td;

  return;
}
