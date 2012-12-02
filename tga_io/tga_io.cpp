# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "tga_io.hpp"

bool tga_byte_swap = true;

//****************************************************************************80

bool tga_byte_swap_get ( void )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_BYTE_SWAP_GET returns the internal value of TGA_BYTE_SWAP.
// 
//  Modified:
// 
//    26 February 2003
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Output, bool TGA_BYTE_SWAP_GET, the internal value of TGA_BYTE_SWAP.
//
{
  return tga_byte_swap;
}
//****************************************************************************80

void tga_byte_swap_set ( bool value )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_BYTE_SWAP_SET sets the internal value of TGA_BYTE_SWAP.
// 
//  Modified:
// 
//    26 February 2003
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, bool VALUE, the new value of TGA_BYTE_SWAP.
//
{
  tga_byte_swap = value;

  return;
}
//****************************************************************************80

void tga_header_print ( tga_header *th )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_HEADER_PRINT prints the header of a TGA file.
// 
//  Modified:
// 
//    07 April 2005
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, tga_header *TH, a pointer to the TGA header.
//
{
  cout << "\n";
  cout << "  ID_LENGTH = " << ( int ) th->id_length << "\n";
  cout << "  COLOR_MAP_TYPE = " << ( int ) th->color_map_type << "\n";
  cout << "  IMAGE_TYPE = " << ( int ) th->image_type << "\n";
  cout << "  COLOR_MAP_ORIGIN = " << th->color_map_origin << "\n";
  cout << "  COLOR_MAP_LENGTH = " << th->color_map_length << "\n";
  cout << "  COLOR_MAP_ENTRY_SIZE = " 
       << ( int ) th->color_map_entry_size << "\n";
  cout << "  IMAGE_X_ORIGIN = " << th->image_x_origin << "\n";
  cout << "  IMAGE_Y_ORIGIN = " << th->image_y_origin << "\n";
  cout << "  IMAGE_WIDTH = " << th->image_width << "\n";
  cout << "  IMAGE_HEIGHT= " << th->image_height << "\n";
  cout << "  IMAGE_PIXEL_DEPTH = " << ( int ) th->image_pixel_depth << "\n";
  cout << "  IMAGE_DESCRIPTOR = " << ( int ) th->image_descriptor << "\n";
  if ( 0 < th->id_length )
  {
    cout << "  IMAGE_ID = \"" << th->image_id << "\"\n";
  }
  else
  {
    cout << "  IMAGE_ID = \"NULL\"\n";
  }

  return;
}
//****************************************************************************80

tga_header *tga_header_read ( ifstream &file_input )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_HEADER_READ reads the header of a TGA file.
// 
//  Modified:
// 
//    05 April 2005
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Output, tga_header *TGA_HEADER_READ, a pointer to the TGA header.
//
{
  tga_header *th;

  th = new tga_header;

  th->id_length = u_char_read ( file_input );
  th->color_map_type = u_char_read ( file_input );
  th->image_type = u_char_read ( file_input );
  th->color_map_origin = u_short_int_read ( file_input );
  th->color_map_length = u_short_int_read ( file_input );
  th->color_map_entry_size = u_char_read ( file_input );

  th->image_x_origin = u_short_int_read ( file_input );
  th->image_y_origin = u_short_int_read ( file_input );
  th->image_width = u_short_int_read ( file_input );
  th->image_height = u_short_int_read ( file_input );
  th->image_pixel_depth = u_char_read ( file_input );
//
//  IMAGE_DESCRIPTOR includes information about ALPHA channel bits.
//  But if IMAGE_DESCRIPTOR is 0, then for sure there are none.
//
  th->image_descriptor = u_char_read ( file_input );
//
//  Read ID_LENGTH_BYTES here.
//
  if ( 0 < th->id_length )
  {
    th->image_id = new unsigned char[(th->id_length)+1];
    file_input.read ( ( char * ) th->image_id, (th->id_length) );
    th->image_id[(th->id_length)] = '\0';
  }
  else
  {
    th->image_id = NULL;
  }
  return th;
}
//****************************************************************************80

void tga_header_write ( ofstream &file_output, tga_header th )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_HEADER_WRITE writes the header of a TGA file.
// 
//  Modified:
// 
//    02 February 2006
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUTPUT, a reference to the output file.
//
//    Input, tga_header TH, a pointer to the TGA header.
//
{
  u_char_write      ( th.id_length,            file_output );
  u_char_write      ( th.color_map_type,       file_output );
  u_char_write      ( th.image_type,           file_output );
  u_short_int_write ( th.color_map_origin,     file_output );
  u_short_int_write ( th.color_map_length,     file_output );
  u_char_write      ( th.color_map_entry_size, file_output );

  u_short_int_write ( th.image_x_origin,       file_output );
  u_short_int_write ( th.image_y_origin,       file_output );
  u_short_int_write ( th.image_width,          file_output );
  u_short_int_write ( th.image_height,         file_output );
  u_char_write      ( th.image_pixel_depth,    file_output );
  u_char_write      ( th.image_descriptor,     file_output );

  if ( 0 < th.id_length )
  {
    file_output.write ( ( char * ) th.image_id, th.id_length );
  }
  return;
}
//****************************************************************************80

tga_color_map *tga_color_map_read ( ifstream &file_input, tga_header *th )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_COLOR_MAP_READ reads the color map of a TGA file.
// 
//  Modified:
// 
//    05 April 2005
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Input, tga_header *TH, a pointer to the TGA header.
//
//    Output, tga_color_map *TGA_COLOR_MAP_READ, a pointer to the TGA color map.
//
{
  tga_color_map *tc;
  int size;

  tc = new tga_color_map;

  size = th->color_map_length * th->color_map_entry_size;

  if ( 0 != th->color_map_type )
  {
    tc->color_map = new unsigned char[size];

    file_input.read ( ( char * ) tc->color_map, size );
  }
  else
  {
    tc->color_map = NULL;
  }
  return tc;
}
//****************************************************************************80

void tga_color_map_write ( ofstream &file_output, tga_header th,
  tga_color_map tc )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_COLOR_MAP_WRITE writes the color map of a TGA file.
// 
//  Modified:
// 
//    03 February 2006
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUTPUT, a reference to the output file.
//
//    Input, tga_header TH, the TGA header.
//
//    Input, tga_color_map TC, the TGA color map.
//
{
  if ( 0 != th.color_map_type )
  {
    file_output.write ( ( char * ) tc.color_map, 
      th.color_map_length * th.color_map_entry_size );
  }
  return;
}
//****************************************************************************80

tga_data *tga_data_read ( ifstream &file_input, tga_header *th, 
  tga_color_map *tc )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_DATA_READ reads the data of a TGA file.
// 
//  Modified:
// 
//    05 April 2005
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Input, tga_header *TH, a pointer to the TGA header.
//
//    Input, tga_color_map *TC, a pointer to the TGA color map.
//
//    Output, tga_data *TGA_DATA_READ, a pointer to the TGA data.
//
{
  int size;
  tga_data *td;
  unsigned char *image_data;

  if ( th->image_type == 0 )
  {
    td = new tga_data;
    td->image_data = NULL;
  }
  else if ( th->image_type == 1 )
  {
    cout << "  This program can't read IMAGE_TYPE 1 files.\n";
    td = NULL;
  } 
  else if ( th->image_type == 2 )
  {
    td = new tga_data;
    size = 3 * th->image_height * th->image_width;
    image_data = new unsigned char[size];
    file_input.read ( ( char * ) image_data, size );
    td->image_data = image_data;
  } 
  else if ( th->image_type == 3 )
  {
    cout << "  This program can't read IMAGE_TYPE 3 files.\n";
    td = NULL;
  } 
  else if ( th->image_type == 9 )
  {
    cout << "  This program can't read IMAGE_TYPE 9 files.\n";
    td = NULL;
  } 
  else if ( th->image_type == 10 )
  {
    cout << "  This program can't read IMAGE_TYPE 10 files.\n";
    td = NULL;
  } 
  else if ( th->image_type == 11 )
  {
    cout << "  This program can't read IMAGE_TYPE 11 files.\n";
    td = NULL;
  } 
  else
  {
    cout << "  Illegal value of IMAGE_TYPE = " 
      << th->image_type << ".\n";
    td = NULL;
  }

  return td;
}
//****************************************************************************80

void tga_data_write ( ofstream &file_output, tga_header th, 
  tga_color_map tc, tga_data td )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_DATA_WRITE writes the data of a TGA file.
// 
//  Modified:
// 
//    04 February 2006
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUTPUT, a reference to the output file.
//
//    Input, tga_header TH, the TGA header.
//
//    Input, tga_color_map TC, the TGA color map.
//
//    Input, tga_data TD, the TGA data.
//
{
  int size;
  unsigned char *image_data;

  if ( th.image_type == 0 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  }
  else if ( th.image_type == 1 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  } 
  else if ( th.image_type == 2 )
  {
    size = 3 * th.image_height * th.image_width;

    file_output.write ( ( char * ) td.image_data, size );
  } 
  else if ( th.image_type == 3 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  } 
  else if ( th.image_type == 9 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  } 
  else if ( th.image_type == 10 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  } 
  else if ( th.image_type == 11 )
  {
    cout << "  This program can't write IMAGE_TYPE " 
      << th.image_type << " files.\n";
  } 
  else
  {
    cout << "  Illegal value of IMAGE_TYPE = " 
      << th.image_type << ".\n";
  }

  return;
}
//****************************************************************************80

bool tga_read ( char *file_input_name, tga_header **th, tga_color_map **tc, 
  tga_data **td )

//****************************************************************************80
//
//  Purpose:
// 
//    TGA_READ reads a TGA file.
// 
//  Modified:
// 
//    31 January 2006
//
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, char *FILE_INPUT_NAME, the name of the file.
//
//    Output, bool TGA_READ, is TRUE if something or other.
//
{
  ifstream file_input;
  
  file_input.open ( file_input_name, ios::binary );

  if ( !file_input.is_open ( ) )
  {
    cout << "\n";
    cout << "TGA_READ - Fatal error!\n";
    cout << "  Could not open the file \"" << file_input_name << "\".\n";
    return false;
  }

  *th = tga_header_read ( file_input );
  *tc = tga_color_map_read ( file_input, *th );
  *td = tga_data_read ( file_input, *th, *tc );

  file_input.close ( );

  return true;
}
//****************************************************************************80

unsigned char u_char_read ( ifstream &file_input )

//****************************************************************************80
//
//  Purpose:
// 
//    U_CHAR_READ reads an unsigned char from a file.
//
//  Modified:
//
//    30 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Output, unsigned char U_CHAR_READ, the value that was read.
//
{
  char c;
  unsigned char clo;

  file_input.read ( &c, 1 );
  clo = ( unsigned char ) c;

  return clo;
}
//****************************************************************************80

void u_char_write ( unsigned char u_char_val, ofstream &file_out )

//****************************************************************************80
//
//  Purpose:
// 
//    U_CHAR_WRITE writes an unsigned char to a file.
//
//  Modified:
//
//    26 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned char U_CHAR_VAL, the value to be written.
//
//    Input, ofstream &FILE_OUT, a reference to the output file.
//
{
  file_out << u_char_val;

  return;
}
//****************************************************************************80

unsigned long int u_long_int_read ( ifstream &file_input )

//****************************************************************************80
//
//  Purpose:
// 
//    U_LONG_INT_READ reads an unsigned long int from a file.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Output, unsigned long int U_LONG_INT_READ, the value that was read.
//
{
  unsigned long int u_long_int_val;
  unsigned short int u_short_int_val_hi;
  unsigned short int u_short_int_val_lo;

  if ( tga_byte_swap )
  {
    u_short_int_val_lo = u_short_int_read ( file_input );
    u_short_int_val_hi = u_short_int_read ( file_input );
  }
  else
  {
    u_short_int_val_hi = u_short_int_read ( file_input );
    u_short_int_val_lo = u_short_int_read ( file_input );
  }

  u_long_int_val = ( u_short_int_val_hi << 16 ) | u_short_int_val_lo;

  return u_long_int_val;
}
//****************************************************************************80

void u_long_int_write ( unsigned long int u_long_int_val, 
  ofstream &file_out )

//****************************************************************************80
//
//  Purpose:
// 
//    U_LONG_INT_WRITE writes an unsigned long int to a file.
//
//  Modified:
//
//    05 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned long int *U_LONG_INT_VAL, the value to be written.
//
//    Input, ofstream &FILE_OUT, a reference to the output file.
//
{
  unsigned short int u_short_int_val_hi;
  unsigned short int u_short_int_val_lo;

  u_short_int_val_hi = ( unsigned short ) ( u_long_int_val / 65536 );
  u_short_int_val_lo = ( unsigned short ) ( u_long_int_val % 65536 );

  if ( tga_byte_swap )
  {
    u_short_int_write ( u_short_int_val_lo, file_out );
    u_short_int_write ( u_short_int_val_hi, file_out );
  }
  else
  {
    u_short_int_write ( u_short_int_val_hi, file_out );
    u_short_int_write ( u_short_int_val_lo, file_out );
  }

  return;
}
//****************************************************************************80

unsigned short int u_short_int_read ( ifstream &file_input )

//****************************************************************************80
//
//  Purpose:
// 
//    U_SHORT_INT_READ reads an unsigned short int from a file.
//
//  Modified:
//
//    30 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ifstream &FILE_INPUT, a reference to the input file.
//
//    Output, unsigned short int U_SHORT_INT_READ, the value that was read.
//
{
  char c;
  unsigned char chi;
  unsigned char clo;
  unsigned short int u_short_int_val;

  if ( tga_byte_swap )
    {
    file_input.read ( &c, 1 );
    clo = ( unsigned char ) c;
    file_input.read ( &c, 1 );
    chi = ( unsigned char ) c;
  }
  else
  {
    file_input.read ( &c, 1 );
    chi = ( unsigned char ) c;
    file_input.read ( &c, 1 );
    clo = ( unsigned char ) c;
  }

  u_short_int_val = ( chi << 8 ) | clo;

  return u_short_int_val;
}
//****************************************************************************80

void u_short_int_write ( unsigned short int u_short_int_val, 
  ofstream &file_out )

//****************************************************************************80
//
//  Purpose:
// 
//    U_SHORT_INT_WRITE writes an unsigned short int to a file.
//
//  Modified:
//
//    26 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, unsigned short int *U_SHORT_INT_VAL, the value to be written.
//
//    Input, ofstream &FILE_OUT, a reference to the output file.
//
{
  unsigned char chi;
  unsigned char clo;

  chi = ( unsigned char ) ( u_short_int_val / 256 );
  clo = ( unsigned char ) ( u_short_int_val % 256 );

  if ( tga_byte_swap )
  {
    file_out << clo << chi;
  }
  else
  {
    file_out << chi << clo;
  }

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
