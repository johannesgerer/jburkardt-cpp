struct tga_header
{
  unsigned char color_map_entry_size;
  unsigned short int color_map_length;
  unsigned short int color_map_origin;
  unsigned char color_map_type;
  unsigned char id_length;
  unsigned char image_descriptor;
  unsigned short int image_height;
  unsigned char *image_id;
  unsigned char image_pixel_depth;
  unsigned char image_type;
  unsigned short int image_width;
  unsigned short int image_x_origin;
  unsigned short int image_y_origin;
};

struct tga_color_map
{
  unsigned char *color_map;
};

struct tga_data
{
  unsigned char *image_data;
};

bool tga_byte_swap_get ( void );
void tga_byte_swap_set ( bool value );

tga_color_map *tga_color_map_read ( ifstream &file_in, tga_header *th );
void tga_color_map_write ( ofstream &file_output, tga_header th,
  tga_color_map tc );

tga_data *tga_data_read ( ifstream &file_in, tga_header *tg, 
  tga_color_map *tc );
void tga_data_write ( ofstream &file_output, tga_header th, 
  tga_color_map tc, tga_data td );

void tga_header_print ( tga_header *th );
tga_header *tga_header_read ( ifstream &file_in );
void tga_header_write ( ofstream &file_output, tga_header th );

bool tga_read ( char *file_name, tga_header **th, tga_color_map **tc, 
  tga_data **td );

void timestamp ( void );
unsigned char u_char_read ( ifstream &file_in );
void u_char_write ( unsigned char u_char_val, ofstream &file_out );
unsigned long int u_long_int_read ( ifstream &file_in );
void u_long_int_write ( unsigned long int u_long_int_val, ofstream &file_out );
unsigned short int u_short_int_read ( ifstream &file_in );
void u_short_int_write ( unsigned short int u_short_int_val, ofstream &file_out );
