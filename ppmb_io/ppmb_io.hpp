char ch_cap ( char ch );

int i4_max ( int i1, int i2 );

bool ppmb_check_data ( int xsize, int ysize, int maxrgb, unsigned char *r,
  unsigned char *g, unsigned char *b );

bool ppmb_example ( int xsize, int ysize, unsigned char *r, 
  unsigned char *g, unsigned char *b );

bool ppmb_read ( string file_in_name, int &xsize, int &ysize, int &maxrgb,
  unsigned char **r, unsigned char **g, unsigned char **b );
bool ppmb_read_data ( ifstream &file_in, int xsize, int ysize, unsigned char *r,
  unsigned char *g, unsigned char *b );
bool ppmb_read_header ( ifstream &file_in, int &xsize, int &ysize, int &maxrgb );
bool ppmb_read_test ( string file_in_name );

bool ppmb_write ( string file_out_name, int xsize, int ysize, unsigned char *r,
  unsigned char *g, unsigned char *b );
bool ppmb_write_data ( ofstream &file_out, int xsize, int ysize, unsigned char *r,
  unsigned char *g, unsigned char *b );
bool ppmb_write_header ( ofstream &file_out, int xsize, int ysize, int maxrgb );
bool ppmb_write_test ( string file_out_name );

bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
