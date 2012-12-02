char ch_cap ( char ch );

bool pgmb_check_data ( int xsize, int ysize, unsigned char maxg, unsigned char *g );
bool pgmb_example ( int xsize, int ysize, unsigned char *g );

bool pgmb_read ( string file_in_name, int &xsize, int &ysize, 
  unsigned char &maxg, unsigned char **g );
bool pgmb_read_data ( ifstream &file_in, int xsize, int ysize, 
  unsigned char *g );
bool pgmb_read_header ( ifstream &file_in, int &xsize, int &ysize, 
  unsigned char &maxg );
bool pgmb_read_test ( string file_in_name );

bool pgmb_write ( string file_out_name, int xsize, int ysize, unsigned char *g );
bool pgmb_write_data ( ofstream &file_out, int xsize, int ysize, 
  unsigned char *g );
bool pgmb_write_header ( ofstream &file_out, int xsize, int ysize, 
  unsigned char maxg );
bool pgmb_write_test ( string file_out_name );

bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
