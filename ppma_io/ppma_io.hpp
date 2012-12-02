char ch_cap ( char ch );
int i4_max ( int i1, int i2 );

bool ppma_check_data ( int xsize, int ysize, int maxrgb, int *rarray,
       int *garray, int *barray );
bool ppma_example ( int xsize, int ysize, int *rarray, int *garray, int *barray );

bool ppma_read ( string file_in_name, int &xsize, int &ysize, int &maxrgb,
       int **rarrary, int **garray, int **barray );
bool ppma_read_data ( ifstream &file_in, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
bool ppma_read_header ( ifstream &file_in, int &xsize, int &ysize, int &maxrgb );
bool ppma_read_test ( string file_in_name );

bool ppma_write ( string file_out_name, int xsize, int ysize, int *rarray, 
      int *garray, int *barray );
bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *rarray,
       int *garray, int *barray );
bool ppma_write_header ( ofstream &file_out, string file_out_name, int xsize, 
       int ysize, int maxrgb );
bool ppma_write_test ( string file_out_name );

bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
