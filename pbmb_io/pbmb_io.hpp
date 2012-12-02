char ch_cap ( char ch );

bool pbmb_check_data ( int xsize, int ysize, int *barray );
bool pbmb_example ( int xsize, int ysize, int *barray );

bool pbmb_read ( string file_in_name, int &xsize, int &ysize, int **barray );
bool pbmb_read_data ( ifstream &file_in, int xsize, int ysize, int *barray );
bool pbmb_read_header ( ifstream &file_in, int &xsize, int &ysize );
bool pbmb_read_test ( string file_in_name );

bool pbmb_write ( string file_out_name, int xsize, int ysize, int *barray );
bool pbmb_write_data ( ofstream &file_out, int xsize, int ysize, int *barray );

bool pbmb_write_header ( ofstream &file_out, int xsize, int ysize );
bool pbmb_write_test ( string file_out_name );

bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
