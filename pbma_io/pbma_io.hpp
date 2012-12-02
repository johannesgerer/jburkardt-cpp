void pbma_check_data ( int xsize, int ysize, int *barray );
void pbma_example ( int xsize, int ysize, int *barray );

void pbma_read ( string file_in_name, int &xsize, int &ysize, int **barrary );
void pbma_read_data ( ifstream &file_in, int xsize, int ysize, int *barray );
void pbma_read_header ( ifstream &file_in, int &xsize, int &ysize );
void pbma_read_test ( string file_in_name );

void pbma_write ( string file_out_name, int xsize, int ysize, int *barray );
void pbma_write_data ( ofstream &file_out, int xsize, int ysize, int *barray );
void pbma_write_header ( ofstream &file_out, string file_out_name, int xsize,
       int ysize );
void pbma_write_test ( string file_out_name );

int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
