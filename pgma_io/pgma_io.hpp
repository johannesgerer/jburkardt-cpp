int i4_min ( int i1, int i2 );

void pgma_check_data ( int xsize, int ysize, int maxg, int *garray );
void pgma_example ( int xsize, int ysize, int *garray );

void pgma_read ( string file_in_name, int &xsize, int &ysize, int &maxg,
       int **garrary );
void pgma_read_data ( ifstream &file_in, int xsize, int ysize, int *garray );
void pgma_read_header ( ifstream &file_in, int &xsize, int &ysize, int &maxg );
void pgma_read_test ( string file_in_name );

void pgma_write ( string file_out_name, int xsize, int ysize, int *garray );
void pgma_write_data ( ofstream &file_out, int xsize, int ysize, int *garray );
void pgma_write_header ( ofstream &file_out, string file_out_name, int xsize, 
       int ysize, int maxg );
void pgma_write_test ( string file_out_name );

int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
