int i4_huge ( );
int *i4mat_histogram ( int m, int n, int a[], int histo_num );
int i4mat_max ( int m, int n, int a[] );
int *news ( int m, int n, int a[] );
void pbma_write ( string file_out_name, int xsize, int ysize, int *b );
void pbma_write_data ( ofstream &file_out, int xsize, int ysize, int *b );
void pbma_write_header ( ofstream &file_out, string file_out_name, int xsize, 
  int ysize );
void pgma_read_data ( ifstream &file_in, int xsize, int ysize, int *g );
void pgma_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxg );
int s_len_trim ( string s );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( );
