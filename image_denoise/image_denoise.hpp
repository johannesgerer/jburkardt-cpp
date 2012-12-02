int *gray_median_news ( int m, int n, int gray[] );
int i4vec_frac ( int n, int a[], int k );
int i4vec_median ( int n, int a[] );
void pgma_read_data ( ifstream &file_in, int xsize, int ysize, int *g );
void pgma_read_header ( ifstream &file_in, int *xsize, int *ysize, int *maxg );
void pgma_write ( string output_name, int xsize, int ysize, int *g );
void pgma_write_data ( ofstream &output, int xsize, int ysize, int *g );
void pgma_write_header ( ofstream &output, string output_name, int xsize, 
  int ysize, int maxg );
void timestamp ( );
