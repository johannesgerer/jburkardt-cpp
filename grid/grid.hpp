char digit_to_ch ( int i );
unsigned long get_seed ( );
void grid_generate ( int m, int n, int center, int *seed, double r[] );
int grid_side ( int m, int n );
int i4_log_10 ( int i );
char *i4_to_s ( int i );
void ksub_random2 ( int n, int k, int *seed, int a[] );
double r8_uniform_01 ( int *seed );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );
void tuple_next_fast ( int m, int n, int rank, int x[] );

