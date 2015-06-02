int main ( );
void chain_init ( int chain_num, double fit[], int gen_num, int par_num, 
  double z[] );
void chain_init_print ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_read_filename, double z[] );
void chain_outliers ( int chain_num, int gen_index, int gen_num, int par_num,
  double fit[], double z[] );
void chain_write ( string chain_filename, int chain_num, double fit[], 
  int gen_num, int par_num, double z[] );
void cr_dis_update ( int chain_index, int chain_num, double cr_dis[], 
  int cr_index, int cr_num, int cr_ups[], int gen_index, int gen_num, 
  int par_num, double z[] );
int cr_index_choose ( int cr_num, double cr_prob[] );
void cr_init ( double cr[], double cr_dis[], int cr_num, double cr_prob[], 
  int cr_ups[] );
void cr_prob_update ( double cr_dis[], int cr_num, double cr_prob[], 
  int cr_ups[] );
double *diff_compute ( int chain_num, int gen_index, int gen_num, 
  int jump_dim[], 
  int jump_num, int pair_num, int par_num, int r[], double z[] );
void dream_algm ( int chain_num, int cr_num, double fit[], int gen_num, 
  double gr[], int &gr_conv, int &gr_count, int gr_num, double gr_threshold,
  double jumprate_table[], int jumpstep, double limits[], int pair_num, 
  int par_num, int printstep, double z[] );
void filename_inc ( string *filename );
void gr_compute ( int chain_num, int gen_index, int gen_num, double gr[], 
  int &gr_conv, int &gr_count, int gr_num, double gr_threshold, int par_num, 
  double z[] );
void gr_init ( double gr[], int &gr_conv, int &gr_count, int gr_num, 
  int par_num );
void gr_write ( double gr[], string gr_filename, int gr_num, int par_num, 
  int printstep );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4mat_print ( int m, int n, int a[], string title );
void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title );
int *i4vec_zero_new ( int n );
void input_print ( string chain_filename, int chain_num, int cr_num, 
  string gr_filename, double gr_threshold, int jumpstep, double limits[],
  int gen_num, int pair_num, int par_num, int printstep, 
  string restart_read_filename, string restart_write_filename );
void jumprate_choose ( double cr[], int cr_index, int cr_num, int gen_index,
  int jump_dim[], int &jump_num, double &jumprate, double jumprate_table[],
  int jumpstep, int par_num );
double *jumprate_table_init ( int pair_num, int par_num );
void jumprate_table_print ( double jumprate_table[], int pair_num, int par_num );
int r8_round_i4 ( double x );
double *r8block_zero_new ( int l, int m, int n );
double *r8mat_zero_new ( int m, int n );
double *r8vec_copy_new ( int n, double a1[] );
void r8vec_heap_d ( int n, double a[] );
double *r8vec_mnor_sample ( double parm[] );
void r8vec_sort_heap_a ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
void r8vec_transpose_print ( int n, double a[], string title );
double *r8vec_zero_new ( int n );
void restart_read ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_read_filename, double z[] );
void restart_write ( int chain_num, double fit[], int gen_num, int par_num, 
  string restart_filename, double z[] );
double *sample_candidate ( int chain_index, int chain_num, double cr[], 
  int cr_index, int cr_num, int gen_index, int gen_num, double jumprate_table[], 
  int jumpstep, double limits[], int pair_num, int par_num, double z[] );
void sample_limits ( double limits[], int par_num, double zp[] );
double *std_compute ( int chain_num, int gen_index, int gen_num, int par_num, 
  double z[] );

