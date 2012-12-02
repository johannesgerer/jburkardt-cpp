char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
int ch_to_digit ( char ch );
double *cluster_energy_compute ( int dim_num, int point_num, int cluster_num, 
  double point[], int cluster[], double cluster_center[] );
double *cluster_initialize_1 ( int dim_num, int point_num, int cluster_num, 
  double point[] );
double *cluster_initialize_2 ( int dim_num, int point_num, int cluster_num, 
  double point[], int *seed );
double *cluster_initialize_3 ( int dim_num, int point_num, int cluster_num, 
  double point[], int *seed );
double *cluster_initialize_4 ( int dim_num, int point_num, int cluster_num, 
  double point[], int *seed );
double *cluster_initialize_5 ( int dim_num, int point_num, int cluster_num, 
  double point[], int *seed );
void cluster_print_summary ( int point_num, int cluster_num, 
  int cluster_population[], double cluster_energy[], double cluster_variance[] );
double *cluster_variance_compute ( int dim_num, int point_num, int cluster_num, 
  double point[], int cluster[], double cluster_center[] );
int file_column_count ( string filename );
int file_row_count ( string input_filename );
void hmeans_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], int cluster[], double cluster_center[], 
  int cluster_population[], double cluster_energy[] );
void hmeans_02 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], int cluster[], double cluster_center[], 
  int cluster_population[], double cluster_energy[], int *seed );
void hmeans_w_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], double weight[], int cluster[], 
  double cluster_center[], int cluster_population[], double cluster_energy[] );
void hmeans_w_02 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], double weight[], int cluster[], 
  double cluster_center[], int cluster_population[], double cluster_energy[], 
  int *seed );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_uniform ( int a, int b, int *seed );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void i4vec_negone ( int n, int a[] );
int *i4vec_negone_new ( int n );
int i4vec_sum ( int n, int a[] );
void i4vec_zero ( int n, int a[] );
int *i4vec_zero_new ( int n );
void kmeans_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], int cluster[], double cluster_center[], 
  int cluster_population[], double cluster_energy[] );
void kmeans_02 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], int cluster[], double cluster_center[], 
  int cluster_population[], double cluster_energy[] );
void kmeans_02_optra ( int dim_num, int point_num, int cluster_num, 
  double point[], double cluster_center[], int cluster[], int cluster2[], 
  int cluster_population[], double an1[], double an2[], int ncp[], 
  double d[], int itran[], int live[], int &indx );
void kmeans_02_qtran ( int dim_num, int point_num, int cluster_num, 
  double point[], double cluster_center[], int cluster[], int cluster2[], 
  int cluster_population[], double an1[], double an2[], int ncp[], double d[], 
  int itran[], int &indx );
void kmeans_03 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], int cluster[], double cluster_center[], 
  int cluster_population[], double cluster_energy[] );
void kmeans_w_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], double weight[], int cluster[], 
  double cluster_center[], int cluster_population[], double cluster_energy[] );
void kmeans_w_03 ( int dim_num, int point_num, int cluster_num, int it_max, 
  int &it_num, double point[], double weight[], int cluster[], 
  double cluster_center[], int cluster_population[], double cluster_energy[] );
int r4_nint ( float x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_uniform_01 ( int *seed );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
bool r8vec_all_nonpositive ( int n, double a[] );
bool r8vec_any_negative ( int n, double a[] );
double r8vec_i4vec_dot_product ( int n, double r8vec[], int i4vec[] );
int r8vec_min_index ( int n, double a[] );
double r8vec_sum ( int n, double a[] );
void r8vec_uniform_01 ( int n, int *seed, double r[] );
double *r8mat_uniform_01_new ( int m, int n, int *seed );
void r8vec_zero ( int n, double a[] );
double *r8vec_zero_new ( int n );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );
