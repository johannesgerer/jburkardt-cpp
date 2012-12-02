bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
void grf_data_print ( int node_num, int edge_num, int edge_pointer[], 
  int edge_data[], double xy[] );
void grf_data_read ( string input_filename, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] );
void grf_data_write ( ofstream &output_unit, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] );
void grf_example ( int node_num, int edge_num, int edge_pointer[], 
  int edge_data[], double xy[] );
void grf_example_size ( int *node_num, int *edge_num );
void grf_header_print ( int node_num, int edge_num );
void grf_header_read ( string input_filename, int *node_num, int *edge_num );
void grf_header_write ( string output_filename, ofstream &output_unit,
  int node_num, int edge_num );
void grf_write ( string output_filename, int node_num, int edge_num, 
  int edge_pointer[], int edge_data[], double xy[] );
void i4vec_copy ( int n, int a1[], int a2[] );
void r8vec_copy ( int n, double a1[], double a2[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
double s_to_r8 ( string s, int *lchar, bool *error );
int s_word_count ( string s );
void timestamp ( void );
