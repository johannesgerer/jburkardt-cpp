bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
void i4vec_copy ( int n, int a1[], int a2[] );
void r8vec_copy ( int n, double a1[], double a2[] );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( void );
void xy_data_print ( int point_num, double xy[] );
void xy_data_read ( string input_filename, int point_num, double xy[] );
void xy_data_write ( ofstream &file_out, int point_num, double xy[] );
void xy_example ( int point_num, double xy[] );
void xy_header_print ( int point_num );
void xy_header_read ( string input_filename, int *point_num );
void xy_header_write ( string file_out_name, ofstream &file_out, int point_num );
void xy_read ( string file_in_name, int *point_num, double *xy[] );
void xy_read_test ( string file_in_name );
void xy_write ( string file_out_name, int point_num, double xy[] );
void xy_write_test ( string file_out_name );
void xyf_data_print ( int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyf_data_read ( string input_filename, int face_num, int face_data_num,
  int face_pointer[], int face_data[] );
void xyf_data_write ( ofstream &output_unit, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyf_example ( int point_num, int face_num, int face_data_num, double xy[],
  int face_pointer[], int face_data[] );
void xyf_example_size ( int *point_num, int *face_num, int *face_data_num );
void xyf_header_print ( int point_num, int face_num, int face_data_num );
void xyf_header_read ( string input_filename, int *face_num, int *face_data_num );
void xyf_header_write ( string output_filename, ofstream &output_unit, 
  int point_num, int face_num, int face_data_num );
void xyf_write ( string output_filename, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyl_data_print ( int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
void xyl_data_read ( string input_filename, int line_num, int line_data_num,
  int line_pointer[], int line_data[] );
void xyl_data_write ( ofstream &file_out, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
void xyl_example ( int point_num, int line_num, int line_data_num, double xy[],
  int line_pointer[], int line_data[] );
void xyl_example_size ( int *point_num, int *line_num, int *line_data_num );
void xyl_header_print ( int point_num, int line_num, int line_data_num );
void xyl_header_read ( string input_filename, int *line_num, int *line_data_num );
void xyl_header_write ( string file_out_name, ofstream &file_out, int point_num,
  int line_num, int line_data_num );
void xyl_write ( string file_out_name, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
