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
void xyz_data_print ( int point_num, double xyz[] );
void xyz_data_read ( string input_filename, int point_num, double xyz[] );
void xyz_data_write ( ofstream &output_unit, int point_num, double xyz[] );
void xyz_example ( int point_num, double xyz[] );
int xyz_example_size ( );
void xyz_header_print ( int point_num );
void xyz_header_read ( string input_filename, int *point_num );
void xyz_header_write ( string output_filename, ofstream &output_unit, 
  int point_num );
void xyz_read ( string input_filename, int *point_num, double *xyz[] );
void xyz_read_test ( string input_filename );
void xyz_write ( string output_filename, int point_num, double xyz[] );
void xyz_write_test ( string output_filename );
void xyzf_data_print ( int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyzf_data_read ( string input_filename, int face_num, int face_data_num,
  int face_pointer[], int face_data[] );
void xyzf_data_write ( ofstream &output_unit, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyzf_example ( int point_num, int face_num, int face_data_num, 
  double xyz[], int face_pointer[], int face_data[] );
void xyzf_example_size ( int *point_num, int *face_num, int *face_data_num );
void xyzf_header_print ( int point_num, int face_num, int face_data_num );
void xyzf_header_read ( string input_filename, int *face_num, 
  int *face_data_num );
void xyzf_header_write ( string output_filename, ofstream &output_unit, 
  int point_num, int face_num, int face_data_num );
void xyzf_write ( string output_filename, int point_num, int face_num,
  int face_data_num, int face_pointer[], int face_data[] );
void xyzl_data_print ( int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
void xyzl_data_read ( string input_filename, int line_num, int line_data_num,
  int line_pointer[], int line_data[] );
void xyzl_data_write ( ofstream &output_unit, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
void xyzl_example ( int point_num, int line_num, int line_data_num, 
  double xyz[], int line_pointer[], int line_data[] );
void xyzl_example_size ( int *point_num, int *line_num, int *line_data_num );
void xyzl_header_print ( int point_num, int line_num, int line_data_num );
void xyzl_header_read ( string input_filename, int *line_num, 
  int *line_data_num );
void xyzl_header_write ( string output_filename, ofstream &output_unit, 
  int point_num, int line_num, int line_data_num );
void xyzl_write ( string output_filename, int point_num, int line_num,
  int line_data_num, int line_pointer[], int line_data[] );
