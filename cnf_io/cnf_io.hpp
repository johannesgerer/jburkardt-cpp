char ch_cap ( char ch );
bool ch_eqi ( char ch1, char ch2 );
bool ch_is_space ( char c );
bool cnf_data_read ( string cnf_file_name, int v_num, int c_num, 
  int l_num, int l_c_num[], int l_val[] );
bool cnf_data_write ( int c_num, int l_num, int l_c_num[], int l_val[], 
  ofstream &output_unit );
bool cnf_evaluate ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[], 
  bool v_val[] );
bool cnf_header_read ( string cnf_file_name, int *v_num, int *c_num, 
  int *l_num );
bool cnf_header_write ( int v_num, int c_num, string output_name, 
  ofstream &output_unit );
void cnf_print ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[] );
bool cnf_write ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[], 
  string output_name );
int i4_power ( int i, int j );
void lvec_next ( int n, bool lvec[] );
string s_adjustl ( string s1 );
string s_blanks_delete ( string s );
bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
void s_word_extract_first ( string s, string &s1, string &s2 );
void timestamp ( void );

