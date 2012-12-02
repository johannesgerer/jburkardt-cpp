# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "cnf_io.hpp"

//****************************************************************************80

char ch_cap ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  return ( ch1 == ch2 );
}
//****************************************************************************80

bool ch_is_space ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_SPACE is TRUE if a character represents "white space".
//
//  Discussion:
//
//    A white space character is a space, a form feed, a newline, a carriage
//    return, a horizontal tab, or a vertical tab.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_SPACE, is TRUE if C is a whitespace character.
//
{
  bool value;

  if ( c == ' ' )
  {
    value = true;
  }
  else if ( c == '\f' )
  {
    value = true;
  }
  else if ( c == '\n' )
  {
    value = true;
  }
  else if ( c == '\r' )
  {
    value = true;
  }
  else if ( c == '\t' )
  {
    value = true;
  }
  else if ( c == '\v' )
  {
    value = true;
  }
  else
  {
    value = false;
  }
  return value;
}
//****************************************************************************80

bool cnf_data_read ( string cnf_file_name, int v_num, int c_num, 
  int l_num, int l_c_num[], int l_val[] )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_DATA_READ reads the data of a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CNF_FILE_NAME, the name of the CNF file.
//
//    Input, int V_NUM, the number of variables.
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, int L_NUM, the number of signed literals.
//
//    Output, int L_C_NUM[C_NUM], the number of signed
//    literals occuring in each clause.
//
//    Output, int L_VAL[L_NUM], a list of all the signed 
//    literals in all the clauses, ordered by clause.
//
//    Output, bool CNF_DATA_READ, is TRUE if there was an error during 
//    the read.
//
{
  int c_num2;
  bool error;
  int ierror;
  ifstream input;
  int l_c_num2;
  int l_num2;
  int l_val2;
  int length;
  string line;
  string rest;
  int v_num2;
  string word;
  
  error = false;

  input.open ( cnf_file_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "CNF_DATA_READ - Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines until you find one that is not blank and does not begin
//  with a "c".  This should be the header line.
//
  while ( 1 )
  {
    getline ( input, line );
    
    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "CNF_DATA_READ - Fatal error!\n";
      cout << "  Error3 while reading the file.\n";
      exit ( 1 );
    }

    if ( line[0] == 'c' || line[0] == 'C' )
    {
      continue;
    }

    if ( 0 < s_len_trim ( line ) )
    {
      break;
    }
  }
//
//  We expect to be reading the line "p cnf V_NUM C_NUM"
//
  if ( line[0] != 'p' && line[0] != 'P' )
  {
    cout << "\n";
    cout << "CNF_DATA_READ - Fatal error!\n";
    cout << "  First non-comment non-blank line does not start\n";
    cout << "  with 'p ' marker.\n";
    exit ( 1 );
  }

  if ( !ch_is_space ( line[1] ) )
  {
    cout << "\n";
    cout << "CNF_DATA_READ - Fatal error!\n";
    cout << "  Character after 'p' must be whitespace.\n";
    exit ( 1 );
  }
//
//  Remove the first two characters and shift left.
//
  line[0] = ' ';
  line[1] = ' ';
  line = s_adjustl ( line );
//
//  Expect the string 'CNF'
//
  if ( ch_eqi ( line[0], 'c' ) && 
       ch_eqi ( line[1], 'n' ) && 
       ch_eqi ( line[2], 'f' ) )
  {
  }
  else
  {
    cout << "\n";
    cout << "CNF_DATA_READ - Fatal error!\n";
    cout << "  First non-comment non-blank line does not start\n";
    cout << "  with 'p cnf' marker.\n";
    exit ( 1 );
  }

  if ( !ch_is_space ( line[3] ) ) 
  {
    cout << "\n";
    cout << "CNF_DATA_READ - Fatal error!\n";
    cout << "  Character after 'p cnf' must be whitespace.\n";
    exit ( 1 );
  }
//
//  Remove the first four characters and shift left.
//
  line[0] = ' ';
  line[1] = ' ';
  line[2] = ' ';
  line[3] = ' ';

  line = s_adjustl ( line );
//
//  Extract the next word, which is the number of variables.
//  You can compare this to V_NUM for an extra check.
//
  sscanf ( line.c_str ( ), "%d  %d", &v_num2, &c_num2 );
//
//  Read remaining lines, counting the literals while ignoring occurrences of '0'.
//
  l_num2 = 0;
  c_num2 = 0;
  l_c_num2 = 0;

  while ( 1 )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == 'c' || line[0] == 'C' )
    {
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    while ( 1 )
    {
      s_word_extract_first ( line, word, rest );
      line = rest;

      if ( s_len_trim ( word ) <= 0 )
      {
        break;
      }

      l_val2 = s_to_i4 ( word, &length, &error );

      if ( error )
      {
        break;
      }

      if ( l_val2 != 0 )
      {
        l_val[l_num2] = l_val2;
        l_num2 = l_num2 + 1;
        l_c_num2 = l_c_num2 + 1;
      }
      else
      {
        l_c_num[c_num2] = l_c_num2;
        c_num2 = c_num2 + 1;
        l_c_num2 = 0;
      }
    }
  }
//
//  At the end:
//
//    C_NUM2 should equal C_NUM.
//    L_NUM2 should equal L_NUM.
//
//  Close file and return.
//
  input.close ( );

  return error;
} 
//****************************************************************************80

bool cnf_data_write ( int c_num, int l_num, int l_c_num[], int l_val[], 
  ofstream &output_unit )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_DATA_WRITE writes data to a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, int L_NUM, the total number of signed literals.
//
//    Input, int L_C_NUM[C_NUM], the number of signed
//    literals occuring in each clause.
//
//    Input, int L_VAL[L_NUM], a list of all the signed 
//    literals in all the clauses, ordered by clause.
//
//    Input, ofstream &OUTPUT_UNIT, the output unit.
//
{
  int c;
  bool error;
  int i1;
  int i2;
  int l;
  int l_c;

  error = false;

  l = 0;

  for ( c = 0; c < c_num; c++ )
  {
    i1 = 1;
    i2 = 10;
    for ( l_c = 0; l_c < l_c_num[c]; l_c++ )
    {
      output_unit << " " << setw(7) << l_val[l];
      l = l + 1;

      if ( ( ( l_c + 1 ) % 10 ) == 0 )
      {
        output_unit << "\n";
      }
    }
    output_unit << " " << setw(7) << 0 << "\n";
  }

  return error;
}
//****************************************************************************80

bool cnf_evaluate ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[], 
  bool v_val[] )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_EVALUATE evaluates a formula in CNF form.
//
//  Discussion:
//
//    The formula is in conjunctive normal form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int V_NUM, the number of variables.
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, int L_NUM, the total number of signed literals.
//
//    Input, int L_C_NUM[C_NUM], the number of signed
//    literals occuring in each clause.
//
//    Input, int L_VAL[L_NUM], a list of all the signed 
//    literals in all the clauses, ordered by clause.
//
//    Input, bool V_VAL[V_NUM], the values assigned to the variables.
//
//    Output, bool CNF_EVALUATE, the value of the CNF formula for the
//    given variable values.
//
{
  int c;
  bool c_val;
  bool f_val;
  int l;
  int l_c;
  bool s_val;
  int v_index;

  f_val = true;

  l = 0;

  for ( c = 0; c < c_num; c++ )
  {
//
//  The clause is false unless some signed literal is true.
//
    c_val = false;
    for ( l_c = 0; l_c < l_c_num[c]; l_c++ )
    {
      s_val = ( 0 < l_val[l] );
      v_index = abs ( l_val[l] );
      l = l + 1;
//
//  The signed literal is true if the sign "equals" the value.
//  Note that we CAN'T exit the loop because we need to run out the 
//  L index!
//
      if ( v_val[v_index-1] == s_val )
      {
        c_val = true;
      }
    }
//
//  The formula is false if any clause is false.
//
    if ( !c_val )
    {
      f_val = false;
      break;
    }
  }

  return f_val;
}
//****************************************************************************80

bool cnf_header_read ( string cnf_file_name, int *v_num, int *c_num, 
  int *l_num )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_HEADER_READ reads the header of a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CNF_FILE_NAME, the name of the CNF file.
//
//    Output, int *V_NUM, the number of variables.
//
//    Output, int *C_NUM, the number of clauses.
//
//    Output, int *L_NUM, the number of signed literals.
//
//    Output, bool CNF_HEADER_READ, is TRUE if there was an error during 
//    the read.
//
{
  ifstream input;
  bool error;
  int ierror;
  int l_val;
  int length;
  string line;
  string rest;
  string word;
  
  error = false;

  input.open ( cnf_file_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Could not open file.\n";
    exit ( 1 );
  }
//
//  Read lines until you find one that is not blank and does not begin
//  with a "c".  This should be the header line.
//
  while ( 1 )
  {
    getline ( input, line );
    
    if ( input.eof ( ) )
    {
      cout << "\n";
      cout << "CNF_HEADER_READ - Fatal error!\n";
      cout << "  Error3 while reading the file.\n";
      exit ( 1 );
    }

    if ( line[0] == 'c' || line[0] == 'C' )
    {
      continue;
    }

    if ( 0 < s_len_trim ( line ) )
    {
      break;
    }
  }
//
//  We expect to be reading the line "p cnf V_NUM C_NUM"
//
  if ( line[0] != 'p' && line[0] != 'P' )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  First non-comment non-blank line does not start\n";
    cout << "  with 'p ' marker.\n";
    exit ( 1 );
  }

  if ( !ch_is_space ( line[1] ) )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Character after 'p' must be whitespace.\n";
    exit ( 1 );
  }
//
//  Remove the first two characters and shift left.
//
  line[0] = ' ';
  line[1] = ' ';
  line = s_adjustl ( line );
//
//  Expect the string 'CNF'
//
  if ( ch_eqi ( line[0], 'c' ) && 
       ch_eqi ( line[1], 'n' ) && 
       ch_eqi ( line[2], 'f' ) )
  {
  }
  else
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  First non-comment non-blank line does not start\n";
    cout << "  with 'p cnf' marker.\n";
    exit ( 1 );
  }

  if ( !ch_is_space ( line[3] ) ) 
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Character after 'p cnf' must be whitespace.\n";
    exit ( 1 );
  }
//
//  Remove the first four characters and shift left.
//
  line[0] = ' ';
  line[1] = ' ';
  line[2] = ' ';
  line[3] = ' ';

  line = s_adjustl ( line );
//
//  Extract the next word, which is the number of variables.
//
  s_word_extract_first ( line, word, rest );
  line = rest;

  if ( s_len_trim ( word ) <= 0 )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Unexpected End of input.\n";
    exit ( 1 );
  }

  *v_num = s_to_i4 ( word, &length, &error );

  if ( error )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Unexpected End of input.\n";
    exit ( 1 );
  }
//
//  Extract the next word, which is the number of clauses.
//
  s_word_extract_first ( line, word, rest );
  line = rest;

  if ( s_len_trim ( word ) <= 0 )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Unexpected End of input.\n";
    exit ( 1 );
  }

  *c_num = s_to_i4 ( word, &length, &error );

  if ( error )
  {
    cout << "\n";
    cout << "CNF_HEADER_READ - Fatal error!\n";
    cout << "  Unexpected End of input.\n";
    exit ( 1 );
  }
//
//  Read remaining lines, counting the literals while ignoring occurrences of '0'.
//
  *l_num = 0;

  while ( 1 )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == 'c' || line[0] == 'C' )
    {
      continue;
    }

    if ( s_len_trim ( line ) < 0 )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    while ( 1 )
    {
      s_word_extract_first ( line, word, rest );
      line = rest;

      if ( s_len_trim ( word ) <= 0 )
      {
        break;
      }

      l_val = s_to_i4 ( word, &length, &error );

      if ( error )
      {
        break;
      }

      if ( l_val != 0 )
      {
        *l_num = *l_num + 1;
      }
    }
  }
//
//  Close file and return.
//
  input.close ( );

  return error;
} 
//****************************************************************************80

bool cnf_header_write ( int v_num, int c_num, string output_name, 
  ofstream &output_unit )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_HEADER_WRITE writes the header for a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int V_NUM, the number of variables.
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, string OUTPUT_NAME, the name of the output file.
//
//    Input, ofstream &OUTPUT_UNIT, the output unit.
//
{
  bool error;

  error = false;

  output_unit << "c " << output_name << "\n";
  output_unit << "c\n";
  output_unit << "p cnf " << v_num << " " << c_num << "\n";
  
  return error;
}
//****************************************************************************80

void cnf_print ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[] )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_PRINT prints CNF information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int V_NUM, the number of variables.
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, int L_NUM, the total number of signed literals.
//
//    Input, int L_C_NUM[C_NUM], the number of signed
//    literals occuring in each clause.
//
//    Input, int L_VAL[L_NUM], a list of all the signed 
//    literals in all the clauses, ordered by clause.
//
{
  int c;
  int l;
  int l_c;

  cout << "\n";
  cout << "CNF data printout:\n";
  cout << "\n";
  cout << "  The number of variables       V_NUM  = " << v_num << "\n";
  cout << "  The number of clauses         C_NUM  = " << c_num << "\n";
  cout << "  The number of signed literals L_NUM  = " << l_num << "\n";
  l = 0;
  for ( c = 0; c < c_num; c++ )
  {
    cout << "\n";
    cout << "  Clause " << c
         << " includes " << l_c_num[c] << " signed literals:\n";
    for ( l_c = 0; l_c < l_c_num[c]; l_c++ )
    {
      cout << setw(4) << l_val[l] << "\n";
      l = l + 1;
    }
  }
  return;
}
//****************************************************************************80

bool cnf_write ( int v_num, int c_num, int l_num, int l_c_num[], int l_val[], 
  string output_name )

//****************************************************************************80
//
//  Purpose:
//
//    CNF_WRITE writes the header and data of a CNF file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 June 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int V_NUM, the number of variables.
//
//    Input, int C_NUM, the number of clauses.
//
//    Input, int L_NUM, the total number of signed literals.
//
//    Input, int L_C_NUM[C_NUM], the number of signed
//    literals occuring in each clause.
//
//    Input, int L_VAL[L_NUM], a list of all the signed 
//    literals in all the clauses, ordered by clause.
//
//    Input, string OUTPUT_NAME, the name of the output file.
//
{
  bool error;
  ofstream output_unit;

  error = false;
//
//  Open the output file.
//
  output_unit.open ( output_name.c_str ( ) );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "CNF_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << output_name << "\".\n";
    error = true;
    return error;
  }
//
//  Write the header.
//
  error = cnf_header_write ( v_num, c_num, output_name, output_unit );

  if ( error )
  {
    cout << "\n";
    cout << "CNF_WRITE - Fatal error!\n";
    cout << "  Cannot write the header for the output file \"" << output_name << "\".\n";
    return error;
  }
//
//  Write the data.
//
  error = cnf_data_write ( c_num, l_num, l_c_num, l_val, output_unit );

  if ( error )
  {
    cout << "\n";
    cout << "CNF_WRITE - Fatal error!\n";
    cout << "  Cannot write the data for the output file \"" << output_name << "\".\n";
    return error;
  }
//
//  Close the file.
//
  output_unit.close ( );

  return error;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void lvec_next ( int n, bool lvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    LVEC_NEXT generates the next logical vector.
//
//  Discussion:
//
//    In the following discussion, we will let '0' stand for FALSE and
//    '1' for TRUE.
//
//    The vectors have the order
//
//      (0,0,...,0),
//      (0,0,...,1), 
//      ...
//      (1,1,...,1)
//
//    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
//    we allow wrap around.
//
//  Example:
//
//    N = 3
//
//    Input      Output
//    -----      ------
//    0 0 0  =>  0 0 1
//    0 0 1  =>  0 1 0
//    0 1 0  =>  0 1 1
//    0 1 1  =>  1 0 0
//    1 0 0  =>  1 0 1
//    1 0 1  =>  1 1 0
//    1 1 0  =>  1 1 1
//    1 1 1  =>  0 0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input/output, bool LVEC[N], on output, the successor to the
//    input vector.  
//
{
  int i;

  for ( i = n - 1; 0 <= i; i-- )
  {
    if ( !lvec[i] )
    {
      lvec[i] = true;
      return;
    }
    lvec[i] = false;
  }
  return;
}
//****************************************************************************80

string s_adjustl ( string s1 )

//****************************************************************************80
//
//  Purpose:
//
//    S_ADJUSTL flushes a string left.
//
//  Discussion:
//
//    Both blanks and tabs are treated as "white space".
//
//    This routine is similar to the FORTRAN90 ADJUSTL routine.
//
//  Example:
//
//    Input             Output
//
//    '     Hello'      'Hello'
//    ' Hi there!  '    'Hi there!'
//    'Fred  '          'Fred'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, the string to be adjusted.
//
//    Output, string S_ADJUSTL, the adjusted string.
//
{
  int i;
  int s2_length;
  string s2;
  int nonb;
  char TAB = 9;

  s2 = s1;
//
//  Check the length of the string to the last nonblank.
//  If nonpositive, return.
//
  s2_length = s2.length ( );

  if ( s2_length <= 0 )
  {
    return s2;
  }
//
//  Find NONB, the location of the first nonblank, nontab.
//
  nonb = 0;

  for ( i = 0; i < s2_length; i++ )
  {
    if ( s1[i] != ' ' && s1[i] != TAB )
    {
      nonb = i;
      break;
    }
  }

  if ( 0 < nonb )
  {
    for ( i = 0; i < s2_length - nonb; i++ )
    {
      s2[i] = s1[i+nonb];
    }
    for ( i = s2_length - nonb; i < s2_length; i++ )
    {
      s2[i] = ' ';
    }

  }
  return s2;
}
//****************************************************************************80

string s_blanks_delete ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_BLANKS_DELETE replaces consecutive blanks by one blank.
//
//  Discussion:
//
//    The remaining characters are left justified and right padded with blanks.
//    TAB characters are converted to spaces.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be transformed.
//
//    Output, string S_BLANKS_DELETE, the transformed string.
//
{
  bool blank;
  char c;
  int get;
  int put;
  int s_length;
  char *s2;
  string s3;

  s_length = s.length ( );
  s2 = new char[s_length+1];
  s2[s_length] = '\0';

  blank = true;
  put = 0;

  for ( get = 0; get < s_length; get++ )
  {
    if ( s[get] != ' ' )
    {
      s2[put] = s[get];
      put = put + 1;
      blank = false;
    }
    else if ( !blank )
    {
      s2[put] = s[get];
      put = put + 1;
      blank = true;
    }
    else
    {
    }
  }
//
//  Suppress a final blank that is not the only character.
//
  if ( 1 < put )
  {
    if ( s2[put-1] == ' ' )
    {
      put = put - 1;
    }
  }
  s2[put] = '\0';
  s3 = string ( s2 );

  return s3;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
    } 
  }

  return true;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( string s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for ( ; ; ) 
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

void s_word_extract_first ( string s, string &s1, string &s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_EXTRACT_FIRST extracts the first word from a string.
//
//  Discussion:
//
//    A "word" is a string of characters terminated by a blank or
//    the end of the string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string.
//
//    Output, string &S1, the first word (initial blanks removed).
//
//    Output, string &S2, the remainder of the string, after removing
//    the first word (initial blanks removed).
//
{
  int i;
  int mode;
  int s_len;

  s_len = s.length ( );
  s1 = "";
  s2 = "";
  mode = 1;

  for ( i = 0; i < s_len; i++ )
  {
    if ( mode == 1 )
    {
      if ( s[i] != ' ' )
      {
         mode = 2;
      }
    }
    else if ( mode == 2 )
    {
      if ( s[i] == ' ' )
      {
        mode = 3;
      }
    }
    else if ( mode == 3 )
    {
      if ( s[i] != ' ' )
      {
        mode = 4;
      }
    }
    if ( mode == 2 )
    {
      s1 = s1 + s[i];
    }
    else if ( mode == 4 )
    {
      s2 = s2 + s[i];
    }
  }

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
