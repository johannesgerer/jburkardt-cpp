# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <string>

using namespace std;

# include "linpack_d.H"
# include "blas1_d.H"

int main ( int argc, char *argv[] );
void basis_write ( string file_out_name, int m, int n, double s, double u[],
  bool comment  );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
bool ch_is_digit ( char c );
int ch_to_digit ( char c );
char digit_inc ( char c );
char digit_to_ch ( int i );
int file_column_count ( string input_filename );
bool file_exist ( string file_name );
void file_name_inc_nowrap ( string *file_name );
int file_row_count ( string input_filename );
int i4_huge ( );
double r8_epsilon ( );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_print ( int m, int n, double a[], string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void singular_vectors ( int m, int n, int basis_num, double a[], double sval[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SVD_BASIS.
//
//  Discussion:
//
//    SVD_BASIS forms a basis from the SVD of a set of data vectors.
//
//    This program uses the singular value decomposition (SVD) to analyze
//    a set of data, and extract a number of dominant modes.
//
//    This program is intended as an intermediate application, in
//    the following situation:
//
//    A) a "high fidelity" or "high resolution" PDE solver is used
//       to determine many (say N = 500) solutions of a discretized
//       PDE at various times, or parameter values.  Each solution
//       may be regarded as an M vector.  Typically, each solution
//       involves an M by M linear system, greatly reduced in
//       complexity because of bandedness or sparsity.
//
//    B) This program is applied to extract L dominant modes from
//       the N solutions.  This is done using the singular value
//       decomposition of the M by N matrix, each of whose columns
//       is one of the original solution vectors.
//
//    C) a "reduced order model" program may then attempt to solve
//       a discretized version of the PDE, using the L dominant
//       modes as basis vectors.  Typically, this means that a dense
//       L by L linear system will be involved.
//
//    Thus, the program might read in 500 files, and write out
//    5 or 10 files of the corresponding size and "shape", representing
//    the dominant solution modes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
{
# define DATA_FILE_BASE_MAX 20

  string basis_file;
  int basis_num;
  bool comment;
  char comment_char;
  int comp_num;
  string data_file;
  int data_file_base_num;
  string data_file_base[DATA_FILE_BASE_MAX];
  int data_file_num;
  int dim_num;
  bool error;
  string file_name;
  int i;
  int ii;
  int j;
  int k;
  int l;
  int node_num;
  double *point;
  int point_num;
  double *sval;
  double *table;

  timestamp ( );

  cout << "\n";
  cout << "SVD_BASIS:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Given a PDE for which:\n";
  cout << "    C is the number of components of the solution\n";
  cout << "      at any single point,\n";
  cout << "    P is the number of points where a solution is given,\n";
  cout << "    N is the number of solution vectors,\n";
  cout << "    L is the number of modes to be extracted.\n";
  cout << "\n";
  cout << "  Then we let M = C*P be the abstract spatial dimension.\n";
  cout << "\n";
  cout << "  Set up A, the M by N matrix of solution vectors,\n";
  cout << "\n";
  cout << "  Get A = U * S * V', the singular value decomposition.\n";
  cout << "\n";
  cout << "  The first L columns of U are our modes.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
//
//  What is the basis size?
//
  cout << "  How many basis vectors (L) are to be extracted?\n";
  cin >> basis_num;

  cout << "\n";
  cout << "  L = " << basis_num << "\n";
//
//  Gather one or more "base" file names.
//
  data_file_base_num = 0;

  for ( ; ; )
  {
    if ( DATA_FILE_BASE_MAX <= data_file_base_num )
    {
      cout << "\n";
      cout << "  No more base file names can be entered.\n";
      return 1;
    }
//
//  Get the next base file name.
//
    cout << "\n";
    cout << "  You specify a consecutive sequence of file names\n";
    cout << "  by giving the first \"base\" file name.\n";
    cout << "\n";
    cout << "  If there are no more sequences to enter,\n";
    cout << "  just hit RETURN.\n";
    cout << "\n";
    cout << "  Enter a new base file name, or $ if done:\n";
//
//  CIN won't allow you to enter a blank line!
//  I just don't have the energy today to replace CIN by GETLINE....
//
    cin >> file_name;

    if ( file_name == "$" )
    {
      file_name = " ";
    }

    if ( s_len_trim ( file_name ) <= 0 )
    {
      cout << "\n";
      cout << "  RETURN was entered.\n";
      cout << "  Presumably, there are no more file sequences.\n";
      break;
    }

    data_file_base[data_file_base_num] = file_name;
    data_file_base_num = data_file_base_num + 1;

    cout << "\n";
    cout << data_file_base_num << ":  \"" << file_name << "\"\n";
//
//  For the very first base file, get the data sizes.
//
    if ( data_file_base_num == 1 )
    {
      r8mat_header_read ( file_name, &comp_num, &node_num );

      dim_num = comp_num * node_num;

      cout << "\n";
      cout << "  According to the first base file,\n";
      cout << "  The number of solution components C =   " << comp_num << "\n";
      cout << "  The number of solution points P =       " << node_num << "\n";
      cout << "  The \"size\" of each solution M = (C*P) = " << dim_num << "\n";
//
//  Idiocy check.  L must be less than or equal to M.
//
      if ( dim_num < basis_num )
      {
        cout << "\n";
        cout << "SVD_BASIS - Fatal error!\n";
        cout << "\n";
        cout << "  M < L.\n";
        cout << "\n";
        cout << "  That is, the number of modes requested (L) is greater\n";
        cout << "  than the spatial dimension (M).\n";
        cout << "  Technically, the program could pad out the answer\n";
        cout << "  with L-M zero vectors, but instead, we will stop\n";
        cout << "  assuming you made an error, or a misapprehension.\n";
        cout << "\n";
        cout << "SVD_BASIS:\n";
        cout << "  Abnormal end of execution.\n";
        cout << "\n";
        timestamp ( );
        return 1;
      }

    }

  }
//
//  Count all the data files.
//
  data_file_num = 0;

  for ( i = 0; i < data_file_base_num; i++ )
  {
    data_file = data_file_base[i];

    for ( ; ; )
    {
      if ( !file_exist ( data_file ) )
      {
        break;
      }
      data_file_num = data_file_num + 1;

      file_name_inc_nowrap ( &data_file );
    }
  }

  if ( data_file_num == 0 )
  {
    cout << "\n";
    cout << "SVD_BASIS - Fatal error!\n";
    cout << "  There do not seem to be any solution files;\n";
    cout << "  that is, files whose names are \"incremented\"\n";
    cout << "  versions of the first file name.\n";
    cout << "\n";
    cout << "  The first file we looked for was \"" << data_file << "\"\n";
    cout << "\n";
    cout << "SVD_BASIS:\n";
    cout << "  Abnormal end of execution.\n";
    cout << "\n";
    timestamp ( );
    return 1;
  }

  cout << "\n";
  cout << "  The number of data files N = " << data_file_num << "\n";
//
//  Set up an array to hold all the data.
//
  point_num = data_file_num;

  cout << "\n";
  cout << "  The data is stored in an M by N matrix A.\n";
  cout << "\n";
  cout << "  The \"spatial\" dimension M is " << dim_num << "\n";
  cout << "  The number of data points N is " << point_num << "\n";
//
//  Allocate space for the POINT array.
//
  point = new double[dim_num*point_num];
//
//  Read the data.
//
  l = 0;

  for ( ii = 0; ii < data_file_base_num; ii++ )
  {
    data_file = data_file_base[ii];

    for ( ; ; )
    {
      if ( !file_exist ( data_file ) )
      {
        break;
      }

      table = r8mat_data_read ( data_file, comp_num, node_num );

      k = 0;
      for ( j = 0; j < node_num; j++ )
      {
        for ( i = 0; i < comp_num; i++ )
        {
          point[k+l*dim_num] = table[i+j*comp_num];
          k = k + 1;
        }
      }

      l = l + 1;
      file_name_inc_nowrap ( &data_file );

      delete [] table;
    }
  }
  cout << "\n";
  cout << "  The data has been read into the matrix A.\n";
//
//----------------------------------------------------------------------------
//
//  Compute the SVD of A.
//
//----------------------------------------------------------------------------
//
  sval = new double[basis_num];

  singular_vectors ( dim_num, point_num, basis_num, point, sval );
//
//----------------------------------------------------------------------------
//
//  Write the first L left singular vectors (columns of U) to files.
//
//----------------------------------------------------------------------------
//
  cout << "\n";
  cout << "SVD_BASIS:\n";
  cout << "  Ready to write the left singular vectors to files.\n";
  cout << "\n";
  cout << "  Do you want comments in the header of the file?\n";
  cout << "  (These begin with the \"#\" character.) (Y/N)\n";
  cout << "\n";
  cout << "  Enter \"Y\" or \"N\":\n";

  cin >> comment_char;

  if ( comment_char == 'Y' || comment_char == 'y' )
  {
    comment = true;
  }
  else
  {
    comment = false;
  }

  basis_file = "svd_000.txt";

  for ( j = 0; j < basis_num; j++ )
  {
    file_name_inc_nowrap ( &basis_file );

    if ( j + 1 == 1 )
    {
      cout << "\n";
      cout << "  Writing first file " << basis_file  << "\n";
    }

    if ( j + 1 == basis_num )
    {
      cout << "  Writing last file  " << basis_file << "\n";
    }
    basis_write ( basis_file, comp_num, node_num, sval[j], 
      point+0+j*dim_num, comment );
  }
//
//  Free memory.
//
  delete [] point;
  delete [] sval;
//
//  Terminate.
//
  cout << "\n";
  cout << "SVD_BASIS:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
# undef DATA_FILE_BASE_MAX
}
//****************************************************************************80

void basis_write ( string file_out_name, int m, int n, double s, double u[],
  bool comment  )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_WRITE writes a basis vector to a file.
//
//  Discussion:
//
//    The initial lines of the file are comments, which begin with a
//    "#" character.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_OUT_NAME, the name of the file to write.
//
//    Input, int M, the number of data components.
//
//    Input, int N, the number of data items.
//
//    Input, double S, the associated singular value.
//
//    Input, double U[M*N], the data values.
//
//    Input, bool COMMENT, is TRUE if comments are to be included.
//
{
  ofstream file_out;
  int i;
  int j;
  int mhi;
  int mlo;

  file_out.open ( file_out_name.c_str ( ) );

  if ( !file_out )
  {
    cout << "\n";
    cout << "BASIS_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  if ( comment )
  {
    file_out << "#  " << file_out_name << "\n";
    file_out << "#  created by routine BASIS_WRITE.C" << "\n";
    file_out << "#  part of SVD_BASIS.C." << "\n";
    file_out << "#\n";
    file_out << "#  Number of components M =  " << setw(12) << m << "\n";
    file_out << "#  Number of items N =       " << setw(12) << n << "\n";
    file_out << "#  Singular value S =        " << setw(14) << s << "\n";
    file_out << "#  EPSILON (unit roundoff) = " << r8_epsilon ( ) << "\n";
    file_out << "#\n";
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      file_out << setw(10) << u[i+j*m] << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  return;
}
//****************************************************************************80*

char ch_cap ( char c )

//****************************************************************************80*
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
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 ) 
  {
    c = c - 32;
  }   

  return c;
}
//****************************************************************************80*

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80*
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
//    Input, char C1, C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 ) 
  {
    c1 = c1 - 32;
  } 
  if ( 97 <= c2 && c2 <= 122 ) 
  {
    c2 = c2 - 32;
  }     

  return ( c1 == c2 );
}
//****************************************************************************80

bool ch_is_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_IS_DIGIT returns TRUE if a character is a decimal digit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to be analyzed.
//
//    Output, bool CH_IS_DIGIT, is TRUE if C is a digit.
//
{
  if ( '0' <= c && c <= '9' )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
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
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

char digit_inc ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_INC increments a decimal digit.
//
//  Example:
//
//    Input  Output
//    -----  ------
//    '0'    '1'
//    '1'    '2'
//    ...
//    '8'    '9'
//    '9'    '0'
//    'A'    'A'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, a digit to be incremented.
//
//    Output, char DIGIT_INC, the incremented digit.
//
{
  if ( '0' <= c && c <= '8' )
  {
    return ( c + 1 );
  }
  else if ( c == '9' ) 
  {
    return '0';
  }
  else
  {
    return c;
  }
}
//****************************************************************************80

char digit_to_ch ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    DIGIT_TO_CH returns the base 10 digit character corresponding to a digit.
//
//  Example:
//
//     I     C
//   -----  ---
//     0    '0'
//     1    '1'
//   ...    ...
//     9    '9'  
//    10    '*'
//   -83    '*'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the digit, which should be between 0 and 9.
//
//    Output, char DIGIT_TO_CH, the appropriate character '0' through '9' or '*'.
//
{
  char c;

  if ( 0 <= i && i <= 9 )
  {
    c = '0' + i;
  }
  else
  {
    c = '*';
  }

  return c;
}
//****************************************************************************80

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words,
//    separated by spaces.  There may also be some blank lines, and some
//    comment lines, which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
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
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  string text;
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    exit ( 1 );
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) <= 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input >> text;

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
      {
        continue;
      }
      got_one = true;
      break;
    }
  }

  input.close ( );

  if ( !got_one )
  {
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

bool file_exist ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_EXIST reports whether a file exists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, bool FILE_EXIST, is TRUE if the file exists.
//
{
  ifstream file;
  bool value;

  file.open ( filename.c_str ( ), ios::in );

  if ( !file )
  {
    value = false;
  }
  else
  {
    value = true;
  }
  return value;
}
//****************************************************************************80

void file_name_inc_nowrap ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_NAME_INC_NOWRAP increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the (nonempty) name contains no digits, or all the digits are
//    9, then the empty string is returned.
//
//    If the empty string is input, the routine stops.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a8to99.txt"     "a9to00.txt"
//      "a9to99.txt"     " "
//      "cat.txt"        " "
//      " "              STOP!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int carry;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILE_NAME_INC_NOWRAP - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;
  carry = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;
      carry = 0;

      if ( c == '9' )
      {
        carry = 1;
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  Unsatisfied carry.  The input digits were all 9.  Return blank.
//
  if ( carry == 1 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

int file_row_count ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  string line;
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4 value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" integer.
//
{
  return 2147483647;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the roundoff unit for R8 arithmetic.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the property 
//    that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but 
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the double precision round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R8MAT_DATA_READ[M*N], the data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  double *table;
  double *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80

void r8mat_header_read ( string input_filename, int *m, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HEADER_READ reads the header from an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    exit ( 1 );
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
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
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
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
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
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
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
}
//****************************************************************************80

void singular_vectors ( int m, int n, int basis_num, double a[], double sval[] )

//****************************************************************************80
//
//  Purpose:
//
//    SINGULAR_VECTORS computes the desired singular values.
//
//  Discussion:
//
//    The LINPACK SVD routine DSVDC is used to compute the singular
//    value decomposition:
//
//      A = U * S * V'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of data points.
//
//    Input, int BASIS_NUM, the number of basis vectors to be extracted.
//
//    Input/output, double A[M*N]; on input, the matrix whose 
//    singular values are to be computed.  On output, A(M,1:BASIS_NUM)
//    contains the first BASIS_NUM left singular vectors.
//
//    Output, double SVAL[BASIS_NUM], the first BASIS_NUM
//    singular values.
//
{
  double *e;
  int i;
  int info;
  int lda;
  int ldu;
  int ldv;
  int job;
  double *s;
  double *u;
  double *v;
  double *work;

  cout << "\n";
  cout << "SINGULAR_VECTORS\n";
  cout << "  For an MxN matrix A in general storage,\n";
  cout << "  The LINPACK routine DSVDC computes the\n";
  cout << "  singular value decomposition:\n";
  cout << "\n";
  cout << "    A = U * S * V'\n";
  cout << "\n";
//
//  Compute the eigenvalues and eigenvectors.
//
  s = new double[ i4_min ( m + 1, n ) ];
  e = new double[n];
  lda = m;
  u = a;
  ldu = m;
  v = NULL;
  ldv = n;
  work = new double[m];
  job = 20;

  info = dsvdc ( a, lda, m, n, s, e, u, ldu, v, ldv, work, job );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "SINGULAR_VECTORS - Warning:\n";
    cout << "  DSVDC returned nonzero INFO = " << info << "\n";;
    return;
  }

  for ( i = 0; i < basis_num; i++ )
  {
    sval[i] = s[i];
  }

  cout << "\n";
  cout << "  The leading singular values:\n";
  cout << "\n";

  for ( i = 0; i < basis_num; i++ )
  {
    cout << "  "
         << setw(4)  << i+1     << "  "
         << setw(16) << sval[i] << "\n";
  }

  delete [] e;
  delete [] s;
  delete [] work;

  return;
}
//****************************************************************************80

void timestamp ( )

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
