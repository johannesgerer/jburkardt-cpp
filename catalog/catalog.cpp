# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <string>
# include <fstream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int get_file_type ( string file_name );
void handle_c ( string file_name );
void handle_cc ( string file_name );
void handle_f90 ( string file_name );
void handle_m ( string file_name );
bool s_begin ( string s1, string s2 );
string s_last_ch ( string s, char ch );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CATALOG.
//
//  Discussion:
//
//    The INDEX program prints those lines of:
//
//      a C file that are two lines after the line "  Purpose:";
//      a C++ file that are two lines after the line "//  Purpose:";
//      a FORTRAN77 file that begin with "cc";
//      a FORTRAN90 file that begin with "!!";
//      a MATLAB file that begin with "%%".
//
//    This allows me to document a program by inserting a one line
//    description of each routine in the body of the code.  Then the INDEX
//    command applied to the file will give me a quick summary of what's there.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  int C = 1;
  int CC = 2;
  int F90 = 4;
  string file_name;
  int file_type;
  int i;
  int M = 5;
  int UNKNOWN = -1;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "CATALOG:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Extract lines in a file that mark each routine.\n";
    cout << "  My C routines have a special Purpose: line.\n";
    cout << "  My C++ routines have a special Purpose: line.\n";
    cout << "  My FORTRAN77 routines have a special cc line.\n";
    cout << "  My FORTRAN90 routines have a special !! line.\n";
    cout << "  My MATLAB routines have a special %% line.\n";
  }
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "  Please enter the name of a file to be indexed.\n";

    cin >> file_name;

    file_type = get_file_type ( file_name );

    if ( file_type == C ) 
    {
      handle_c ( file_name );
    }
    else if ( file_type == CC ) 
    {
      handle_cc ( file_name );
    }
    else if ( file_type == F90 )
    {
      handle_f90 ( file_name );
    }
    else if ( file_type == M )
    {
      handle_m ( file_name );
    }
    else if ( file_type == UNKNOWN )
    {
      cerr << "\n";
      cerr << "CATALOG - Fatal error!\n" ;
      cerr << "  The file type of \"" << file_name 
           << "\" could not be determined.\n";
    }  
  }
//
//  Otherwise, get the file(s) from the argument list. 
//
  else 
  {

    for ( i = 1 ; i < argc ; ++i ) 
    {
      file_name = argv[i];

      file_type = get_file_type ( file_name );

      if ( file_type == C ) 
      {
        handle_c ( file_name );
      }
      else if ( file_type == CC ) 
      {
        handle_cc ( file_name );
      }
      else if ( file_type == F90 )
      {
        handle_f90 ( file_name );
      }
      else if ( file_type == M )
      {
        handle_m ( file_name );
      }
      else if ( file_type == UNKNOWN )
      {
        cerr << "\n";
        cerr << "CATALOG - Fatal error!\n" ;
        cerr << "  The file type of \"" << file_name 
             << "\" could not be determined.\n";
      }
    }
  } 
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "CATALOG:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

int get_file_type ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    GET_FILE_TYPE determines the type of a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file.
//
//    Output, int FILE_TYPE, is C, CC, F90, M, or UNKNOWN.
//
{
  int C = 1;
  int CC = 2;
  string ext;
  int F90 = 4;
  int file_type;
  int M = 5;
  int UNKNOWN = -1;

  ext = s_last_ch ( file_name, '.' );
 
  if ( ext.empty ( ) ) 
  {
    file_type = UNKNOWN;
  }
  else if ( ext.compare ( ".f" ) == 0 )
  {
    file_type = F90;
  }
  else if ( ext.compare ( ".f77" ) == 0 )
  {
    file_type = F90;
  }
  else if ( ext.compare ( ".f90" ) == 0 )
  {
    file_type = F90;
  }
  else if ( ext.compare ( ".f95" ) == 0 )
  {
    file_type = F90;
  }
  else if ( ext.compare ( ".for" ) == 0 )
  {
    file_type = F90;
  }
  else if ( ext.compare ( ".c" ) == 0 )
  {
    file_type = C;
  }
  else if ( ext.compare ( ".cc" ) == 0 )
  {
    file_type = CC;
  }
  else if ( ext.compare ( ".C" ) == 0 )
  {
    file_type = CC;
  }
  else if ( ext.compare ( ".cxx" ) == 0 )
  {
    file_type = CC;
  }
  else if ( ext.compare ( ".cpp" ) == 0 )
  {
    file_type = CC;
  }
  else if ( ext.compare ( ".m" ) == 0 )
  {
    file_type = M;
  }
  else
  {
    file_type = UNKNOWN;
  }

  return file_type;
}
//****************************************************************************80

void handle_c ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_C searches a single C file for lines beginning with '  Purpose:'.
//
//  Discussion:
//
//    If the routine finds a line beginning with "  Purpose:", then it prints
//    out the line that comes two lines later, assuming this to be the one-line
//    text corresponding to that header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to index.
//
{
  ifstream file_in;
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  string line;
  int line_length;
//
//  Open the file.
//
  file_in.open ( file_name.c_str ( ) );

  if ( !file_in ) 
  {
    cerr << "\n";
    cerr << "HANDLE_C - Fatal error!\n";
    cerr << "  Cannot open \"" << file_name << "\".\n";
    return;
  }
     
  cout << "\n";
  cout << "Index of \"" << file_name << "\":\n";
  cout << "\n";

  while ( 1 ) 
  {
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) )
    {
      break;
    }

    i0 = i1;
    i1 = i2;

    if ( s_begin ( line, "  Purpose:" ) )
    {
      i2 = 1;
    }
    else
    {
      i2 = 0;
    }

    if ( i0 == 1 )
    {
      cout << line.substr ( 2 ) << "\n";
    }

  }
//
//  Close the file.
//
  file_in.close ( );

  return;
}
//****************************************************************************80

void handle_cc ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_CC searches a single C++ file for lines beginning with '//  Purpose:'.
//
//  Discussion:
//
//    If the routine finds a line beginning with "  Purpose:", then it prints
//    out the line that comes two lines later, assuming this to be the one-line
//    text corresponding to that header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to index.
//
{
  ifstream file_in;
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  string line;
  int line_length;
//
//  Open the file.
//
  file_in.open ( file_name.c_str ( ) );

  if ( !file_in ) 
  {
    cerr << "\n";
    cerr << "HANDLE_CC - Fatal error!\n";
    cerr << "  Cannot open \"" << file_name << "\".\n";
    return;
  }
     
  cout << "\n";
  cout << "Index of \"" << file_name << "\":\n";
  cout << "\n";

  while ( 1 ) 
  {  
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) )
    {
      break;
    }

    i0 = i1;
    i1 = i2;

    if ( s_begin ( line, "//  Purpose:" ) )
    {
      i2 = 1;
    }
    else
    {
      i2 = 0;
    }
    if ( i0 == 1 )
    {
      cout <<  line.substr ( 4 ) << "\n";
    }
  }
//
//  Close the file.
//
  file_in.close ( );   

  return;
}
//****************************************************************************80

void handle_f90 ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_F90 searches a single FORTRAN file for lines beginning with '!!'.
//
//  Discussion;
//
//    FORTRAN77 files might use the 'CC' or 'cc' marker instead.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to index.
//
{
  ifstream file_in;
  string line;
  int line_length;
//
//  Open the file.
//
  file_in.open ( file_name.c_str ( ) );

  if ( !file_in ) 
  {
    cerr << "\n";
    cerr << "HANDLE_F90 - Fatal error!\n";
    cerr << "  Cannot open \"" << file_name << "\".\n";
    return;
  }
     
  cout << "\n";
  cout << "  Index of \"" << file_name << "\":\n";
  cout << "\n";

  while ( 1 ) 
  {  
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) )
    {
      break;
    }

    if ( s_begin ( line, "!!" ) || 
         s_begin ( line, "cc" ) || 
         s_begin ( line, "CC" ) )
    {
      cout << " " << line.substr ( 2 ) << "\n";
    }
  }
//
//  Close the file.
//
  file_in.close ( );

  return;
}
//****************************************************************************80

void handle_m ( string file_name )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_M searches a single MATLAB file for lines beginning with '%%'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to index.
//
{
  ifstream file_in;
  string line;
  int line_length;
//
//  Open the file.
//
  file_in.open ( file_name.c_str ( ) );

  if ( !file_in ) 
  {
    cerr << "\n";
    cerr << "HANDLE_M - Fatal error!\n";
    cerr << "  Cannot open \"" << file_name << "\".\n";
    return;
  }
     
  cout << "\n";
  cout << "  Index of \"" << file_name << "\":\n";
  cout << "\n";

  while ( 1 ) 
  {  
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) )
    {
      break;
    }

    if ( s_begin ( line, "%%" ) )
    {
      cout << " " << line.substr ( 2 ) << "\n";
    }
  }
//
//  Close the file.
//
  file_in.close ( );

  return;
}
//****************************************************************************80

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
{
  int i;
  int n1;
  int n2;

  n1 = s1.length ( );
  n2 = s2.length ( );

  if ( n1 < n2 )
  {
    return false;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( s1[i] != s2[i] ) 
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

string s_last_ch ( string s, char ch )

//****************************************************************************80
//
//  Purpose:
//
//    S_LAST_CH points to the last occurrence of a character in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Input, char CH, a character.
//
//    Output, string S_LAST_CH, the substring beginning with the last occurrence
//    of the given character.
//
{
  int position;
  int s_length;
  string t;
  int t_length;

  s_length = s.length ( );
//
//  Find the last occurrence.
//
  for ( position = s_length - 1; 0 <= position; position-- )
  {
    if ( s[position] == ch )
    {
      t_length = s_length - position;
      t = s.substr ( position, t_length );
      return t;
    }
  }

  t.clear ( );

  return t;
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
//    02 October 2003
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
