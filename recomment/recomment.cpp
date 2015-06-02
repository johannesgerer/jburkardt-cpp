# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char ch );
void handle ( string file_in_name, string file_out_name );
bool s_eqi ( string s1, string s2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    RECOMMENT makes a copy of a C file, converting to C++ comment style.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  string file_in_name;
  string file_out_name;
  bool VERBOSE = false;

  if ( VERBOSE )
  {
    timestamp ( );
    cout << "\n";
    cout << "RECOMMENT:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Replace C comments by C++ comments.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  }
//
//  If the input file was not specified on the command line, get it now.
//
  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "RECOMMENT:\n";
    cout << "  Please enter the name of the file to be recommented.\n";

    cin >> file_in_name;
  }
  else
  {
    file_in_name = argv[1];
  }
//
//  If the output file was not specified on the command line, get it now.
//
  if ( argc < 3 )
  {
    cout << "\n";
    cout << "RECOMMENT:\n";
    cout << "  Please enter the OUTPUT file name:\n";

    cin >> file_out_name;
  }
  else
  {
    file_out_name = argv[2];
  }

  if ( s_eqi ( file_in_name, file_out_name ) )
  {
    cout << "\n";
    cout << "RECOMMENT - Fatal error!\n";
    cout << "  As a safety precaution, the input and output files\n";
    cout << "  MUST differ by more than just case.\n";
    cout << "\n";
    cout << "  Abnormal end of execution\n";
    timestamp ( );
    return 1;
  }

  handle ( file_in_name, file_out_name );
//
//  Terminate.
//
  if ( VERBOSE )
  {
    cout << "\n";
    cout << "RECOMMENT:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
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

void handle ( string file_in_name, string file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE recomments a single C file.
//
//  Discussion:
//
//    Each line of the input file is copied to the output file.
//
//    However, lines that are C comments are reformatted to the
//    C++ comment style.
//
//    This routine incorporates suggestions and coding provided on
//    28 April 2005 by Steven M Martin of JDS Uniphase, Melbourne Florida.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_IN_NAME, the name of the input file.
//
//    Input, string FILE_OUT_NAME, the name of the output file.
//
{
  char c0;
  char c1;
  char c2;
  int COMMAND = 1;
  int C_COMMENT = 2;
  int CPP_COMMENT = 3;
  ifstream file_in;
  ofstream file_out;
  int line_length;
  int mode;
  int STRING = 4;
  bool terminate;
//
//  Open the input file.
//
  file_in.open ( file_in_name.c_str ( ) );

  if ( !file_in )
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open the input file \"" << file_in_name << "\"\n";
    return;
  }
//
//  Open the output file.
//
  file_out.open ( file_out_name.c_str ( ) );

  if ( !file_out )
  {
    cout << "\n";
    cout << "HANDLE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\"\n";
    return;
  }

  file_out << "// File recommented by recomment.cpp\n";
  file_out << "// on " << __DATE__ << " at " << __TIME__ << ".\n";
  file_out << "//\n";
//
//  For each character read in, can we output this character now?
//
  line_length = 0;
  c0 = '\0';
  c1 = '\0';
  c2 = '\0';
  mode = COMMAND;
  terminate = false;
//
//  C0 has been written.
//  C1 has been read, and we're thinking about writing it.
//  C2 is the next character.
//
  while ( !terminate )
  {
    c0 = c1;
    c1 = c2;
    file_in.get ( c2 );

    if ( file_in.eof ( ) )
    {
      terminate = true;
    }
//
//  Processing "ordinary" text.
//
    if ( mode == COMMAND )
    {
      if ( !c1 )
      {
      }
      else if ( c1 == '\n' )
      {
        file_out.put ( c1 );
        line_length = 0;
      }
      else if ( c1 == '/' && c2 == '*' )
      {
        if ( 0 < line_length )
        {
          file_out.put ( '\n' );
          line_length = 0;
        }
        mode = C_COMMENT;
        file_out.put ( '/' );
        file_out.put ( '/' );
        line_length = line_length + 2;
        c0 = '/';
        c1 = '/';
        c2 = '\0';
      }
      else if ( c1 == '/' && c2 == '/' )
      {
        if ( 0 < line_length )
        {
          file_out.put ( '\n' );
          line_length = 0;
        }
        mode = CPP_COMMENT;
        file_out.put ( '/' );
        file_out.put ( '/' );
        line_length = line_length + 2;
        c0 = '/';
        c1 = '/';
        c2 = '\0';
      }
//
//  Skip over the double quote inside a single quote.
//
      else if ( c1 == '\'' && c2 == '"' )
      {
        file_out.put ( '\'' );
        file_out.put ( '"' );
        line_length = line_length + 2;
        c0 = '\'';
        c1 = '"';
        c2 = '\0';
      }
//
//  Starting a string.
//
      else if ( c1 == '"' )
      {
        file_out.put ( '"' );
        line_length = line_length + 1;
        mode = STRING;
      }
      else
      {
        file_out.put ( c1 );
        line_length = line_length + 1;
      }
    }
//
//  Processing a C comment.
//
    else if ( mode == C_COMMENT )
    {
      if ( !c1 )
      {
      }
      else if ( c1 == '\n' )
      {
        file_out.put ( c1 );
        file_out.put ( '/' );
        file_out.put ( '/' );
        c0 = '/';
        c1 = '/';
        line_length = 2;
      }
//
//    if ( 0 < line_length )
//  CHANGED TO
//    if ( 2 < line_length )
//  07 July 2014
//
      else if ( c1 == '*' && c2 == '/' )
      {
        if ( 2 < line_length )
        {
          file_out.put ( '\n' );
        }
        c1 = '\0';
        c2 = '\0';
        mode = COMMAND;
      }
      else
      {
        file_out.put ( c1 );
        line_length = line_length + 1;
      }
    }
//
//  Processing a C++ comment.
//
    else if ( mode == CPP_COMMENT )
    {
      if ( !c1 )
      {
      }
      else if ( c1 == '\n' )
      {
        file_out.put ( c1 );
        line_length = 0;
        mode = COMMAND;
      }
      else
      {
        file_out.put ( c1 );
        line_length = line_length + 1;
      }
    }
//
//  Processing a string.
//
    else if ( mode == STRING )
    {

      if ( c1 == '\\' && c2 == '"' )
      {
        file_out.put ( c1 );
        file_out.put ( c2 );
        line_length = line_length + 2;
        c2 = '\0';
      }
      else if ( c1 == '"' )
      {
        file_out.put ( c1 );
        line_length = line_length + 1;
        mode = COMMAND;
      }
      else
      {
        file_out.put ( c1 );
        line_length = line_length + 1;
      }
    }
  }
//
// Close the files.
//
  file_in.close ( );
  file_out.close ( );

  return;
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
//    04 October 2003
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
