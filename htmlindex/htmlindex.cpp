# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
char ch_cap ( char c );
int ch_index_first ( string s, char c );
int get_file_type ( string file_name );
void handle_c ( string file_name, string s1, string s2 );
void handle_cc ( string file_name, string s1, string s2 );
void handle_f77 ( string file_name, string s1, string s2 );
void handle_f90 ( string file_name, string s1, string s2 );
void handle_m ( string file_name, string s1, string s2 );
int i4_min ( int i1, int i2 );
bool s_begin ( string s1, string s2 );
string s_cap ( string s );
string s_last_ch ( string s, char ch );
void timestamp ( );
// 
//  GLOBAL DATA
//

char date_string[30];

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HTMLINDEX.
//
//  Discussion:
//
//    HTMLINDEX makes a skeleton HTML page from a marked up source code file.
//
//    HTMLINDEX writes a skeleton HTML page to standard output, including
//    stubs for the software name, description, copying option, an
//    unnumbered list of the double exclamation lines (typically reporting
//    the routines contained in the program), and some clean up lines
//    at the end.   
//
//    In my FORTRAN files, each function and subroutine includes a single
//    line that lists its name and purpose.  This line always has the form:
//
//      !! NAME purpose
//
//    and HTMLINDEX makes an HTML page suitable for describing the file,
//    including a list of functions and subroutines compiled from these
//    lines.
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
//  Usage:
//
//    htmlindex program.f90 > program.html
//
{
  int C = 1;
  int CC = 2;
  int F77 = 3;
  int F90 = 4;
  string file_name;
  int file_type;
  string head;
  string head_cap;
  int i;
  char* ipoint;
  int j;
  int M = 5;
  time_t seconds;
  struct tm  *time_struct;
  int UNKNOWN = -1;
  bool VERBOSE = false;

  if ( VERBOSE ) 
  {
    timestamp ( );
    cout << "\n";
    cout << "HTMLINDEX:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Produce a skeleton HTML page for a marked up file.\n";
  }
//
//  Retrieve the date.
//  (This was surprisingly unpleasant!)
//
  time ( &seconds );
  time_struct = localtime ( &seconds );
  ipoint = asctime ( time_struct );
  strcpy ( date_string, ipoint );
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "HTLMINDEX:\n";
    cout << "  Please enter the name of a file to be indexed.\n";

    cin >> file_name;

    j = ch_index_first ( file_name, '.' );
    head = file_name.substr ( 0, j );
    head_cap = s_cap ( head );

    file_type = get_file_type ( file_name );

    if ( file_type == C ) 
    {
      handle_c ( file_name, head, head_cap );
    }
    else if ( file_type == CC ) 
    {
      handle_cc ( file_name, head, head_cap );
    }
    else if ( file_type == F77 )
    {
      handle_f77 ( file_name, head, head_cap );
    }
    else if ( file_type == F90 )
    {
      handle_f90 ( file_name, head, head_cap );
    }
    else if ( file_type == M )
    {
      handle_m ( file_name, head, head_cap );
    }
    else if ( file_type == UNKNOWN )
    {
      cerr << "\n";
      cerr << "HTMLINDEXINDEX - Error!\n" ;
      cerr << "  The file type of \"" << file_name 
           << "\" could not be determined.\n";
    }
  }
//
//  Otherwise, open each target file, split it and close it. 
//
  else 
  {
    for ( i = 1; i < argc; i++ ) 
    {
      file_name = argv[i];
      j = ch_index_first ( file_name, '.' );
      head = file_name.substr ( 0, j );
      head_cap = s_cap ( head );

      file_type = get_file_type ( file_name );

      if ( file_type == C ) 
      {
        handle_c ( file_name, head, head_cap );
      }
      else if ( file_type == CC ) 
      {
        handle_cc ( file_name, head, head_cap );
      }
      else if ( file_type == F77 )
      {
        handle_f77 ( file_name, head, head_cap );
      }
      else if ( file_type == F90 )
      {
        handle_f90 ( file_name, head, head_cap );
      }
      else if ( file_type == M )
      {
        handle_m ( file_name, head, head_cap );
      }
      else if ( file_type == UNKNOWN )
      {
        cerr << "\n";
        cerr << "HTMLINDEX - Error!\n" ;
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
    cout << "HTMLINDEX:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( ) ;
  }

  return 0;
}
//****************************************************************************80

char ch_cap ( char c )

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
//****************************************************************************80

int ch_index_first ( string s, char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_INDEX_FIRST finds the first occurrence of a character in a string.
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
//    Input, string S, a string to be searched.
//
//    Input, char C, the character to be searched for.
//
//    Output, int CH_INDEX_FIRST, the index of the first occurrence 
//    of the character, or -1 if it does not occur.
//
{
  int i;
  int nchar;

  nchar = s.length ( );

  for ( i = 0; i < nchar; i++ ) 
  {
    if ( s[i] == c ) 
    {
      return i;
    }
  }

  return -1;
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
//    12 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file.
//
//    Output, int FILE_TYPE, is C, CC, F, or UNKNOWN.
//
{
  int C = 1;
  int CC = 2;
  string ext;
  int F77 = 3;
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
    file_type = F77;
  }
  else if ( ext.compare ( ".f77" ) == 0 )
  {
    file_type = F77;
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
    file_type = F77;
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

void handle_c ( string file_name, string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_C processes a single C file.
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
//    Input, string FILE_NAME, the name of the file to be processed.
//
//    Input, string S1, the name of the file, minus extension.
//
//    Input, string S2, the same as S1, but capitalized.
//
{
  ifstream file_in;
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  string line;
  int line_length;
  int inc;
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
//
//  Write the header.
//
  cout << "<html>\n";
  cout << "\n";
  cout << "  <head>\n";
  cout << "    <title>\n";
  cout << "      ";
  cout << s2;
  cout << " - title goes here.\n";
  cout << "    </title>\n";
  cout << "  </head>\n";
  cout << "\n";
  cout << "  <body bgcolor=\"#EEEEEE\" link=\"#CC0000\"";
  cout << " alink=\"#FF3300\" vlink=\"#000055\">\n";
  cout << "\n";
  cout << "    <h1 align = \"center\">\n";
  cout << "      ";
  cout << s2;
  cout << " <br> heading goes here.\n";
  cout << "    </h1>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b>\n";
  cout << "      is a C (program/library) which\n";
  cout << "      (description goes here).\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Usage:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Licensing:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      The computer code and data files made available on this\n";
  cout << "      web page are distributed under\n";
  cout << "      <a href = \"../../txt/gnu_lgpl.txt\">the GNU LGPL license.</a>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Languages:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b> is available in\n";
  cout << "      <a href = \"../../c_src/" << s1 << "/" << s1 << ".html\">a C version</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Related Data and Programs:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Reference:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ol>\n";
  cout << "        <li>\n";
  cout << "          \n";
  cout << "        </li>\n";
  cout << "      </ol>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Source Code:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".c\">";
  cout << s1;
  cout << ".c</a>, the source code.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".h\">";
  cout << s1;
  cout << ".h</a>, the include file.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".sh\">";
  cout << s1;
  cout << ".sh</a>,\n";
  cout << "          BASH commands to compile the source code.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Examples and Tests:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.c\">";
  cout << s1;
  cout << "_prb.c</a>\n";
  cout << "          a sample calling program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.sh\">";
  cout << s1;
  cout << "_prb.sh</a>,\n";
  cout << "          BASH commands to compile and run the sample program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb_output.txt\">";
  cout << s1;
  cout << "_prb_output.txt</a>,\n";
  cout << "          the output file.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
//
//  Start the HTML Paragraph, and "Unnumbered List".
//
  cout << "    <h3 align = \"center\">\n";
  cout << "      List of Routines:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
//
//  Get another line of text.
//
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
//
//  If we find the marker, then the text we want is two lines later.
//
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
//
//  Begin the HTML List Item, and go into Bold face.
//
      cout << "        <li>\n";
      cout << "          <b>";
//
//  Print out all the characters, until you see a space.
//
      inc = 4;
      while ( line[inc] != ' ' && line[inc] && inc < line_length ) 
      {
        cout << line[inc];
        inc = inc + 1;
      }
//
//  Terminate Bold face, then print the rest of the characters,
//  and end the List Item.
//
      cout << "</b>";
      while ( inc < line_length )
      {
        cout << line[inc];
        inc = inc + 1;
      }
      cout << "\n";
      cout << "        </li>\n";
    }
  }
//
//  End the HTML "Unnumbered List" and Paragraph.
//
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      You can go up one level to <a href = \"../c_src.html\">\n";
  cout << "      the C source codes</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <i>\n";
  cout << "      Last revised on " << date_string;
  cout << "    </i>\n";
  cout << "\n";
  cout << "    <!-- John Burkardt -->\n";
  cout << "\n";
  cout << "  </body>\n";
  cout << "\n";
  cout << "  <!-- Initial HTML skeleton created by HTMLINDEX. -->\n";
  cout << "\n";
  cout << "</html>\n";
//
//  Close the file.
//
  file_in.close ( );

  return;    
}
//****************************************************************************80

void handle_cc ( string file_name, string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_CC processes a single C++ file.
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
//    Input, string FILE_NAME, the name of the file to be processed.
//
//    Input, string S1, the name of the file, minus extension.
//
//    Input, string S2, the same as S1, but capitalized.
//
{
  ifstream file_in;
  int i0 = 0;
  int i1 = 0;
  int i2 = 0;
  string line;
  int line_length;
  int inc;
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
//
//  Write the header.
//
  cout << "<html>\n";
  cout << "\n";
  cout << "  <head>\n";
  cout << "    <title>\n";
  cout << "      ";
  cout << s2;
  cout << " - title goes here.\n";
  cout << "    </title>\n";
  cout << "  </head>\n";
  cout << "\n";
  cout << "  <body bgcolor=\"#EEEEEE\" link=\"#CC0000\"";
  cout << " alink=\"#FF3300\" vlink=\"#000055\">\n";
  cout << "\n";
  cout << "    <h1 align = \"center\">\n";
  cout << "      ";
  cout << s2;
  cout << " <br> heading goes here.\n";
  cout << "    </h1>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b>\n";
  cout << "      is a C++ (program/library) which\n";
  cout << "      (description goes here).\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Usage:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Licensing:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      The computer code and data files made available on this\n";
  cout << "      web page are distributed under\n";
  cout << "      <a href = \"../../txt/gnu_lgpl.txt\">the GNU LGPL license.</a>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Languages:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b> is available in\n";
  cout << "      <a href = \"../../cpp_src/" << s1 << "/" << s1 << ".html\">a C++ version</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Related Data and Programs:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Reference:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ol>\n";
  cout << "        <li>\n";
  cout << "          \n";
  cout << "        </li>\n";
  cout << "      </ol>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Source Code:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".cpp\">";
  cout << s1;
  cout << ".cpp</a>, the source code.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".hpp\">";
  cout << s1;
  cout << ".hpp</a>, the include file.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".sh\">";
  cout << s1;
  cout << ".sh</a>,\n";
  cout << "          BASH commands to compile the source code.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Examples and Tests:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.cpp\">";
  cout << s1;
  cout << "_prb.cpp</a>,\n";
  cout << "          a sample calling program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.sh\">";
  cout << s1;
  cout << "_prb.sh</a>,\n";
  cout << "          BASH commands to compile and run the sample program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb_output.txt\">";
  cout << s1;
  cout << "_prb_output.txt</a>,\n";
  cout << "          the output file.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
//
//  Start the HTML Paragraph, and "Unnumbered List".
//
  cout << "    <h3 align = \"center\">\n";
  cout << "      List of Routines:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
//
//  Get another line of text.
//
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
//
//  If we find the marker, then the text we want is two lines later.
//
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
//
//  Begin the HTML List Item, and go into Bold face.
//
      cout << "        <li>\n";
      cout << "          <b>";
//
//  Print out all the characters, until you see a space.
//
      inc = 6;
      while ( line[inc] != ' ' && line[inc] && inc < line_length ) 
      {
        cout << line[inc];
        inc = inc + 1;
      }
//
//  Terminate Bold face, then print the rest of the characters,
//  and end the List Item.
//
      cout << "</b>";
      while ( inc < line_length )
      {
        cout << line[inc];
        inc = inc + 1;
      }
      cout << "\n";
      cout << "        </li>\n";
    }
  }
//
//  End the HTML "Unnumbered List" and Paragraph.
//
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      You can go up one level to <a href = \"../cpp_src.html\">\n";
  cout << "      the C++ source codes</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <i>\n";
  cout << "      Last revised on " << date_string;
  cout << "    </i>\n";
  cout << "\n";
  cout << "    <!-- John Burkardt -->\n";
  cout << "\n";
  cout << "  </body>\n";
  cout << "\n";
  cout << "  <!-- Initial HTML skeleton created by HTMLINDEX. -->\n";
  cout << "\n";
  cout << "</html>\n";
//
//  Close the file.
//
  file_in.close ( );

  return;    
}
//****************************************************************************80

void handle_f77 ( string file_name, string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_F77 processes a single FORTRAN77 file.
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
//    Input, string FILE_NAME, the name of the file to be processed.
//
//    Input, string S1, the name of the file, minus extension.
//
//    Input, string S2, the same as S1, but capitalized.
//
{
  ifstream file_in;
  string line;
  int line_length;
  int inc;
//
//  Open the file.
//
  file_in.open ( file_name.c_str ( ) );

  if ( !file_in ) 
  {
    cerr << "\n";
    cerr << "HANDLE_F77 - Fatal error!\n";
    cerr << "  Cannot open \"" << file_name << "\".\n";
    return;
  }
//
//  Write the header.
//
  cout << "<html>\n";
  cout << "\n";
  cout << "  <head>\n";
  cout << "    <title>\n";
  cout << "      " << s2 << " - title goes here.\n";
  cout << "    </title>\n";
  cout << "  </head>\n";
  cout << "\n";
  cout << "  <body bgcolor=\"#EEEEEE\" link=\"#CC0000\"";
  cout << " alink=\"#FF3300\" vlink=\"#000055\">\n";
  cout << "\n";
  cout << "    <h1 align = \"center\">\n";
  cout << "      ";
  cout << s2;
  cout << " <br> heading goes here.\n";
  cout << "    </h1>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b>\n";
  cout << "      is a FORTRAN77 (program/library) which\n";
  cout << "      (description goes here).\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Usage:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Licensing:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      The computer code and data files made available on this\n";
  cout << "      web page are distributed under\n";
  cout << "      <a href = \"../../txt/gnu_lgpl.txt\">the GNU LGPL license.</a>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Languages:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b> is available in\n";
  cout << "      <a href = \"../../f77_src/" << s1 << "/" << s1 << ".html\">a FORTRAN77 version</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Related Data and Programs:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Reference:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ol>\n";
  cout << "        <li>\n";
  cout << "          \n";
  cout << "        </li>\n";
  cout << "      </ol>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Source Code:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".f\">";
  cout << s1;
  cout << ".f</a>, the source code.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".sh\">";
  cout << s1;
  cout << ".sh</a>,\n";
  cout << "          BASH commands to compile the source code.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Examples and Tests:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.f\">";
  cout << s1;
  cout << "_prb.f</a>,\n";
  cout << "          a sample calling program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.sh\">";
  cout << s1;
  cout << "_prb.sh</a>,\n";
  cout << "          BASH commands to compile and run the sample program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb_output.txt\">";
  cout << s1;
  cout << "_prb_output.txt</a>,\n";
  cout << "          the output file.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
//
//  Start the HTML Paragraph, and "Unnumbered List".
//
  cout << "    <h3 align = \"center\">\n";
  cout << "      List of Routines:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
//
//  Get another line of text.
//
  while ( 1 ) 
  {  

    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) ) 
    {
      break;
    }
//
//  If the text begins with "!!" then we want to print it,
//  but skipping the first three characters.
//
    if ( s_begin ( line, "!!" ) || 
         s_begin ( line, "cc" ) || 
         s_begin ( line, "CC" ) )
    {
//
//  Begin the HTML List Item, and go into Bold face.
//
      cout << "        <li>\n";
      cout << "          <b>";
//
//  Print out all the characters, until you see a space.
//
      inc = 3;
      while ( line[inc] != ' ' && line[inc] && inc < line_length ) 
      {
        cout << line[inc];
        inc = inc + 1;
      }
//
//  Terminate Bold face, then print the rest of the characters,
//  and end the List Item.
//
      cout << "</b>";
      while ( inc < line_length )
      {
        cout << line[inc];
        inc = inc + 1;
      }
      cout << "\n";
      cout << "        </li>\n";
    }
  }
//
//  End the HTML "Unnumbered List" and Paragraph.
//
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      You can go up one level to <a href = \"../f77_src.html\">\n";
  cout << "      the FORTRAN77 source codes</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <i>\n";
  cout << "      Last revised on " << date_string;
  cout << "    </i>\n";
  cout << "\n";
  cout << "    <!-- John Burkardt -->\n";
  cout << "\n";
  cout << "  </body>\n";
  cout << "\n";
  cout << "  <!-- Initial HTML skeleton created by HTMLINDEX. -->\n";
  cout << "\n";
  cout << "</html>\n";
//
//  Close the file.
//
  file_in.close ( );

  return;    
}
//****************************************************************************80

void handle_f90 ( string file_name, string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_F90 processes a single FORTRAN90 file.
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
//    Input, string FILE_NAME, the name of the file to be processed.
//
//    Input, string S1, the name of the file, minus extension.
//
//    Input, string S2, the same as S1, but capitalized.
//
{
  ifstream file_in;
  string line;
  int line_length;
  int inc;
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
//
//  Write the header.
//
  cout << "<html>\n";
  cout << "\n";
  cout << "  <head>\n";
  cout << "    <title>\n";
  cout << "      " << s2 << " - title goes here.\n";
  cout << "    </title>\n";
  cout << "  </head>\n";
  cout << "\n";
  cout << "  <body bgcolor=\"#EEEEEE\" link=\"#CC0000\"";
  cout << " alink=\"#FF3300\" vlink=\"#000055\">\n";
  cout << "\n";
  cout << "    <h1 align = \"center\">\n";
  cout << "      ";
  cout << s2;
  cout << " <br> heading goes here.\n";
  cout << "    </h1>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b>\n";
  cout << "      is a FORTRAN90 (program/library) which\n";
  cout << "      (description goes here).\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Usage:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Licensing:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      The computer code and data files made available on this\n";
  cout << "      web page are distributed under\n";
  cout << "      <a href = \"../../txt/gnu_lgpl.txt\">the GNU LGPL license.</a>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Languages:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b> is available in\n";
  cout << "      <a href = \"../../f_src/" << s1 << "/" << s1 << ".html\">a FORTRAN90 version</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Related Data and Programs:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Reference:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ol>\n";
  cout << "        <li>\n";
  cout << "          \n";
  cout << "        </li>\n";
  cout << "      </ol>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Source Code:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".f90\">";
  cout << s1;
  cout << ".f90</a>, the source code.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".sh\">";
  cout << s1;
  cout << ".sh</a>,\n";
  cout << "          BASH commands to compile the source code.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Examples and Tests:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.f90\">";
  cout << s1;
  cout << "_prb.f90</a>,\n";
  cout << "          a sample calling program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb.sh\">";
  cout << s1;
  cout << "_prb.sh</a>,\n";
  cout << "          BASH commands to compile and run the sample program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_prb_output.txt\">";
  cout << s1;
  cout << "_prb_output.txt</a>,\n";
  cout << "          the output file.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
//
//  Start the HTML Paragraph, and "Unnumbered List".
//
  cout << "    <h3 align = \"center\">\n";
  cout << "      List of Routines:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
//
//  Get another line of text.
//
  while ( 1 ) 
  {  
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) ) 
    {
      break;
    }
//
//  If the text begins with "!!" then we want to print it,
//  but skipping the first three characters.
//
    if ( s_begin ( line, "!!" ) || 
         s_begin ( line, "cc" ) || 
         s_begin ( line, "CC" ) )
    {
//
//  Begin the HTML List Item, and go into Bold face.
//
      cout << "        <li>\n";
      cout << "          <b>";
//
//  Print out all the characters, until you see a space.
//
      inc = 3;
      while ( line[inc] != ' ' && line[inc] && inc < line_length ) 
      {
        cout << line[inc];
        inc = inc + 1;
      }
//
//  Terminate Bold face, then print the rest of the characters,
//  and end the List Item.
//
      cout << "</b>";
      while ( inc < line_length )
      {
        cout << line[inc];
        inc = inc + 1;
      }
      cout << "\n";
      cout << "        </li>\n";
    }
  }
//
//  End the HTML "Unnumbered List" and Paragraph.
//
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      You can go up one level to <a href = \"../f_src.html\">\n";
  cout << "      the FORTRAN90 source codes</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <i>\n";
  cout << "      Last revised on " << date_string;
  cout << "    </i>\n";
  cout << "\n";
  cout << "    <!-- John Burkardt -->\n";
  cout << "\n";
  cout << "  </body>\n";
  cout << "\n";
  cout << "  <!-- Initial HTML skeleton created by HTMLINDEX. -->\n";
  cout << "\n";
  cout << "</html>\n";
//
//  Close the file.
//
  file_in.close ( );

  return;    
}
//****************************************************************************80

void handle_m ( string file_name, string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE_M processes a single MATLAB file.
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
//    Input, string FILE_NAME, the name of the file to be processed.
//
//    Input, string S1, the name of the file, minus extension.
//
//    Input, string S2, the same as S1, but capitalized.
//
{
  ifstream file_in;
  string line;
  int line_length;
  int inc;
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
//
//  Write the header.
//
  cout << "<html>\n";
  cout << "\n";
  cout << "  <head>\n";
  cout << "    <title>\n";
  cout << "      " << s2 << " - title goes here.\n";
  cout << "    </title>\n";
  cout << "  </head>\n";
  cout << "\n";
  cout << "  <body bgcolor=\"#EEEEEE\" link=\"#CC0000\"";
  cout << " alink=\"#FF3300\" vlink=\"#000055\">\n";
  cout << "\n";
  cout << "    <h1 align = \"center\">\n";
  cout << "      ";
  cout << s2;
  cout << " <br> heading goes here.\n";
  cout << "    </h1>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b>\n";
  cout << "      is a MATLAB (program/library) which\n";
  cout << "      (description goes here).\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Usage:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Licensing:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      The computer code and data files made available on this\n";
  cout << "      web page are distributed under\n";
  cout << "      <a href = \"../../txt/gnu_lgpl.txt\">the GNU LGPL license.</a>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Languages:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <b>" << s2 << "</b> is available in\n";
  cout << "      <a href = \"../../m_src/" << s1 << "/" << s1 << ".html\">a Matlab version</a>.\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Related Data and Programs:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Reference:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ol>\n";
  cout << "        <li>\n";
  cout << "          \n";
  cout << "        </li>\n";
  cout << "      </ol>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Source Code:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << ".m\">";
  cout << s1;
  cout << ".m</a>, the source code.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <h3 align = \"center\">\n";
  cout << "      Examples and Tests:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_test.m\">";
  cout << s1;
  cout << "_test.m</a>,\n";
  cout << "          a sample calling program.\n";
  cout << "        </li>\n";
  cout << "        <li>\n";
  cout << "          <a href = \"";
  cout << s1;
  cout << "_test_output.txt\">";
  cout << s1;
  cout << "_test_output.txt</a>,\n";
  cout << "          the output file.\n";
  cout << "        </li>\n";
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
//
//  Start the HTML Paragraph, and "Unnumbered List".
//
  cout << "    <h3 align = \"center\">\n";
  cout << "      List of Routines:\n";
  cout << "    </h3>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      <ul>\n";
//
//  Get another line of text.
//
  while ( 1 ) 
  {  
    getline ( file_in, line );
    line_length = line.length ( );

    if ( file_in.eof ( ) ) 
    {
      break;
    }
//
//  If the text begins with "%%" then we want to print it,
//  but skipping the first three characters.
//
    if ( s_begin ( line, "%%" ) )
    {
//
//  Begin the HTML List Item, and go into Bold face.
//
      cout << "        <li>\n";
      cout << "          <b>";
//
//  Print out all the characters, until you see a space.
//
      inc = 3;
      while ( line[inc] != ' ' && line[inc] && inc < line_length ) 
      {
        cout << line[inc];
        inc = inc + 1;
      }
//
//  Terminate Bold face, then print the rest of the characters,
//  and end the List Item.
//
      cout << "</b>";
      while ( inc < line_length )
      {
        cout << line[inc];
        inc = inc + 1;
      }
      cout << "\n";
      cout << "        </li>\n";
    }
  }
//
//  End the HTML "Unnumbered List" and Paragraph.
//
  cout << "      </ul>\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <p>\n";
  cout << "      You can go up one level to <a href = \"../m_src.html\">\n";
  cout << "      the MATLAB source codes</a>.\n";
  cout << "    </p>\n";
  cout << "\n";
  cout << "    <hr>\n";
  cout << "\n";
  cout << "    <i>\n";
  cout << "      Last revised on " << date_string;
  cout << "    </i>\n";
  cout << "\n";
  cout << "    <!-- John Burkardt -->\n";
  cout << "\n";
  cout << "  </body>\n";
  cout << "\n";
  cout << "  <!-- Initial HTML skeleton created by HTMLINDEX. -->\n";
  cout << "\n";
  cout << "</html>\n";
//
//  Close the file.
//
  file_in.close ( );

  return;    
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
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

string s_cap ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_CAP capitalizes all the characters in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be capitalized.
//
//    Output, string S_CAP, the capitalized string.
//
{
  int i;
  int s_length;
  string s2;

  s_length = s.length ( );
  s2 = s;

  for ( i = 0; i < s_length; i++ ) 
  {
    s2[i] = ch_cap ( s2[i] );
  }

  return s2;
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
