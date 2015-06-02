# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cstring>
# include <ctime>

using namespace std;

int main ( int argc, char** argv );
void test01 ( void );
void test02 ( void );
void test03 ( void );
void timestamp ( void );

//****************************************************************************80

int main ( int argc, char** argv )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MICE_READER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "MICE_READER:\n";
  cout << "  C++ version\n";
  cout << "  Read a file of strings, one line at a time.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "MICE_READER:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 reads the file using GETLINE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  ifstream input;
  string input_filename = "mice_file.txt";
  string line;
  int line_num;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Read a line at a time, using getline().\n";
  cout << "  Input terminates on EOF().\n";
  cout << "  getline() sees initial and trailing blanks.\n";
  cout << "  It sees blank lines.\n";

  input.open ( input_filename.c_str ( ) );

  cout << "\n";
  line_num = 0;

  while ( 1 )
  {
    getline ( input, line );
    if ( input.eof ( ) )
    {
      break;
    }
    line_num = line_num + 1;
    cout << line_num << ": \"" << line << "\"\n";
  }

  input.close ( );

  return;
}
//****************************************************************************80

void test02 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 reads the file one character at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  char c;
  char c_old;
  ifstream input;
  string input_filename = "mice_file.txt";

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Read one character at a time.\n";
  cout << "  The << operator won't return a blank or newline\n";
  cout << "  when asked to retrieve a character.\n";

  input.open ( input_filename.c_str ( ) );

  cout << "\n";
  cout << "\"";

  c = '\n';

  while ( 1 )
  {
    c_old = c;

    input >> c;
    
    if ( input.eof ( ) )
    {
      cout << "\"\n";
      break;
    }

    cout << c;
  }

  input.close ( );

  return;
}
//****************************************************************************80

void test03 ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 reads the file one string (word) at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  ifstream input;
  string input_filename = "mice_file.txt";
  string line;
  int word_num;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Read one string (word) at a time.\n";
  cout << "  The << operator won't return blanks\n";
  cout << "  when asked to retrieve a string.\n";
  cout << "\n";

  input.open ( input_filename.c_str ( ) );

  word_num = 1;

  while ( 1 )
  {
    input >> line;
    
    if ( input.eof ( ) )
    {
      break;
    }

    cout << word_num << ": \"" << line << "\"\n";

    word_num = word_num + 1;
  }

  input.close ( );

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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
