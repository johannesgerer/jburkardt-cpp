# include <cstdlib>
# include <cmath>
# include <iostream>
# include <fstream>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int pcl_to_pure ( char *file_in_name, char *file_out_name );
void s_tab_blank ( char *s );
void timestamp ( void );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    PCL_READ reads a PCL file for Arabidopsis and writes out a pure data file.
//
//  Discussion:
//
//    The PCL file has a (long) first line containing titles.  There follow
//    N lines, each containing 3 labels, and 14 floating values, separated
//    by TAB characters.
//
//    The output file contains N lines, each containing 14 floating values,
//    separated by spaces.
//
//  Usage:
//
//    pcl_read
//     in which case the input and output files will be prompted for;
//
//    pcl_read file1
//      in which case the output file will be prompted for;
//
//    pcl_read file1 file2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2003
//
//  Author:
//
//    John Burkardt
//
{
  char file_in_name[81];
  char file_out_name[81];
  int iarg;
  int status;
  char *string;

  timestamp ( );
  cout << "\n";
  cout << "PCL_READ:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Read an Arabidopsis PCL file.\n";
  cout << "  Write a copy containing just the numeric data.\n";

  iarg = 0;
//
//  Determination of the input file name.
//
  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "PCL_READ:\n";
    cout << "  Please enter the input PCL file name.\n";

    cin.getline ( file_in_name, sizeof ( file_in_name ) );

  }
  else
  {
    strcpy ( file_in_name, argv[1] );
  }
  cout << "\n";
  cout << "PCL_READ:\n";
  cout << "  Data will be read from '" << file_in_name << "'.\n";
//
//  Determine the output file name.
//
  if ( argc <= 2 )
  {
    cout << "\n";
    cout << "PCL_READ:\n";
    cout << "  Please enter the output file name.\n";

    cin.getline ( file_out_name, sizeof ( file_out_name ) );
  }
  else
  {
    strcpy ( file_out_name, argv[2] );
  }

  cout << "\n";
  cout << "PCL_READ:\n";
  cout << "  Data will be written to '" << file_out_name << "'.\n";

  status = pcl_to_pure ( file_in_name, file_out_name );

  cout << "\n";
  cout << "PCL_READ:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return status;
}
//****************************************************************************80

int pcl_to_pure ( char *file_in_name, char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    PCL_TO_PURE copies numeric data from a PCL file to another file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *file_in_name, the name of the input file.
//
//    Input, char *file_out_name, the name of the output file.
//
{
  ifstream file_in;
  ofstream file_out;
  int i;
  char *s;
  int skip;
  int status;
  char string[1024];
  int text_num = 0;
  int width;
  float x[14];
//
//  Open the input file.
//
  file_in.open ( file_in_name );

  if ( !file_in )
  {
    cout << "\n";
    cout << "PCL_READ - Fatal error!\n";
    cout << "  Could not open the input file '"
      << file_in_name << "'!\n";
    return 1;
  }
//
//  Open the output file.
//
  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PCL_READ - Fatal error!\n";
    cout << "  Could not open the output file '"
      << file_out_name << "'!\n";
    return 1;
  }

  for ( ;; )
  {
//
//  Read a line from the file.
//
    file_in.getline ( string, sizeof ( string ) );

    if ( file_in.eof ( ) )
    {
      break;
    }

    text_num = text_num + 1;
//
//  First line is the title.
//
    if ( text_num == 1 )
    {
      cout << "\n";
      cout << "PCL_TO_PURE: Title:\n";
      cout << "  " << string << "\n";
      continue;
    }
//
//  All other lines have the form of 17 fields separated by tabs.
//  The first 3 fields are identifiers.  
//  The following 14 are the numeric fields we want.
//
//  Easiest thing to do is skip over the first three tabs.
//
    skip = 0;
    s = string;
    for ( i = 0; i < 3; i++ )
    {
      while ( *s != '\t' )
      {
        s = s + 1;
        skip = skip + 1;
      }

      s = s + 1;
      skip = skip + 1;
    }
//
//  Before you write the thing out, convert each tab to a space.
//
    s_tab_blank ( s );

    file_out << s << "\n";
  }

  file_in.close ( );
  file_out.close ( );

  cout << "\n";
  cout << "PCL_TO_PURE:\n";
  cout << "  Number of input lines read was " 
    << text_num << ".\n";

  return 0;
}
//****************************************************************************80

void s_tab_blank ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//   S_TAB_BLANK replaces each TAB character by a space.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, char *s, the string to be transformed.
//
{
  while ( *s != 0 )
  {
    if ( *s == '\t' )
    {
      *s = ' ';
    }
    s++;
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
//    21 August 2002
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
# define TIME_SIZE 29

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  if ( len != 0 ) 
  {
    cout << time_buffer << "\n";
  }

  return;
# undef TIME_SIZE
}
