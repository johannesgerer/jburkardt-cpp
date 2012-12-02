# include <cstdlib>
# include <iostream>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
bool ch_is_digit ( char c );
bool includes_date ( char *line, int monnum, int daynum );
bool includes_weekday ( char *line, int weekdaynum );
bool is_leap_year ( int year );
void next_day ( int dmy[3] );
int month_length ( int month, int year );
bool parse_command_line ( int argc, char *argv[], int *lookahead );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MY_CALENDAR.
//
//  Discussion:
//
//    MY_CALENDAR prints out "important" dates from a calendar file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 January 2011
//
//  Author:
//
//    John Burkardt
//
{
# define DMAX 100

  int daynum[DMAX];
  bool DEBUG = false;
  int dmy[3];
  bool error;
  ifstream file_in;
  string file_in_name = "/Users/jburkardt/public_html/calendar";
  int i;
  char line[255];
  struct tm *loctim;
  int lookahead = 0;
  int monnum[DMAX];
  time_t systim;
  bool VERBOSE = false;
  int w;
  int wdaynum;
  int yearnum[DMAX];
//
//  Print hello.
//
  if ( VERBOSE || DEBUG )
  {
    timestamp ( );
    cout << "\n";
    cout << "MY_CALENDAR:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  Simple calendar reminder utility.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
    cout << "\n";
  }
//
//  Get the system time.
//
  systim = time ( NULL );
//
//  Convert system time to local time.
//
  loctim = localtime ( &systim );
//
//  Pick out useful items from the local time structure.
//
  daynum[0] = (*loctim).tm_mday;
  monnum[0] = (*loctim).tm_mon + 1;
  yearnum[0] = 1900 + (*loctim).tm_year;

  wdaynum = (*loctim).tm_wday;

  if ( DEBUG )
  {
    cout << "\n";
    cout << "DEBUG: WDAYNUM = " << wdaynum << "\n";
  }
//
//  Sunday
//
  if ( wdaynum == 0 )
  {
    lookahead = 3;
  }
//
//  Monday through Thursday.
//
  else if ( 1 <= wdaynum && wdaynum <= 4 )
  {
    lookahead = 2;
  }
//
//  Friday through Saturday.
//
  else
  {
    lookahead = 4;
  }

  if ( DEBUG )
  {
    cout << "\n";
    cout << "DEBUG: Lookahead = " << lookahead << "\n";
  }
//
//  See if user wants to change lookahead.
//
  error = !parse_command_line ( argc, argv, &lookahead );

  if ( error )
  {
    cout << "\n";
    cout << "MY_CALENDAR - Fatal error!\n";
    cout << "  Problems reading the command line arguments.\n";
    return 1;
  }

  if ( DMAX < lookahead )
  {
    cout << "\n";
    cout << "MY_CALENDAR - Fatal error!\n";
    cout << "  Can't look that far ahead!\n";
  }
//
//  Compute the D/M/Y of the next LOOKAHEAD days.
//
  for ( i = 1; i < lookahead; i++ )
  {
    dmy[0] = daynum[i-1];
    dmy[1] = monnum[i-1];
    dmy[2] = yearnum[i-1];

    next_day ( dmy );

    daynum[i] = dmy[0];
    monnum[i] = dmy[1];
    yearnum[i] = dmy[2];
  }
//
//  DEBUG: Print out the date arrays.
//
  if ( DEBUG )
  {
    for ( i = 0; i < lookahead; i++ )
    {
      cout << "  D/M/Y = " << daynum[i] << "/"
        << monnum[i] << "/" << yearnum[i] << "\n";
    }
  }
//
//  Open the calendar file.
//
  file_in.open ( file_in_name.c_str ( ) );

  if ( !file_in )
  {
    cout << "\n";
    cout << "MY_CALENDAR - Fatal error!\n";
    cout << "  Could not open the user calendar file:\""
      << file_in_name << "\"\n";
    return 1;
  }
//
//  While not the end-of-file,
//
  while ( 1 )
  {
//
//  Get the next line of input,
//
    file_in.getline ( line, sizeof ( line ) );

    if ( file_in.eof ( ) )
    {
      break;
    }
//
//  Search the line for variations on dates between today and
//  today + LOOKAHEAD working days.
//
    i = 0;
    w = ( wdaynum + i ) % 7;
    if ( includes_weekday ( line, w ) )
    {
      cout << line << "\n";
    }

    for ( i = 0; i < lookahead; i++ )
    {
      if ( includes_date ( line, monnum[i], daynum[i] ) )
      {
        cout << line << "\n";
        break;
      }
    }

  }

  file_in.close ( );

  if ( VERBOSE )
  {
    cout << "\n";
    cout << "MY_CALENDARr:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }
  return 0;
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

bool includes_date ( char *line, int monnum, int daynum )

//****************************************************************************80
//
//  Purpose:
//
//    INCLUDES_DATE reports if a line contains a given date.
//
//  Discussion:
//
//    The date is specified as a day and  month.
//
//    For instance, given DAYNUM=5 and MONNUM=3, it accepts
//    lines containing:
//
//      5 March, March 5, 05 March and March 05.
//
//    Big bug:
//
//      If March 1 is being sought, the code accepts March 10,
//      March 11, ..., March 19 as well!
//
//    Try  4 Feb, Feb 4, 04 Feb, Feb 04.
//    Could allow 04/02, 04/2, 4/02, 4/2 or the reversal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *LINE, the line to check.
//
//    Input, int MONNUM, the month number, between 1 and 12.
//
//    Input, int DAYNUM, the day number, between 1 and 31.
//
//    Output, bool INCLUDES_DATE, is true if the line contains the date.
{
  char date[4][7];
  int i;
  int ihi;
  int l;
  char *match_add;
  char month[12][4] =
  {
    "Jan", "Feb", "Mar", "Apr", "May", "Jun",
    "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
  };
//
//  Store in the DATE array all the acceptable variations
//  of DAYNUM, MONNUM that represent the date.
//
  sprintf ( &date[0][0], "%d %s",daynum, month[monnum-1] );
  sprintf ( &date[1][0], "%s %d",month[monnum-1], daynum );

  if ( daynum < 10 )
  {
    sprintf ( &date[2][0], "0%d %s",daynum, month[monnum-1] );
    sprintf ( &date[3][0], "%s 0%d",month[monnum-1], daynum );
    ihi = 3;
  }
  else
  {
    ihi = 1;
  }
//
//  Do any of these variations occur in the LINE?
//
  for ( i = 0; i <= ihi; i = i+1 )
  {

    match_add = strstr ( line, &date[i][0] );

    if ( match_add != NULL )
    {
//
//  Weed out cases where "May 3" is matched by "May 30, I went home".
//
      l = strlen ( &date[i][0] );

      if ( !ch_is_digit( *(match_add+l) ) )
      {
        return true;
      }
    }

  }

  return false;
}
//****************************************************************************80

bool includes_weekday ( char *line, int weekdaynum )

//****************************************************************************80
//
//  Purpose:
//
//    INCLUDES_WEEKDAY checks if a line refers to a "generic" weekday.
//
//  Discussion:
//
//    For instance, given WEEKDAYNUM = 3 (meaning Wednesday) this
//    routine will return 1 if the line contains the string
//
//      "Wednesdays"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *LINE, the line to be examined.
//
//    Input, int WEEKDAYNUM, the day of the week to be checked for.
//
//    Output, bool INCLUDES_WEEKDAY, is true if the line contains a
//    reference to the weekday.
{
  char *match_add;
  char  weekdays[7][11] =
  {
    "Sundays",
    "Mondays",
    "Tuesdays",
    "Wednesdays",
    "Thursdays",
    "Fridays",
    "Saturdays"
  };

  weekdaynum = weekdaynum % 7;

  match_add = strstr ( line, &weekdays[weekdaynum][0] );

  if ( match_add == NULL )
  {
    return false;
  }
  else
  {
    return true;
  }
}
//****************************************************************************80

bool is_leap_year ( int year )

//****************************************************************************80
//
//  Purpose:
//
//    IS_LEAP_YEAR determines if a given year was a leap year.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int YEAR, the year.
//
//    Output, bool IS_LEAP_YEAR, is true if the year was a leap year.
//
{
  if ( year%400 == 0 )
  {
    return true;
  }
  else if ( year%100 == 0 )
  {
    return false;
  }
  else if ( year%4 == 0 )
  {
    return true;
  }
  else
  {
    return false;
  }
}
//****************************************************************************80

int month_length ( int month, int year )

//****************************************************************************80
//
//  Purpose:
//
//    MONTH_LENGTH returns the number of days in a given month.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int MONTH, the month.
//
//    Input, int YEAR, the year.
//
//    Output, int MONTH_LENGTH, the number of days in the month.
//
{
  int mdays[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
//
//  Adjust number of days in February leap years;
//
  if ( month == 2 && is_leap_year ( year ) )
  {
    return ( mdays[month-1]+1 );
  }
  else
  {
    return ( mdays[month-1] );
  }
}
//****************************************************************************80

void next_day ( int dmy[3] )

//****************************************************************************80
//
//  Purpose:
//
//    NEXT_DAY increments the date by one day;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int DMY[3], the day, month and year of a date.
//
{
  dmy[0] = dmy[0] + 1;
//
//  If we've exceeded the days in the month, decrement the day
//  and increment the month.
//
  if ( dmy[0] > month_length ( dmy[1], dmy[2] ) )
  {
    dmy[0] = dmy[0] - month_length ( dmy[1], dmy[2] );
    dmy[1] = dmy[1] + 1;
//
//  If we've exceeded the months in the year, decrement the month
//  and increment the year.
//
    if ( dmy[1] > 12 )
    {
      dmy[1] = dmy[1] - 12;
      dmy[2] = dmy[2] + 1;
    }
  }
}
//****************************************************************************80

bool parse_command_line ( int argc, char *argv[], int *lookahead )

//****************************************************************************80
//
//  Purpose:
//
//    PARSE_COMMAND_LINE parses the command line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ARGC, the command line argument counter.
//
//    Input, char *ARGV[], the command line argument array.
//
//    Output, int *LOOKAHEAD, the calendar lookahead, if specified.
//
//    Output, bool PARSE_COMMAND_LINE, is true if there was no error.
{
//
//  This macro sets the character pointer x to the address of the string
//  containing the parameter's argument.  The argument may be specified as
//  either: -xargument or -x argument.
//
  while ( --argc > 0 && **++argv == '-' )

    {
    switch( (*argv)[1] )
      {
//
//  N: Specifying number of days of lookahead.
//
      case 'n':
      case 'N':

        if ( (*argv)[2] )
        {
          sscanf ( (*argv)+2, "%d", lookahead );

        }
        else if ( argc > 1 )
        {
          argc = argc - 1;
          argv = argv + 1;
          sscanf ( *argv, "%d", lookahead );
        }

        break;

      default:
        cout << "\n";
        cout << "CALENDAR - Fatal error!\n";
        cout << "  Invalid arguments.\n";
        return false;
      }
    }

  return true;
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
//    03 October 2003
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
