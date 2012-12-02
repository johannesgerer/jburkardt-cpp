# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <strings.h>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
#include <X11/Xaw/Text.h>

using namespace std;

char *concat_args ( int n, char *words[] );
int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MEMO displays its command line arguments in an X window.
//
//  Discussion:
//
//    The routine given in the book doesn't work on my SGI.  I had to switch
//    from the Xw library to the Xaw library, and God knows whose fault it is
//    now that the program is not working!
//
//  Compilation:
//
//    IRIX: cc memo.c -lXaw -lXt
//
//  Reference:
//
//    Douglas Young,
//    X Window Systems Programming and Applications with Xt,
//    Prentice Hall, 1989
//
//  Modified:
//
//    30 September 2002
//
{
  char *message;
  Widget msg_widget;
  int n;
  Widget toplevel;
  Arg wargs[1];
//
//  Initialize the intrinsics, and create a top level shell widget;
//
  toplevel = XtInitialize ( argv[0], "Memo", NULL, 0, &argc, argv );
//
//  Concatenate the remaining command line arguments.
//
  message = concat_args ( argc, argv );

  cout << "Program was invoked with argument string: '" 
       << message << "'\n";
//
//  If a message is given on the command line, use it in the
//  XtNstring argument for the widget.
//
  n = 0;
  if ( message != NULL ) 
  {
    XtSetArg ( wargs[n], XtNstring, message );
    n = n + 1;
  }
  cout << "Got past call to XtSetArg\n";
//
//  Create the message widget.
//
  msg_widget = XtCreateManagedWidget ( "message", textWidgetClass,
    toplevel, wargs, n );

  cout << "Got past call to XtCreateManagedWidget\n";
//
//  Realize the widgets.
//
//  This is the call this is causing the segmentation fault on the SGI.
//
  XtRealizeWidget ( toplevel );

  cout << "Got past call to XtRealizeWidget\n";
//
//  Enter the event loop.
//
  XtMainLoop ( );
  cout << "Got past call to XtMainLoop.\n";

  return 0;
}
//******************************************************************************

char *concat_args ( int n, char *words[] ) 

//******************************************************************************
//
//  Purpose:
//
//    CONCAT_ARGS concatenates a set of strings into one blank separated string.
//
//  Discussion:
//
//    Actually, this concatenates words 1 through N-1, skipping word 0!
//
//  Modified:
//
//    30 September 2002
//
//  Parameters:
//
//    Input, int N, the number of strings to concatenate.
//
//    Input, char *WORDS[], pointers to the strings.
//
//    Output, char *CONCAT_ARGS, a pointer to the concatenated string.
//
{
  char *buffer;
  int i;
  int len = 0;
//
//  If there are no arguments other than the program name, just return an empty string.
//
  if ( n <= 1 ) 
  {
    cout << "CONCAT_ARGS - Returning empty string.\n";
    return ( "" );
  }
//
//  Figure out the total length of the string, plus a space between each pair of words.
//
  for ( i = 1; i < n; i++ )
  {
    len = len + strlen ( words[i] );
  }

  len = len + ( n - 1 );

  cout << "CONCAT_ARGS: LEN = " << len << ".\n";
//
//  Allocate the buffer and initialize.
//
  buffer = ( char * ) malloc ( len + 1 );

  if ( buffer == NULL ) 
  {
    cout << "\n";
    cout << "CONCAT_ARGS - Fatal error!\n";
    cout << "  Unable to allocate memory.\n";
    exit ( 1 );
  }
//
//  BUFFER starts out as an empty string.
//
  buffer[0] = '\0';
//
//  Add each of the words to BUFFER.
//
  for ( i = 1; i < n; i++ )
  {
    if ( i > 1 ) 
    {
      strcat ( buffer, " " );
    }
    strcat ( buffer, words[i] );
  }

  return ( buffer );
}
