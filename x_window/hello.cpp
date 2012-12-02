# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <Xm/Xm.h>
# include <Xm/Label.h>

using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************

int main ( int argc, char *argv[] )

//****************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for HELLO.
//
//  Discussion:
//
//    HELLO displays the string "Hello, world".
//
//  Compilation:
//
//    IRIX: cc hello.c -lXaw -lXt
//
//  Reference:
//
//    Douglas Young,
//    Object Oriented Programming with C++ and OSF/Motif,
//    Prentice Hall, 1995
//
//  Modified:
//
//    31 January 2003
//
{
  XtAppContext app;
  Arg args[10];
  Widget label;
  int n;
  Widget shell;
  XmString xmstr;
//
//  Initialize Xt and create a top level shell widget;
//
  shell = XtAppInitialize ( &app, "Hello", NULL,
    0, &argc, argv, NULL, NULL, 0 );
//
//  Create a compound string to dsiplay the Hello message.
//
  xmstr = XmStringCreateLocalized ( "Hello World!" );
//
//  Create a label widget to display the string.
//
  n = 0;
  XtSetArg ( args[n], XmNlabelString, xmstr );
  n = n + 1;
  label = XtCreateManagedWidget ( "label", xmLabelWidgetClass,
    shell, args, n );
//
//  Free the compound string.
//
  XmStringFree ( xmstr );
//
//  Realize the widgets.
//
  XtRealizeWidget ( shell );
//
//  Enter the event loop.
//
  XtAppMainLoop ( app );

  return ( 0 );
}
