# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <X11/Intrinsic.h>
# include <X11/StringDefs.h>
# include <X11/Xaw/Box.h>
# include <X11/Xaw/Label.h>
# include <X11/Xm/PushB.h>

using namespace std;

int main ( int argc, char *argv[] );
void Farewell ( Widget widget, XtPointer clientData, XtPointer callData );

//******************************************************************************

int main ( int argc, char *argv[] )

//******************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for GOODBYE.
//
//  Discussion:
//
//    GOODBYE puts up an X window with a labeled button. 
//
//  Compilation:
//
//    IRIX: cc goodbye.c -lXaw -lXm -lXt
//
//  Reference:
//
//    Paul Asente and Ralph Swick,
//    X Window System Toolkit,
//    Digital Press, 1990
//
//  Modified:
//
//    01 October 2002
//
{
  XtAppContext app;
  Widget box;
  Widget button;
  Arg buttonArgs[] = 
  {
    {XtNx, 10},
    {XtNy, 40},
    {XtNlabel, "Click and die!"} 
  };
  Widget label;
  Arg labelArgs[] = 
  {
    { XtNx, 10},
    { XtNy, 10},
    {XtNlabel, "Goodbye, world!"}
  };
  Widget toplevel;
//
//  Create a top level shell widget;
//
  toplevel = XtAppInitialize ( &app, "Goodbye", (XrmOptionDescList) NULL,
    0, &argc, argv, (String *) NULL, (ArgList) NULL, 0 );
//
//  Create the box widget.
//
  box = XtCreateManagedWidget ( "box", boxWidgetClass,
    toplevel, (Arg *) NULL, 0 );
//
//  Create the label widget.
//
  label = XtCreateManagedWidget ( "label", labelWidgetClass,
    box, labelArgs, XtNumber(labelArgs) );
//
//  Create the button widget.
//
  button = XtCreateManagedWidget ( "button", xmPushButtonWidgetClass,
    box, buttonArgs, XtNumber(buttonArgs) );
//
//  Associate the callback function.
//
  XtAddCallback ( button, XtNcallback, Farewell, NULL );
//
//  Realize the widgets.
//
  XtRealizeWidget ( toplevel );
//
//  Enter the event loop.
//
  XtAppMainLoop ( app );
}
//******************************************************************************

void Farewell ( Widget widget, XtPointer clientData, 
  XtPointer callData )

//******************************************************************************
//
//  Purpose:
//
//    FAREWELL is a callback function.
//
//  Modified:
//
//    01 October 2002
// 
//  Parameters:
//
//    ?put, Widget widget, ?
//
//    ?put, XtPointer clientData, ?
//
//    ?put, XtPointer callData, ?
//
{
  exit ( 0 );
}
