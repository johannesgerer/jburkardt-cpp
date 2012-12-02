# include <cstdlib>
# include <iostream>
//
//  UNISTD is needed for the SLEEP command.
//
# include <unistd.h>

# include <GLUT/glut.h>
//# include <GL/glx.h>
//# include <GL/gl.h>

using namespace std;
//
// Use GLX_RED_SIZE 1 to get the deepest buffer with 1 red bit.
//
static int attributeListSgl[]  =
{
  GLX_RGBA,
  GLX_RED_SIZE,   1,
  GLX_GREEN_SIZE, 1,
  GLX_BLUE_SIZE,  1,
  None 
};
//
// Use GLX_DOUBLEBUFFER in case single buffering is not supported.
//
static int attributeListDbl[]  = 
{
  GLX_RGBA,
  GLX_DOUBLEBUFFER,
  GLX_RED_SIZE,    1,
  GLX_GREEN_SIZE, 1,
  GLX_BLUE_SIZE,  1,
  None 
};

int main ( int argc, char *argv[] );
static Bool WaitForNotify ( Display *d,  XEvent *e, char  *arg );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
// 
//  Purpose:
//
//    YELLOW_WINDOW opens a window and sets it to yellow/..
//
//  Discussion:
//
//    This is the minimum code required to create an RGBA-format, OpenGL-
//    compatible X window and clear it to yellow. 
//
//    The OpenGL application interface for X Windows is known as GLX.
//
{
  Display *dpy;
  XVisualInfo *vi;
  Colormap cmap;
  XSetWindowAttributes swa;
  Window win;
  GLXContext cx;
  XEvent event;
  int swap_flag = GL_FALSE;
//
// Get a connection. 
//
  dpy = XOpenDisplay ( 0 );
// 
// Get an appropriate visual.
//
  vi = glXChooseVisual ( dpy, DefaultScreen(dpy), attributeListSgl );

  if ( vi == NULL )
  {
    vi = glXChooseVisual ( dpy, DefaultScreen(dpy), attributeListDbl );
    swap_flag = GL_TRUE;
  }
//
// Create a GLX context.
//
  cx = glXCreateContext ( dpy, vi, 0,  GL_TRUE );
//
// Create a color map.
//
  cmap = XCreateColormap ( dpy, RootWindow(dpy, vi->screen ),
    vi->visual,  AllocNone );
// 
// Create a window. 
//
  swa.colormap = cmap;
  swa.border_pixel = 0;

  swa.event_mask = StructureNotifyMask;

  win = XCreateWindow ( dpy, RootWindow(dpy, vi->screen), 0, 0, 100, 100,
    0, vi->depth,  InputOutput, vi->visual,
    CWBorderPixel|CWColormap|CWEventMask,  &swa );

  XMapWindow ( dpy, win );
  XIfEvent ( dpy, &event, WaitForNotify, (char*)win );
//
// Connect the context to the window.
//
  glXMakeCurrent ( dpy, win, cx );
//
// Clear the buffer.
//
  glClearColor ( 1, 1, 0, 1 );
  glClear ( GL_COLOR_BUFFER_BIT );
  glFlush ( );

  if ( swap_flag )
  { 
    glXSwapBuffers ( dpy, win );
  }
// 
// Wait a while. 
//
  sleep ( 10 );
}
//****************************************************************************80

static Bool WaitForNotify ( Display *d,  XEvent *e, char  *arg ) 

//****************************************************************************80
{
  return ( e->type == MapNotify ) && ( e->xmap.window == (Window)arg );
}
