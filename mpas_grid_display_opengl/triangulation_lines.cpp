# include <cstdlib>
# include <ctime>
# include <cmath>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "netcdfcpp.h"
//
//  This is the include statement I need for Mac OS X.
//
# include <GLUT/glut.h>

//# include <GL/glut.h>
//# include <GL/freeglut.h>

int main ( int argc, char *argv[] );
void display ( );
void mouse ( int btn, int state, int x, int y );
void myinit ( );
void myReshape ( int w, int h );
void netcdf_mpas_read_cellsonedge ( string filename, int nedges,
  int cellsonedge[] );
int netcdf_mpas_read_ncells ( string filename );
int netcdf_mpas_read_nedges ( string filename );
void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[],
  double ycell[], double zcell[] );
double r8_huge ( );
double r8_max ( double x, double y );
double r8vec_max ( int n, double r8vec[] );
double r8vec_min ( int n, double r8vec[] );
void spin_image ( );
void timestamp ( );
//
//  Global data.
//
  static GLint axis = 2;
  int *cellsonedge;
  int ncells;
  int nedges;
  int pixel_height;
  int pixel_width;
  bool spinning = true;
  static GLfloat theta[3] = { 0.0, 0.0, 0.0 };
  double theta_speed = 0.020;
  double *xcell;
  double xyz_center[3];
  double xyz_max[3];
  double xyz_min[3];
  double xyz_range[3];
  double xyz_scale;
  double *ycell;
  double *zcell;

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGULATION_LINES.
//
//  Discussion:
//
//    This program reads certain information from an MPAS NETCDF grid file,
//    and displays the lines of the triangulation.
//
//  Usage:
//
//    triangulation_lines file.nc
//
//    where
//
//    * file.nc is an MPAS NETCDF grid file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  int cell;
  int edge;
  int i;
  string filename;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TRIANGULATION_LINES:\n";
  cout << "  C++ version\n";
  cout << "  Read an MPAS NETCDF grid file\n";
  cout << "  Display the faces of the triangulation.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 )
  {
    cout << "\n";
    cout << "TRIANGULATION_LINES:\n";
    cout << "  Please enter the MPAS NETCDF grid filename.\n";

    cin >> filename;
  }
  else
  {
    filename = argv[1];
  }
//
//  Get sizes.
//
  ncells = netcdf_mpas_read_ncells ( filename );
  nedges = netcdf_mpas_read_nedges ( filename );

  cout << "\n";
  cout << "  The number of cells NCELLS = " << ncells << "\n";
  cout << "  The number of edges EDGES =  " << nedges << "\n";
//
//  Get the cell center coordinates.
//
  xcell = new double[ncells];
  ycell = new double[ncells];
  zcell = new double[ncells];

  netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );
//
//  Get the cell-to-cell connectivity in cellsOnEdge.
//
  cellsonedge = new int[2*nedges];

  netcdf_mpas_read_cellsonedge ( filename, nedges, cellsonedge );
//
//  Worry about data scaling.
//
  xyz_min[0] = r8vec_min ( ncells, xcell );
  xyz_max[0] = r8vec_max ( ncells, xcell );
  xyz_min[1] = r8vec_min ( ncells, ycell );
  xyz_max[1] = r8vec_max ( ncells, ycell );
  xyz_min[2] = r8vec_min ( ncells, zcell );
  xyz_max[2] = r8vec_max ( ncells, zcell );

  xyz_range[0] = xyz_max[0] - xyz_min[0];
  xyz_range[1] = xyz_max[1] - xyz_min[1];
  xyz_range[2] = xyz_max[2] - xyz_min[2];

  cout << "\n";
  cout << "  Minimum: " << xyz_min[0]
       << "  " << xyz_min[1]
       << "  " << xyz_min[2] << "\n";
  cout << "  Maximum: " << xyz_max[0]
       << "  " << xyz_max[1]
       << "  " << xyz_max[2] << "\n";
  cout << "  Range:   " << xyz_range[0]
       << "  " << xyz_range[1]
       << "  " << xyz_range[2] << "\n";

  if ( xyz_range[0] == 0.0 )
  {
    cout << "\n";
    cout << "TRIANGULATION_LINES: - Fatal error!\n";
    cout << "  The X data range is 0.\n";
    exit ( 1 );
  }

  if ( xyz_range[1] == 0.0 )
  {
    cout << "\n";
    cout << "TRIANGULATION_LINES: - Fatal error!\n";
    cout << "  The Y data range is 0.\n";
    exit ( 1 );
  }
  if ( xyz_range[2] == 0.0 )
  {
    cout << "\n";
    cout << "TRIANGULATION_LINES: - Fatal error!\n";
    cout << "  The Z data range is 0.\n";
    exit ( 1 );
  }

  xyz_scale = 0.0;
  for (i = 0; i < 3; i++ )
  {
    xyz_center[i] = ( xyz_min[i] + xyz_max[i] ) / 2.0;
    xyz_scale = r8_max ( xyz_scale, ( xyz_max[i] - xyz_min[i] ) / 2.0 );
  }
//
//  A sphere doesn't need this rescaling.
//  A box does!
//
//xyz_scale = sqrt ( 3.0 ) * xyz_scale;
//
//  Translate the data so it is centered.
//  Scale the data so it fits in the unit cube.
//
  for ( cell = 0; cell < ncells; cell++ )
  {
    xcell[cell] = ( xcell[cell] - xyz_center[0] ) / xyz_scale;
    ycell[cell] = ( ycell[cell] - xyz_center[1] ) / xyz_scale;
    zcell[cell] = ( zcell[cell] - xyz_center[2] ) / xyz_scale;
  }
//
//  CELLSONEDGE is 1 based.  Fix that.
//
  for ( edge = 0; edge < nedges; edge++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      cellsonedge[i+edge*2] = cellsonedge[i+edge*2] - 1;
    }
  }
//
//  Hand things over to OpenGL.
//
  glutInit ( &argc, argv );
  glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
  glutInitWindowSize ( 1000, 1000 );
  glutInitWindowPosition ( 0, 0 );
  glutCreateWindow ( filename.c_str ( ) );
  glutReshapeFunc ( myReshape );
  glutDisplayFunc ( display );
  glutIdleFunc ( spin_image );
  glutMouseFunc ( mouse );
//
//  Enable hidden surface removal.
//
  glEnable ( GL_DEPTH_TEST );

  myinit ( );

  glutMainLoop ( );
//
//  Things that won't actually happen because we never return from glutMainLoop:
//
  delete [] cellsonedge;
  delete [] xcell;
  delete [] ycell;
  delete [] zcell;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGULATION_LINES::\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void display ( )

//****************************************************************************80
//
//  Purpose:
//
//    DISPLAY generates the graphics output.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  int cell;
  int edge;
  int i;
  float p[3];
//
//  Clear the window.
//
  glClear ( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glLoadIdentity ( );

  glRotatef ( theta[0], 1.0, 0.0, 0.0 );
  glRotatef ( theta[1], 0.0, 1.0, 0.0 );
  glRotatef ( theta[2], 0.0, 0.0, 1.0 );
//
//  Draw lines in BLUE.
//  The linewidth was set in myinit.
//
  glColor3f ( 0.0, 0.0, 1.0 );

  for ( edge = 0; edge < nedges; edge++ )
  {
    glBegin ( GL_LINE_STRIP );

    for ( i = 0; i < 2; i++ )
    {
      cell = cellsonedge[i+edge*2];
      p[0] = ( float ) xcell[cell];
      p[1] = ( float ) ycell[cell];
      p[2] = ( float ) zcell[cell];
      glVertex3fv ( p );
    }
    glEnd ( );
  }
//
//  Draw points in RED.
//  The point size was set in myinit.
//
  glColor3f ( 1.0, 0.0, 0.0 );

  for ( cell = 0; cell < ncells; cell++ )
  {
    glBegin ( GL_POINTS );

    p[0] = ( float ) xcell[cell];
    p[1] = ( float ) ycell[cell];
    p[2] = ( float ) zcell[cell];
    glVertex3fv ( p );

    glEnd ( );
  }
//
//  Clear all the buffers.
//
  glFlush ( );
//
//  Switch between the two buffers for fast animation.
//
  glutSwapBuffers ( );

  return;
}
//****************************************************************************80

void mouse ( int btn, int state, int x, int y )

//****************************************************************************80
//
//  Purpose:
//
//    MOUSE determines the response to mouse input.
//
//  Discussion:
//
//    The original routine assumed the user had a three button mouse, and
//    dedicated one axis to each.
//
//    Since Apple prefers the esthetics of a one button mouse, we're forced
//    to live with that choice.  This routine alternately pauses rotation,
//    or increments the rotation axis by 1, no matter which button is pushed.
//
//  Modified:
//
//    30 December 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  if ( ( btn == GLUT_LEFT_BUTTON   && state == GLUT_DOWN ) ||
       ( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) ||
       ( btn == GLUT_RIGHT_BUTTON  && state == GLUT_DOWN ) )
  {
    if ( spinning )
    {
      spinning = false;
      theta_speed = 0.0;
    }
    else
    {
      spinning = true;
      theta_speed = 0.020;
      axis = axis + 1;
    }
  }

  axis = axis % 3;

  return;
}
//****************************************************************************80

void myinit ( )

//****************************************************************************80
//
//  Purpose:
//
//    MYINIT initializes OpenGL state variables dealing with viewing and attributes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
{
  GLfloat line_width;
  GLfloat point_size;
//
//  Set the background to WHITE.
//
  glClearColor ( 1.0, 1.0, 1.0, 1.0 );
//
//  Antialiasing.
//
  glEnable ( GL_POINT_SMOOTH );
  glEnable ( GL_LINE_SMOOTH );
  glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  glHint ( GL_LINE_SMOOTH_HINT, GL_DONT_CARE );
//
//  The default point size is 1.0.
//
  if ( ncells <= 400 )
  {
    point_size = 8.0;
  }
  else if ( ncells <= 800 )
  {
    point_size = 4.0;
  }
  else if ( ncells <= 1600 )
  {
    point_size = 2.0;
  }
  else
  {
    point_size = 1.0;
  }
  glPointSize ( point_size );
//
//  The default line width is 1.0.
//
  if ( ncells <= 1600 )
  {
    line_width = 2.0;
  }
  else
  {
    line_width = 1.0;
  }

  glLineWidth ( line_width );

  return;
}
//****************************************************************************80

void myReshape ( int w, int h )

//****************************************************************************80
//
//  Purpose:
//
//    MYRESHAPE determines the window mapping.
//
//  Modified:
//
//    30 December 2008
//
//  Author:
//
//    Edward Angel
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  glViewport ( 0, 0, w, h );
  glMatrixMode ( GL_PROJECTION );
  glLoadIdentity ( );

  if ( w <= h )
  {
    glOrtho (
      -1.05, +1.05,
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,
      -10.0, 10.0 );
  }
  else
  {
    glOrtho (
      - 1.05 * ( GLfloat ) h / ( GLfloat ) w, + 1.05 * ( GLfloat ) h / ( GLfloat ) w,
      - 1.05, + 1.05,
      -10.0, 10.0 );
  }

  glMatrixMode ( GL_MODELVIEW );

  return;
}
//****************************************************************************80

void netcdf_mpas_read_cellsonedge ( string filename, int nedges,
  int cellsonedge[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_CELLSONEDGE gets the cellsOnEdge information.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NEDGES, the number of edges.
//
//    Output, int CELLSONEDGE[2*NEDGES];
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "cellsOnEdge" );
  (*var_id).get ( &cellsonedge[0], nedges, 2 );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

int netcdf_mpas_read_ncells ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NCELLS gets the number of cells.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_NCELLS, the value of NCELLS.
//
{
  int ncells;
  long int ncells_size;
  NcDim *ncells_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get NCELLS, which is a NETCDF dimension.
//
  ncells_id = ncid.get_dim ( "nCells" );

  ncells_size = (*ncells_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  ncells = ( int ) ncells_size;

  return ncells;
}
//****************************************************************************80

int netcdf_mpas_read_nedges ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_NEDGES gets the number of edges.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Output, int NETCDF_MPAS_READ_NEDGES, the value of NEDGES.
//
{
  int nedges;
  long int nedges_size;
  NcDim *nedges_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get NCELLS, which is a NETCDF dimension.
//
  nedges_id = ncid.get_dim ( "nEdges" );

  nedges_size = (*nedges_id).size ( );
//
//  Close the file.
//
  ncid.close ( );

  nedges = ( int ) nedges_size;

  return nedges;
}
//****************************************************************************80

void netcdf_mpas_read_xyzcell ( string filename, int ncells, double xcell[],
  double ycell[], double zcell[] )

//****************************************************************************80
//
//  Purpose:
//
//    NETCDF_MPAS_READ_XYZCELL reads xCell, yCell, zCell.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Russ Rew, Glenn Davis, Steve Emmerson, Harvey Davies, Ed Hartne,
//    The NETCDF User's Guide,
//    Unidata Program Center, March 2009.
//
//  Parameters:
//
//    Input, string NC_FILENAME, the name of the NETCDF file to examine.
//
//    Input, int NCELLS, the number of nodes.
//
//    Output, double XCELL[NCELLS], YCELL[NCELLS], ZCELL[NCELLS], the
//    coordinates of the nodes.
//
{
  NcVar *var_id;
//
//  Open the file.
//
  NcFile ncid ( filename.c_str ( ), NcFile::ReadOnly );
//
//  Get the variable values.
//
  var_id = ncid.get_var ( "xCell" );
  (*var_id).get ( &xcell[0], ncells );
  var_id = ncid.get_var ( "yCell" );
  (*var_id).get ( &ycell[0], ncells );
  var_id = ncid.get_var ( "zCell" );
  (*var_id).get ( &zcell[0], ncells );
//
//  Close the file.
//
  ncid.close ( );

  return;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  if ( n <= 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_MAX - Fatal error!\n";
    cerr << "  Vector size N <= 0.\n";
    exit ( 1 );
  }

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

void spin_image ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPIN_IMAGE adjusts the angle of rotation and redisplays the picture.
//
//  Modified:
//
//    15 December 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Angel,
//    Interactive Computer Graphics:
//    A Top-Down Approach with OpenGL,
//    Second Edition,
//    Addison Wesley, 2000.
//
{
  theta[axis] = theta[axis] + theta_speed;

  if ( 360.0 < theta[axis] )
  {
    theta[axis] = theta[axis] - 360.0;
  }
  glutPostRedisplay ( );

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
