#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#include "constants.hpp"
#include "vec_utils.hpp"
#include "netcdf_utils.hpp"

using namespace std;

#ifdef _MACOS
	#include <GLUT/glut.h>
#elif _LINUX
	#include <GL/glut.h>
#endif

int main ( int argc, char *argv[] );
void display ( );
void mouse ( int btn, int state, int x, int y );
void gl_init ( );
void myReshape ( int w, int h );
void spin_image ( );
void timestamp ( );

void build_connectivity();
void setup_ranges();
void build_range(int id);

void rescale_cells();
void rescale_vertices();

void draw_cells();
void draw_triangles();
void draw_edges();

void draw_cell_lines();
void draw_triangle_lines();
void draw_edge_lines();

void color_cells();
void color_triangles();
void color_edges();

void arrowKeys( int a_keys, int x, int y );
void keyPressed( unsigned char key, int x, int y );
void translateView ( double updown, double leftright);
void polarView( double distance, double twist, double elevation, double azimuth );
void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b);
//
//  Global data.
//
string filename;
static GLint axis = 1;
GLint window;
int drawing = 0;
int draw_lines = 1;

double line_factor = 1.002;
double range_factor = 0.50;

double projUpDown = 0.0;
double projLeftRight = 0.0;
double projDistance	= 3.0;
double projTwist		= 0.0;
double projElevation	= 0.0;
double projAzimuth		= 0.0;

bool spinning = false;
static GLfloat theta[3] = { 0.0, 0.0, 0.0 };
double theta_speed = 0.020;

int ntime;
int nvertlevels;
int ncells;
int nvertices;
int nedges;
int maxedges;
int pixel_height;
int pixel_width;

int edge_field, cell_field, vertex_field;
int cur_level = 0;
int cur_time = 0;

double *xcell;
double *ycell;
double *zcell;
double *xvertex;
double *yvertex;
double *zvertex;

int *verticesoncell;
int *verticesonedge;
int *cellsonvertex;
int *cellsonedge;
int *nedgesoncell;

double *cell_values;
double *triangle_values;
double *edge_values;

vector<GLfloat> cell_cells;
vector<GLfloat> vertex_cells;
vector<GLfloat> edge_cells;
vector<GLfloat> cell_lines;
vector<GLfloat> vertex_lines;
vector<GLfloat> edge_lines;
vector<GLfloat> cell_colors;
vector<GLfloat> vertex_colors;
vector<GLfloat> edge_colors;

vector< vector<double> > ranges; // 0 - min, 1 - max

double xyz_center[3];
double xyz_max[3];
double xyz_min[3];
double xyz_range[3];
double xyz_scale;

int main ( int argc, char *argv[] ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    MAIN is the main program for TRIANGULATION_FACES.
	//
	//  Discussion:
	//
	//    This program reads certain information from an MPAS NETCDF grid file,
	//    and displays the faces of the triangulation.
	//
	//  Usage:
	//
	//    triangulation_faces file.nc
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
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	//
	int cell;
	int edge;
	int i;
	int v;

	cout << "\n";
	timestamp ( );

	cout << "\n";
	cout << "MPAS_DRAW:\n";
	cout << "  C++ version\n";
	cout << "  Read an MPAS NETCDF grid file\n";
	cout << "  Visualize the mpas grid/output file.\n";
	cout << "\n";
	cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
	//
	//  If the input file was not specified, get it now.
	//
	if ( argc <= 1 )
	{
		cout << "\n";
		cout << "MPAS_DRAW:\n";
		cout << "  Please enter the MPAS NETCDF grid filename.\n";

		cin >> filename;
	}
	else
	{
		filename = argv[1];
	}

	build_connectivity();
	setup_ranges();
	color_cells();
	color_triangles();
	color_edges();

	//
	//  Hand things over to OpenGL.
	//
	glutInit ( &argc, argv );
	glutInitDisplayMode ( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );
	glutInitWindowSize ( 600, 600 );
	glutInitWindowPosition ( 0, 0 );
	window = glutCreateWindow ( filename.c_str ( ) );
	glutReshapeFunc ( myReshape );
	glutDisplayFunc ( display );
	glutIdleFunc ( spin_image );
	glutMouseFunc ( mouse );
	glutKeyboardFunc( keyPressed );
	glutSpecialFunc( arrowKeys );

	//
	//  Enable hidden surface removal.
	//
	glEnable ( GL_DEPTH_TEST );

	gl_init ( );

	glutMainLoop ( );

	//
	//  Things that won't actually happen because we never return from glutMainLoop:
	//
	//
	//  Terminate.
	//
	cout << "\n";
	cout << "TRIANGULATION_FACES::\n";
	cout << "  Normal end of execution.\n";

	cout << "\n";
	timestamp ( );

	return 0;
}/*}}}*/
void display ( ){/*{{{*/

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
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();						// moves to center of screen
	translateView ( projUpDown, projLeftRight );
	polarView( projDistance, projTwist, projElevation, projAzimuth );

	switch(drawing){
		case 0:
			draw_triangles();
			if(draw_lines)
				draw_triangle_lines();
			break;
		case 1:
			draw_cells();
			if(draw_lines)
				draw_cell_lines();
			break;
		case 2:
			draw_edges();
			if(draw_lines)
				draw_edge_lines();
			break;
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
}/*}}}*/
void mouse ( int btn, int state, int x, int y ){/*{{{*/
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
	if ( ( btn == GLUT_LEFT_BUTTON   && state == GLUT_DOWN ) ||
			( btn == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN ) ||
			( btn == GLUT_RIGHT_BUTTON  && state == GLUT_DOWN ) ) {
		if ( spinning )	{
			spinning = false;
			theta_speed = 0.0;
		} else {
			spinning = true;
			theta_speed = 0.020;
			axis = axis + 1;
		}
	}

	axis = axis % 3;

	return;
}/*}}}*/
void gl_init ( ){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    gl_init initializes OpenGL state variables dealing with viewing and attributes.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license.
	//
	//  Modified:
	//
	//    01 January 2010
	//
	//  Author:
	//
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	GLfloat line_width;
	GLfloat point_size;
	//
	//  Set the background to WHITE.
	//
	glClearColor ( 1.0, 1.0, 1.0, 1.0 );
	//
	//  Antialiasing.
	//
	glDepthFunc( GL_LESS );
	glEnable(GL_DEPTH_TEST);
	glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glHint ( GL_LINE_SMOOTH_HINT, GL_DONT_CARE );
	glShadeModel( GL_SMOOTH );

	glEnableClientState( GL_VERTEX_ARRAY );
	glEnableClientState( GL_COLOR_ARRAY );


	glEnable ( GL_POINT_SMOOTH );
	glEnable ( GL_LINE_SMOOTH );
	//
	//  The default point size is 1.0.
	//
	if ( ncells <= 400 ) {
		point_size = 16.0;
	} else if ( ncells <= 800 ) {
		point_size = 8.0;
	} else if ( ncells <= 1600 ) {
		point_size = 4.0;
	} else if ( ncells <= 3200 ) {
		point_size = 2.0;
	} else {
		point_size = 1.0;
	}
	glPointSize ( point_size );
	//
	//  The default line width is 1.0.
	//
	if ( ncells <= 1600 ){
		line_width = 1.0;
	} else {
		line_width = 0.1;
	}
	glLineWidth ( line_width );

	return;
}/*}}}*/
void myReshape ( int w, int h ){/*{{{*/

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
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//
	//  Reference:
	//
	//    Edward Angel,
	//    Interactive Computer Graphics:
	//    A Top-Down Approach with OpenGL,
	//    Second Edition,
	//    Addison Wesley, 2000.
	//
	//

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

	glViewport ( 0, 0, w, h );
	glMatrixMode ( GL_PROJECTION );
	glLoadIdentity ( );

	gluPerspective( 45.0f, (GLfloat)w / (GLfloat)h, 0.1f, 5.0f );
	//glMatrixMode ( GL_MODELVIEW );

	return;
}/*}}}*/
void spin_image ( ){/*{{{*/

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
	theta[axis] = theta[axis] + theta_speed;

	if ( 360.0 < theta[axis] )
	{
		theta[axis] = theta[axis] - 360.0;
	}
	glutPostRedisplay ( );

	return;
}/*}}}*/
void timestamp ( ){/*{{{*/

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
}/*}}}*/

void build_connectivity(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    BUILD_CONNECTIVITY builds the connectivity arrays for Display to draw
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
	//    John Burkardt, Geoff Womeldorff, Doug Jacobsen
	//

	int i, j;
	int v1, v2;
	int c1, c2, c3, dc;


	//
	//  Get sizes.
	//
	ntime = netcdf_mpas_read_dim(filename, "Time" );
	nvertlevels = netcdf_mpas_read_dim(filename, "nVertLevels");
	ncells = netcdf_mpas_read_dim (filename, "nCells");
	nvertices = netcdf_mpas_read_dim(filename, "nVertices");
	nedges = netcdf_mpas_read_dim(filename, "nEdges");
	maxedges = netcdf_mpas_read_dim(filename, "maxEdges");

	cout << "\n";
	cout << "  The number of time steps NTIME   = " << ntime << "\n";
	cout << "  The number of vertical levels NVERTLEVELS = " << nvertlevels << "\n";
	cout << "  The number of cells NCELLS       = " << ncells << "\n";
	cout << "  The number of vertices NVERTICES =  " << nvertices << "\n";
	cout << "  The number of edges NEDGES =  " << nedges << "\n";

	xcell = new double[ncells];
	ycell = new double[ncells];
	zcell = new double[ncells];

	xvertex = new double[nvertices];
	yvertex = new double[nvertices];
	zvertex = new double[nvertices];

	verticesoncell = new int[maxedges*ncells];
	verticesonedge = new int[2*nedges];
	cellsonvertex = new int[3*nvertices];
	cellsonedge = new int[2*nedges];
	nedgesoncell = new int[ncells];

	netcdf_mpas_read_xyzcell ( filename, ncells, xcell, ycell, zcell );
	netcdf_mpas_read_xyzvertex ( filename, nvertices, xvertex, yvertex, zvertex );

	rescale_cells();
	rescale_vertices();

	netcdf_mpas_read_verticesoncell ( filename, maxedges, ncells, verticesoncell );
	netcdf_mpas_read_verticesonedge ( filename, nedges, verticesonedge );
	netcdf_mpas_read_cellsonvertex ( filename, nvertices, cellsonvertex );
	netcdf_mpas_read_cellsonedge ( filename, nedges, cellsonedge );
	netcdf_mpas_read_nedgesoncell ( filename, ncells, nedgesoncell );
	//
	//  connectivity is 1 based.  Fix that.
	//
	for ( i = 0; i < nvertices; i++ )
	{
		for ( j = 0; j < 3; j++ )
		{
			cellsonvertex[i*3+j]--;
		}
	}

	for( i = 0; i < ncells; i++){
		for ( j = 0; j < maxedges; j++){
			verticesoncell[i*maxedges + j]--;
		}
	}

	for( i = 0; i < nedges; i++){
		for( j = 0; j < 2; j++){
			cellsonedge[i*2 + j]--;
			verticesonedge[i*2 + j]--;
		}
	}


	for(i = 0; i < nvertices; i++){
		c1 = cellsonvertex[i*3];
		c2 = cellsonvertex[i*3 + 1];
		c3 = cellsonvertex[i*3 + 2];

		if(c1 == -1){
			if(c2 == -1){
				dc = c3;
			} else {
				dc = c2;
			}
		} else {
			dc = c1;
		}

		if(c1 == -1){
			c1 = dc;
		}
		if(c2 == -1){
			c2 = dc;
		}
		if(c3 == -1){
			c3 = dc;
		}
		vertex_cells.push_back(xcell[c1]);
		vertex_cells.push_back(ycell[c1]);
		vertex_cells.push_back(zcell[c1]);

		vertex_cells.push_back(xcell[c2]);
		vertex_cells.push_back(ycell[c2]);
		vertex_cells.push_back(zcell[c2]);

		vertex_cells.push_back(xcell[c3]);
		vertex_cells.push_back(ycell[c3]);
		vertex_cells.push_back(zcell[c3]);

		vertex_lines.push_back(xcell[c1]*line_factor);
		vertex_lines.push_back(ycell[c1]*line_factor);
		vertex_lines.push_back(zcell[c1]*line_factor);
		vertex_lines.push_back(xcell[c2]*line_factor);
		vertex_lines.push_back(ycell[c2]*line_factor);
		vertex_lines.push_back(zcell[c2]*line_factor);

		vertex_lines.push_back(xcell[c2]*line_factor);
		vertex_lines.push_back(ycell[c2]*line_factor);
		vertex_lines.push_back(zcell[c2]*line_factor);
		vertex_lines.push_back(xcell[c3]*line_factor);
		vertex_lines.push_back(ycell[c3]*line_factor);
		vertex_lines.push_back(zcell[c3]*line_factor);

		vertex_lines.push_back(xcell[c3]*line_factor);
		vertex_lines.push_back(ycell[c3]*line_factor);
		vertex_lines.push_back(zcell[c3]*line_factor);
		vertex_lines.push_back(xcell[c1]*line_factor);
		vertex_lines.push_back(ycell[c1]*line_factor);
		vertex_lines.push_back(zcell[c1]*line_factor);
	}

	for(i = 0; i < ncells; i++){
		for(j = 0; j < nedgesoncell[i]; j++){
			v1 = verticesoncell[i*maxedges + j%nedgesoncell[i]];
			v2 = verticesoncell[i*maxedges + (j+1)%nedgesoncell[i]];

			cell_cells.push_back(xcell[i]);
			cell_cells.push_back(ycell[i]);
			cell_cells.push_back(zcell[i]);
			cell_cells.push_back(xvertex[v1]);
			cell_cells.push_back(yvertex[v1]);
			cell_cells.push_back(zvertex[v1]);
			cell_cells.push_back(xvertex[v2]);
			cell_cells.push_back(yvertex[v2]);
			cell_cells.push_back(zvertex[v2]);

			cell_lines.push_back(xvertex[v1]*line_factor);
			cell_lines.push_back(yvertex[v1]*line_factor);
			cell_lines.push_back(zvertex[v1]*line_factor);

			cell_lines.push_back(xvertex[v2]*line_factor);
			cell_lines.push_back(yvertex[v2]*line_factor);
			cell_lines.push_back(zvertex[v2]*line_factor);
		}
	}

	for(i = 0; i < nedges; i++){
		c1 = cellsonedge[i*2];
		c2 = cellsonedge[i*2+1];
		v1 = verticesonedge[i*2];
		v2 = verticesonedge[i*2+1];

		if(c1 == -1){
			c1 = c2;
		} else if(c2 == -1){
			c2 = c1;
		}

		if(v1 == -1){
			v1 = v2;
		} else if (v2 == -1){
			v2 = v1;
		}
		edge_cells.push_back(xcell[c1]);
		edge_cells.push_back(ycell[c1]);
		edge_cells.push_back(zcell[c1]);

		edge_cells.push_back(xvertex[v1]);
		edge_cells.push_back(yvertex[v1]);
		edge_cells.push_back(zvertex[v1]);

		edge_cells.push_back(xvertex[v2]);
		edge_cells.push_back(yvertex[v2]);
		edge_cells.push_back(zvertex[v2]);

		edge_cells.push_back(xcell[c2]);
		edge_cells.push_back(ycell[c2]);
		edge_cells.push_back(zcell[c2]);

		edge_cells.push_back(xvertex[v2]);
		edge_cells.push_back(yvertex[v2]);
		edge_cells.push_back(zvertex[v2]);

		edge_cells.push_back(xvertex[v1]);
		edge_cells.push_back(yvertex[v1]);
		edge_cells.push_back(zvertex[v1]);

		edge_lines.push_back(xcell[c1]*line_factor);
		edge_lines.push_back(ycell[c1]*line_factor);
		edge_lines.push_back(zcell[c1]*line_factor);
		edge_lines.push_back(xvertex[v1]*line_factor);
		edge_lines.push_back(yvertex[v1]*line_factor);
		edge_lines.push_back(zvertex[v1]*line_factor);

		edge_lines.push_back(xvertex[v1]*line_factor);
		edge_lines.push_back(yvertex[v1]*line_factor);
		edge_lines.push_back(zvertex[v1]*line_factor);
		if(c1 == c2){
			edge_lines.push_back(xvertex[v2]*line_factor);
			edge_lines.push_back(yvertex[v2]*line_factor);
			edge_lines.push_back(zvertex[v2]*line_factor);
		} else {
			edge_lines.push_back(xcell[c2]*line_factor);
			edge_lines.push_back(ycell[c2]*line_factor);
			edge_lines.push_back(zcell[c2]*line_factor);

			edge_lines.push_back(xcell[c2]*line_factor);
			edge_lines.push_back(ycell[c2]*line_factor);
			edge_lines.push_back(zcell[c2]*line_factor);
			edge_lines.push_back(xvertex[v2]*line_factor);
			edge_lines.push_back(yvertex[v2]*line_factor);
			edge_lines.push_back(zvertex[v2]*line_factor);
		}

		edge_lines.push_back(xvertex[v2]*line_factor);
		edge_lines.push_back(yvertex[v2]*line_factor);
		edge_lines.push_back(zvertex[v2]*line_factor);
		edge_lines.push_back(xcell[c1]*line_factor);
		edge_lines.push_back(ycell[c1]*line_factor);
		edge_lines.push_back(zcell[c1]*line_factor);
	}

	cell_values = new double[ncells];
	triangle_values = new double[nvertices];
	edge_values = new double[nedges];

	color_cells();
	color_triangles();
	color_edges();

	delete [] xcell;
	delete [] ycell;
	delete [] zcell;
	delete [] xvertex;
	delete [] yvertex;
	delete [] zvertex;
	delete [] verticesoncell;
	delete [] verticesonedge;
	delete [] cellsonvertex;
	delete [] cellsonedge;
}/*}}}*/
void setup_ranges(){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    SETUP_RANGES...
    //
	int num_vars;

	num_vars = netcdf_mpas_read_num_vars(filename);

	ranges.resize(num_vars);

	return;
}/*}}}*/
void build_range(int id){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    BUILD_RANGE...
    //
	int num_items;
	double max, min;
	vector<double> temp_data;

	if(ranges[id].size() == 0){
		cout << "Building min-max range for coloring." << endl;
		num_items = netcdf_mpas_field_num_items(filename, id);

		temp_data.resize(num_items);

		netcdf_mpas_read_full_field(filename, id, &temp_data[0]);

		r8vec_min_max(num_items, &temp_data[0], min, max);

		ranges[id].push_back(min);
		ranges[id].push_back(max);
		cout << "Range build complete." << endl;
	}
}/*}}}*/

void rescale_cells(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    RESCALE_CELLS scales xcell, ycell, and zcell arrays to have a norm of 1
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
	//    John Burkardt, Doug Jacobsen
	//

	int i;
	int cell;
	xyz_min[0] = r8vec_min ( ncells, xcell );
	xyz_max[0] = r8vec_max ( ncells, xcell );
	xyz_min[1] = r8vec_min ( ncells, ycell );
	xyz_max[1] = r8vec_max ( ncells, ycell );
	xyz_min[2] = r8vec_min ( ncells, zcell );
	xyz_max[2] = r8vec_max ( ncells, zcell );

	xyz_range[0] = xyz_max[0] - xyz_min[0];
	xyz_range[1] = xyz_max[1] - xyz_min[1];
	xyz_range[2] = xyz_max[2] - xyz_min[2];

	if ( xyz_range[0] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_cells(): - Fatal error!\n";
		cout << "  The X data range is 0.\n";
		exit ( 1 );
	}

	if ( xyz_range[1] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_cells(): - Fatal error!\n";
		cout << "  The Y data range is 0.\n";
		exit ( 1 );
	}
	if ( xyz_range[2] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_cells(): - Fatal error!\n";
		cout << "  The Z data range is 0.\n";
		exit ( 1 );
	}

	xyz_scale = 0.0;
	for (i = 0; i < 3; i++ )
	{
		xyz_center[i] = ( xyz_min[i] + xyz_max[i] ) / 2.0;
		xyz_scale = std::max ( xyz_scale, ( xyz_max[i] - xyz_min[i] ) / 2.0 );
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

}/*}}}*/
void rescale_vertices(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    RESCALE_VERTICES scales xvertex, yvertex, and zvertex arrays to have a norm of 1
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
	//    John Burkardt, Doug Jacobsen
	//

	int i;
	int vertex;
	xyz_min[0] = r8vec_min ( nvertices, xvertex );
	xyz_max[0] = r8vec_max ( nvertices, xvertex );
	xyz_min[1] = r8vec_min ( nvertices, yvertex );
	xyz_max[1] = r8vec_max ( nvertices, yvertex );
	xyz_min[2] = r8vec_min ( nvertices, zvertex );
	xyz_max[2] = r8vec_max ( nvertices, zvertex );

	xyz_range[0] = xyz_max[0] - xyz_min[0];
	xyz_range[1] = xyz_max[1] - xyz_min[1];
	xyz_range[2] = xyz_max[2] - xyz_min[2];

	if ( xyz_range[0] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_vertices(): - Fatal error!\n";
		cout << "  The X data range is 0.\n";
		exit ( 1 );
	}

	if ( xyz_range[1] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_vertices(): - Fatal error!\n";
		cout << "  The Y data range is 0.\n";
		exit ( 1 );
	}
	if ( xyz_range[2] == 0.0 )
	{
		cout << "\n";
		cout << "rescale_vertices(): - Fatal error!\n";
		cout << "  The Z data range is 0.\n";
		exit ( 1 );
	}

	xyz_scale = 0.0;
	for (i = 0; i < 3; i++ )
	{
		xyz_center[i] = ( xyz_min[i] + xyz_max[i] ) / 2.0;
		xyz_scale = std::max ( xyz_scale, ( xyz_max[i] - xyz_min[i] ) / 2.0 );
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
	for ( vertex = 0; vertex < nvertices; vertex++ )
	{
		xvertex[vertex] = ( xvertex[vertex] - xyz_center[0] ) / xyz_scale;
		yvertex[vertex] = ( yvertex[vertex] - xyz_center[1] ) / xyz_scale;
		zvertex[vertex] = ( zvertex[vertex] - xyz_center[2] ) / xyz_scale;
	}

}/*}}}*/

void draw_triangles ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_TRIANGLES draws the Delaunay triangles
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glColorPointer( 3, GL_FLOAT, 0, &vertex_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &vertex_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, vertex_cells.size()/3);

	return;
}/*}}}*/
void draw_cells ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_CELLS draws the Voronoi cells
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glColorPointer( 3, GL_FLOAT, 0, &cell_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &cell_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, cell_cells.size()/3);

	return;
}/*}}}*/
void draw_edges ( ){/*{{{*/

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_EDGES draws quadrilaterals to represent edge values
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glColorPointer( 3, GL_FLOAT, 0, &edge_colors[0] );
	glVertexPointer( 3, GL_FLOAT, 0, &edge_cells[0] );
	glDrawArrays( GL_TRIANGLES, 0, edge_cells.size()/3);

	return;
}/*}}}*/

void draw_triangle_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_TRIANGLE_LINES draws the lines for the Delaunay triangle edges
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &vertex_lines[0] );
	glDrawArrays( GL_LINES, 0, vertex_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/
void draw_cell_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_CELL_LINES draws the lines for the Voronoi cell edges.
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &cell_lines[0] );
	glDrawArrays( GL_LINES, 0, cell_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/
void draw_edge_lines(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    DRAW_EDGE_LINES draws the lines for the edge "edges".
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	glDisableClientState( GL_COLOR_ARRAY );
	glColor3f ( LINE_R, LINE_G, LINE_B );
	glVertexPointer( 3, GL_FLOAT, 0, &edge_lines[0] );
	glDrawArrays( GL_LINES, 0, edge_lines.size()/3 );
	glEnableClientState( GL_COLOR_ARRAY);

	return;
}/*}}}*/

void color_cells(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_CELLS builds the color array used to display the Voronoi cells
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
	//    Doug Jacobsen
	//

	int i, j, k, o;
	double max, min;
	long *dims;
	int num_dims = 0;
	int cell_dim;
	int vert_dim;
	int time_dim;
	int num_items;
	float h, s, v;
	float r, g, b;

	cell_colors.clear();

	s = 1.0;
	v = 1.0;
	if(cell_field == 0){
		for(i = 0; i < ncells; i++){
			o = nedgesoncell[i];

			for(j = 0; j < 3*o; j++){
				if(o < 5){
					cell_colors.push_back(LIL_R);
					cell_colors.push_back(LIL_G);
					cell_colors.push_back(LIL_B);
				} else if(o == 5){
					cell_colors.push_back(PENTA_R);
					cell_colors.push_back(PENTA_G);
					cell_colors.push_back(PENTA_B);
				} else if(o == 6){
					cell_colors.push_back(HEXA_R);
					cell_colors.push_back(HEXA_G);
					cell_colors.push_back(HEXA_B);
				} else if(o == 7){
					cell_colors.push_back(HEPTA_R);
					cell_colors.push_back(HEPTA_G);
					cell_colors.push_back(HEPTA_B);
				} else {
					cell_colors.push_back(OTHER_R);
					cell_colors.push_back(OTHER_G);
					cell_colors.push_back(OTHER_B);
				}
			}
		}
	} else {
		netcdf_mpas_print_field_info(filename,cell_field);
		netcdf_mpas_read_field(filename, cell_field, cell_values, cur_time, cur_level);

		build_range(cell_field);

		min = ranges[cell_field].at(0);
		max = ranges[cell_field].at(1);

		cout << "Range: " << min << " to " << max << endl;

		for(i = 0; i < ncells; i++){
			if((max-min) != 0.0){
				h = (cell_values[i] - min)/(max-min) * range_factor;
			} else {
				h = (cell_values[i] - min)/1.0 * range_factor;
			}

			o = nedgesoncell[i];

			hsv_to_rgb(h, s, v, r, g, b);

			for(j = 0; j < 3*o; j++){
				cell_colors.push_back(r);
				cell_colors.push_back(g);
				cell_colors.push_back(b);
			}
		}
	}

	return;
}/*}}}*/
void color_triangles(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_TRIANGLES builds the color array used to display the Delaunay triangles
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
	//    Doug Jacobsen
	//

	int i, j;
	double max, min;
	float h, s, v;
	float r, g, b;

	vertex_colors.clear();

	s = 1.0;
	v = 1.0;
	if(vertex_field == 0){
		for(i = 0; i < nvertices; i++){
			for(j = 0; j < 3; j++){
				vertex_colors.push_back(HEXA_R);
				vertex_colors.push_back(HEXA_G);
				vertex_colors.push_back(HEXA_B);
			}
		}
	} else {
		netcdf_mpas_print_field_info(filename,vertex_field);
		netcdf_mpas_read_field(filename, vertex_field, triangle_values, cur_time, cur_level);

		build_range(vertex_field);

		max = ranges[vertex_field][0];
		min = ranges[vertex_field][1];

		cout << "Range: " << min << " to " << max << endl;

		for(i = 0; i < nvertices; i++){
			if(max-min != 0.0){
				h = (triangle_values[i] - min)/(max-min)*range_factor;
			}else{
				h = (triangle_values[i] - min)/1.0 * range_factor;
			}

			hsv_to_rgb(h, s, v, r, g, b);

			for(j = 0; j < 3; j++){
				vertex_colors.push_back(r);
				vertex_colors.push_back(g);
				vertex_colors.push_back(b);
			}
		}
	}

	return;
}/*}}}*/
void color_edges(){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    COLOR_EDGES builds the color array used to display the edge quads
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
	//    Doug Jacobsen
	//
	int i, j, o;
	double max, min;
	float h, s, v;
	float r, g, b;

	edge_colors.clear();

	s = 1.0;
	v = 1.0;
	if(edge_field == 0){
		for(i = 0; i < nedges; i++){
			for(j = 0; j < 6; j++){
				edge_colors.push_back(HEXA_R);
				edge_colors.push_back(HEXA_G);
				edge_colors.push_back(HEXA_B);
			}
		}
	} else {
		netcdf_mpas_print_field_info(filename, edge_field);
		netcdf_mpas_read_field(filename, edge_field, edge_values, cur_time, cur_level);
		r8vec_min_max(nedges, edge_values, min, max);

		build_range(edge_field);

		min = ranges[edge_field][0];
		max = ranges[edge_field][1];

		cout << "Range: " << min << " to " << max << endl;

		for(i = 0; i < nedges; i++){
			if(max-min != 0.0){
				h = (edge_values[i] - min)/(max-min)*range_factor;
			} else {
				h = (edge_values[i] - min)/1.0 * range_factor;
			}

			hsv_to_rgb(h, s, v, r, g, b);

			for(j = 0; j < 6; j++){
				edge_colors.push_back(r);
				edge_colors.push_back(g);
				edge_colors.push_back(b);
			}
		}
	}
}/*}}}*/

void arrowKeys( int a_keys, int x, int y ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    ARROWKEYS processes special keys, for example the arrow keys, and the f-keys
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
	//    Geoff Womeldorff
	//

	switch( a_keys ) {
		case GLUT_KEY_UP:
			projElevation = (double)( ( (int)projElevation + 6 ) % 360 );
			break;

		case GLUT_KEY_DOWN:
			projElevation = (double) ( ( (int)projElevation - 6 ) % 360 );
			break;

		case GLUT_KEY_RIGHT:
			projAzimuth = (double) ( ( (int)projAzimuth + 6 ) % 360 );
			break;

		case GLUT_KEY_LEFT:
			projAzimuth = (double) ( ( (int)projAzimuth - 6 ) % 360 );
			break;

		case GLUT_KEY_F1:
			glutFullScreen();
			break;

		case GLUT_KEY_F2:
			glutReshapeWindow( kWindowWidth, kWindowHeight );
			break;

		default:
			break;

	}

	return;
}/*}}}*/
void keyPressed( unsigned char key, int x, int y ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    KEYPRESSED processes non-special keys, for example the letter keys
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
	//    Geoff Womeldorff, Doug Jacobsen
	//


	usleep( 100 );

	// Used to get codes of keys you don't know.
//	cout << "Key pressed: " << key << endl;
//	cout << "Code is: " << (int)key << endl;

	switch( key ) {

		case KEY_ESCAPE:
			glutDestroyWindow( window );
			glDisableClientState( GL_VERTEX_ARRAY );
			glDisableClientState( GL_COLOR_ARRAY );
			//glDisableClientState( GL_NORMAL_ARRAY );
			delete [] nedgesoncell;

			cout << "End of normal execution. Exiting." << endl;

			exit( 0 );
		case KEY_l:
			cur_level = (cur_level + 1) % nvertlevels;
			cout << "Current vertical level: " << cur_level+1 << " out of " << nvertlevels << endl;
			switch(drawing){
				case 0:
					color_triangles();
					break;
				case 1:
					color_cells();
					break;
				case 2:
					color_edges();
					break;
				default:
					break;
			}
			break;
		case KEY_t:
			cur_time = (cur_time + 1) % ntime;
			cout << "Current time level: " << cur_time+1 << " out of " << ntime << endl;
			switch(drawing){
				case 0:
					color_triangles();
					break;
				case 1:
					color_cells();
					break;
				case 2:
					color_edges();
					break;
				default:
					break;
			}
			break;
		case KEY_s:
			projLeftRight -= 0.1;
			break;
		case KEY_f:
			projLeftRight += 0.1;
			break;
		case KEY_d:
			projUpDown -= 0.1;
			break;
		case KEY_e:
			projUpDown += 0.1;
			break;
		case KEY_w:
			draw_lines = (draw_lines + 1) % 2;
			break;
		case KEY_c:
			drawing = ( drawing + 1 ) % 3;
			break;
		case KEY_COMMA:
			projDistance += 0.1;
//			if ( projDistance > 3.0 )
//				projDistance = 3.0;
			break;
		case KEY_PERIOD:
			projDistance -= 0.1;
//			if ( projDistance < 2.0 ) projDistance = 2.0;
			break;
		case KEY_v:
			switch(drawing){
				case 0:
					vertex_field = netcdf_mpas_list_nvertex_fields(filename);
					color_triangles();
					break;
				case 1:
					cell_field = netcdf_mpas_list_ncell_fields(filename);
					color_cells();
					break;
				case 2:
					edge_field = netcdf_mpas_list_nedge_fields(filename);
					color_edges();
					break;
				default:
					break;
			}
			break;
		case KEY_r:
			cout << "Resetting to default parameters." << endl;
			cur_time = 0;
			cur_level = 0;
			cell_field = 0;
			edge_field = 0;
			vertex_field = 0;
			projDistance = 3.0;
			projUpDown = 0.0;
			projLeftRight = 0.0;
			projTwist = 0.0;
			projElevation = 0.0;
			projAzimuth = 0.0;
			switch(drawing){
				case 0:
					color_triangles();
					break;
				case 1:
					color_cells();
					break;
				case 2:
					color_edges();
					break;
				default:
					break;
			}
			break;
		default:
			break;
	}

	return;
}/*}}}*/

void translateView ( double updown, double leftright){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    TRANSLATEVIEW Translates the drawing left, right, up, or down
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
	//    Doug Jacobsen
	//

	glTranslated ( leftright, updown, 0.0);
}/*}}}*/
void polarView( double distance, double twist, double elevation, double azimuth ) {/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    POLARVIEW Rotates the drawing (in two directions), and translates in and out
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
	//    Geoff Womeldorff
	//

	glTranslated( 0.0, 0.0, -1.0 * distance );
	glRotated( -1.0 * twist, 0.0, 0.0, 1.0 );
	glRotated( -1.0 * elevation, 1.0, 0.0, 0.0 );
	glRotated( azimuth, 0.0, 0.0, 1.0 );
}/*}}}*/

void hsv_to_rgb(float h, float s, float v, float& r, float& g, float& b){/*{{{*/
	//****************************************************************************80
	//
	//  Purpose:
	//
	//    HSV_TO_RGB convertes a hue-saturation-value triplet to a red-green-blue triplet
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
	//    Geoff Womeldorff, Doug Jacobsen
	//

	int i;
	float f;
	float m;
	float n;

	h *= 6.0;

	i = floor( h );
	f = h - i;

	if(!(i & 1)) f = 1 - f; // if i is even
	m = v * (1 - s);
	n = v * (1 - s * f);

	switch (i) {
		case 6:
		case 0:
			r = v; g = n; b = m;
			break;
		case 1:
			r = n; g = v; b = m;
			break;
		case 2:
			r = m; g = v; b = n;
			break;
		case 3:
			r = m; g = n; b = v;
			break;
		case 4:
			r = n; g = m; b = v;
			break;
		case 5:
			r = v; g = m; b = n;
			break;
	}

	return;
}/*}}}*/
