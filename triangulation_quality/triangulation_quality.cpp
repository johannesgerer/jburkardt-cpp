# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
void alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *alpha_min, double *alpha_ave, 
  double *alpha_area );
double arc_cosine ( double c );
void area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double &area_min, double &area_max, double &area_ratio,
  double &area_ave, double &area_std, int &area_negative, int &area_zero );
void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
int file_column_count ( string filename );
int file_row_count ( string filename );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int *i4mat_data_read ( string input_filename, int m, int n );
void i4mat_header_read ( string input_filename, int *m, int *n );
void i4mat_transpose_print ( int m, int n, int a[], string title );
void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
void mesh_base_zero ( int node_num, int element_order, 
  int element_num, int element_node[] );
void q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *q_min, double *q_max, double *q_ave, 
  double *q_area );
double r8_abs ( double x );
double r8_huge ( );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double *r8mat_data_read ( string input_filename, int m, int n );
void r8mat_header_read ( string input_filename, int *m, int *n );
void r8mat_transpose_print ( int m, int n, double a[], string title );
void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title );
int s_len_trim ( string s );
int s_to_i4 ( string s, int *last, bool *error );
bool s_to_i4vec ( string s, int n, int ivec[] );
double s_to_r8 ( string s, int *lchar, bool *error );
bool s_to_r8vec ( string s, int n, double rvec[] );
int s_word_count ( string s );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGULATION_QUALITY.
//
//  Discussion:
//
//    TRIANGULATION_QUALITY determines quality measures for a triangulation.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//    The three quality measures are:
//
//      ALPHA_MEASURE
//      AREA_MEASURE
//      Q_MEASURE
//
//    In each case, the ideal value of the quality measure is 1, and
//    the worst possible value is 0.
//
//    The program also prints out the geometric bandwidth, which is the
//    bandwidth of the adjacency matrix of the nodes.
//
//  Usage:
//
//    triangulation_quality prefix
//
//    where 'prefix' is the common filename prefix:
//
//    * prefix_nodes.txt contains the node coordinates,
//    * prefix_elements.txt contains the element definitions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Local parameters:
//
//    Local, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    lists the nodes that make up each element.
//
//    Local, int ELEMENT_NUM, the number of elements.
//
//    Local, int ELEMENT_ORDER, the order of the elements, either 3 or 6.
//
//    Local, int NODE_DIM, the spatial dimension.
//
//    Local, int NODE_NUM, the number of nodes.
//
//    Local, double NODE_XY[NODE_DIM*NODE_NUM], the point set.
//
{
  double alpha_area;
  double alpha_ave;
  double alpha_min;
  double area_ave;
  double area_max;
  double area_min;
  int area_negative;
  double area_ratio;
  double area_std;
  int area_zero;
  string element_filename;
  int *element_node;
  int element_num;
  int element_order;
  int iarg;
  int iargc;
  int m;
  int ml;
  int mu;
  int node_dim;
  string node_filename;
  int node_num;
  double *node_xy;
  string prefix;
  double q_area;
  double q_ave;
  double q_max;
  double q_min;
  double value;

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TRIANGULATION_QUALITY:\n";
  cout << "  C++ version:\n";
  cout << "  Compute triangulation quality measures.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
//
//  The first commandline argument is the filename prefix.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TRIANGULATION_QUALITY:\n";
    cout << "  Please enter the filename prefix.\n";

    cin >> prefix;
  }
  else 
  {
    prefix = argv[1];
  }
//
//  Create the filenames.
//
  node_filename = prefix + "_nodes.txt";
  element_filename = prefix + "_elements.txt";
//
//  Read the node data.
//
  r8mat_header_read ( node_filename, &node_dim, &node_num );

  cout << "\n";
  cout << "  Read the header of \"" << node_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension NODE_DIM = " << node_dim << "\n";
  cout << "  Number of nodes NODE_NUM  = " << node_num << "\n";

  node_xy = r8mat_data_read ( node_filename, node_dim, node_num );

  cout << "\n";
  cout << "  Read the data in \"" << node_filename << "\".\n";

  r8mat_transpose_print_some ( node_dim, node_num, node_xy, 1, 1, node_dim, 5, 
    "  First 5 nodes:" );
//
//  Read the triangulation data.
//
  i4mat_header_read ( element_filename, &element_order, 
    &element_num );

  cout << "\n";
  cout << " Read the header of \"" << element_filename << "\".\n";
  cout << "\n";
  cout << "  Element order      = " << element_order << "\n";
  cout << "  Number of elements = " << element_num << "\n";

  element_node = i4mat_data_read ( element_filename, 
    element_order, element_num );

  cout << "\n";
  cout << "  Read the data in \"" << element_filename << "\".\n";

  i4mat_transpose_print_some ( element_order, element_num, element_node, 
    1, 1, element_order, 10, "  First 10 elements:" );
//
//  Detect and correct 1-based node indexing.
//
  mesh_base_zero ( node_num, element_order, element_num, element_node );
//
//  Compute the quality measures.
//
  alpha_measure ( node_num, node_xy, element_order, element_num, 
    element_node, &alpha_min, &alpha_ave, &alpha_area );

  cout << "\n";
  cout << "  ALPHA compares the smallest angle against 60 degrees.\n";
  cout << "  Values of ALPHA range from 0 (extremely poor) to 1 (excellent).\n";
  cout << "  (The second figure is the same number in degrees.)\n";
  cout << "\n";
  cout << "  ALPHA_MIN  : minimum over all triangles = " << alpha_min 
       << "  " << 60.0 * alpha_min << "\n";
  cout << "  ALPHA_AVE  : average over all triangles = " << alpha_ave 
       << "  " << 60.0 * alpha_ave << "\n";
  cout << "  ALPHA_AREA : average weighted by area   = " << alpha_area 
       << "  " << 60.0 * alpha_area << "\n";

  area_measure ( node_num, node_xy, element_order, element_num, 
    element_node, area_min, area_max, area_ratio, area_ave, area_std,
    area_negative, area_zero );

  cout << "\n";
  cout << "  AREA compares the areas of the triangles.\n";
  cout << "  Values of AREA_RATIO range from 0 (extremely poor) to 1 (excellent).\n";
  cout << "\n";
  cout << "  AREA_MIN   : minimum area         = " << area_min << "\n";
  cout << "  AREA_MAX   : maximum area         = " << area_max << "\n";
  cout << "  AREA_RATIO : minimum/maximum area = " << area_ratio << "\n";
  cout << "  AREA_AVE   : average area         = " << area_ave << "\n";
  cout << "  AREA_STD   : standard deviation   = " << area_std << "\n";
  cout << "  AREA_NEG   : area < 0             = " << area_negative << "\n";
  cout << "  AREA_ZERO  : area = 0             = " << area_zero << "\n";

  q_measure ( node_num, node_xy, element_order, element_num, element_node,
    &q_min, &q_max, &q_ave, &q_area );

  cout << "\n";
  cout << "  Q is the ratio of 2 * inradius to outradius.\n";
  cout << "  Values of Q range from 0 (extremely poor) to 1 (excellent).\n";
  cout << "\n";
  cout << "  Q_MIN  : minimum Q                  = " << q_min << "\n";
  cout << "  Q_MAX  : maximum Q                  = " << q_max << "\n";
  cout << "  Q_AVE  : average Q                  = " << q_ave << "\n";
  cout << "  Q_AREA : average Q weighted by area = " << q_area << "\n";

  bandwidth_mesh ( element_order, element_num, element_node, &ml, &mu, &m );

  cout << "\n";
  cout << "  The geometric bandwidth          M = " << m << "\n";
//
//  Free memory.
//
  delete [] node_xy;
  delete [] element_node;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGULATION_QUALITY:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *alpha_min, double *alpha_ave, 
  double *alpha_area )

//****************************************************************************80
//
//  Purpose:
//
//    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
//
//  Discusion:
//
//    The ALPHA measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    We compute the minimum angle among all the triangles in the triangulated
//    dataset and divide by the maximum possible value (which, in degrees,
//    is 60).  The best possible value is 1, and the worst 0.  A good
//    triangulation should have an ALPHA score close to 1.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, real ( kind = 8 ) Z(2,N), the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
//    the triangulation.
//
//    Output, double *ALPHA_MIN, the minimum value of ALPHA over all
//    triangles.
//
//    Output, double *ALPHA_AVE, the value of ALPHA averaged over
//    all triangles.
//
//    Output, double *ALPHA_AREA, the value of ALPHA averaged over
//    all triangles and weighted by area.
//
{
  double a_angle;
  int a_index;
  double a_x;
  double a_y;
  double ab_len;
  double alpha;
  double area;
  double area_total;
  double b_angle;
  int b_index;
  double b_x;
  double b_y;
  double bc_len;
  double c_angle;
  int c_index;
  double c_x;
  double c_y;
  double ca_len;
  double pi = 3.141592653589793;
  int triangle;
  double value;

  *alpha_min = r8_huge ( );
  *alpha_ave = 0.0;
  *alpha_area = 0.0;
  area_total = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*triangle_order];
    b_index = triangle_node[1+triangle*triangle_order];
    c_index = triangle_node[2+triangle*triangle_order];

    a_x = z[0+(a_index-1)*2];
    a_y = z[1+(a_index-1)*2];
    b_x = z[0+(b_index-1)*2];
    b_y = z[1+(b_index-1)*2];
    c_x = z[0+(c_index-1)*2];
    c_y = z[1+(c_index-1)*2];

    area = 0.5 * r8_abs ( a_x * ( b_y - c_y ) 
                        + b_x * ( c_y - a_y ) 
                        + c_x * ( a_y - b_y ) );

    ab_len = sqrt ( pow ( a_x - b_x, 2 ) + pow ( a_y - b_y, 2 ) );
    bc_len = sqrt ( pow ( b_x - c_x, 2 ) + pow ( b_y - c_y, 2 ) );
    ca_len = sqrt ( pow ( c_x - a_x, 2 ) + pow ( c_y - a_y, 2 ) );
//
//  Take care of a ridiculous special case.
//
    if ( ab_len == 0.0 && bc_len == 0.0 && ca_len == 0.0 )
    {
      a_angle = 2.0 * pi / 3.0;
      b_angle = 2.0 * pi / 3.0;
      c_angle = 2.0 * pi / 3.0;
    }
    else
    {
      if ( ca_len == 0.0 || ab_len == 0.0 )
      {
        a_angle = pi;
      }
      else
      {
        a_angle = arc_cosine ( 
          ( ca_len * ca_len + ab_len * ab_len - bc_len * bc_len ) 
          / ( 2.0 * ca_len * ab_len ) );
      }

      if ( ab_len == 0.0 || bc_len == 0.0 )
      {
        b_angle = pi;
      }
      else
      {
        b_angle = arc_cosine ( 
          ( ab_len * ab_len + bc_len * bc_len - ca_len * ca_len ) 
          / ( 2.0 * ab_len * bc_len ) );
      }

      if ( bc_len == 0.0 || ca_len == 0.0 )
      {
        c_angle = pi;
      }
      else
      {
        c_angle = arc_cosine ( 
          ( bc_len * bc_len + ca_len * ca_len - ab_len * ab_len ) 
          / ( 2.0 * bc_len * ca_len ) );
      }
    }
    *alpha_min = r8_min ( *alpha_min, a_angle );
    *alpha_min = r8_min ( *alpha_min, b_angle );
    *alpha_min = r8_min ( *alpha_min, c_angle );

    *alpha_ave = *alpha_ave + *alpha_min;

    *alpha_area = *alpha_area + area * *alpha_min;

    area_total = area_total + area;
  }
  *alpha_ave = *alpha_ave / ( double ) ( triangle_num );
  *alpha_area = *alpha_area / area_total;
//
//  Normalize angles from [0,pi/3] radians into qualities in [0,1].
//
  *alpha_min = *alpha_min * 3.0 / pi;
  *alpha_ave = *alpha_ave * 3.0 / pi;
  *alpha_area = *alpha_area * 3.0 / pi;

  return;
}
//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
# define PI 3.141592653589793

  if ( c <= -1.0 )
  {
    return PI;
  } 
  else if ( 1.0 <= c )
  {
    return 0.0;
  }
  else
  {
    return acos ( c );
  }
# undef PI
}
//****************************************************************************80

void area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double &area_min, double &area_max, double &area_ratio,
  double &area_ave, double &area_std, int &area_negative, int &area_zero )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_MEASURE determines the area ratio quality measure.
//
//  Discusion:
//
//    This measure computes the area of every triangle, and returns
//    the ratio of the minimum to the maximum triangle.  A value of
//    1 is "perfect", indicating that all triangles have the same area.
//    A value of 0 is the worst possible result.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double &AREA_MIN, &AREA_MAX, the minimum and maximum 
//    areas.
//
//    Output, double &AREA_RATIO, the ratio of the minimum to the
//    maximum area.
//
//    Output, double &AREA_AVE, the average area.
//
//    Output, double &AREA_STD, the standard deviation of the areas.
//
//    Output, int &AREA_NEGATIVE, the number of triangles with negative area.
//
//    Output, int &AREA_ZERO, the number of triangles with zero area.
{
  double area;
  int triangle;
  double value;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  area_max = 0.0;
  area_min = r8_huge ( );
  area_ave = 0.0;
  area_negative = 0;
  area_zero = 0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 ) 
                        + x2 * ( y3 - y1 ) 
                        + x3 * ( y1 - y2 ) );

    area_min = r8_min ( area_min, r8_abs ( area ) );
    area_max = r8_max ( area_max, r8_abs ( area ) );
    area_ave = area_ave + r8_abs (area );

    if ( area < 0.0 )
    {
      area_negative = area_negative + 1;
    }
    if ( area == 0.0 )
    {
      area_zero = area_zero + 1;
    }
  }

  area_ave = area_ave / ( double ) ( triangle_num );
  area_std = 0.0;
  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 ) 
                        + x2 * ( y3 - y1 ) 
                        + x3 * ( y1 - y2 ) );

    area_std = area_std + pow ( area - area_ave, 2 );
  }
  area_std = sqrt ( area_std / ( double ) ( triangle_num ) );

  if ( 0.0 < area_max )
  {
    area_ratio = area_min / area_max;
  }
  else
  {
    area_ratio = 0.0;
  }
  return;
}
//****************************************************************************80

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
//
//  Discussion:
//
//    The quantity computed here is the "geometric" bandwidth determined
//    by the finite element mesh alone.
//
//    If a single finite element variable is associated with each node
//    of the mesh, and if the nodes and variables are numbered in the
//    same way, then the geometric bandwidth is the same as the bandwidth
//    of a typical finite element matrix.
//
//    The bandwidth M is defined in terms of the lower and upper bandwidths:
//
//      M = ML + 1 + MU
//
//    where 
//
//      ML = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but earlier column,
//
//      MU = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but later column.
//
//    Because the finite element node adjacency relationship is symmetric,
//    we are guaranteed that ML = MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
//
//    Output, int *M, the bandwidth of the matrix.
//
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
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

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 ) 
  {
    c1 = c1 - 32;
  } 
  if ( 97 <= c2 && c2 <= 122 ) 
  {
    c2 = c2 - 32;
  }     

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

int file_column_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file are presumed to consist of COLUMN_NUM words, separated
//    by spaces.  There may also be some blank lines, and some comment lines,
//    which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char text[255];
//
//  Open the file.
//
  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( text ) == 0 )
    {
      continue;
    }

    if ( text[0] == '#' )
    {
      continue;
    }
    got_one = true;
    break;
  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( filename.c_str ( ) );

    for ( ; ; )
    {
      input.getline ( text, sizeof ( text ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( text ) == 0 )
      {
        continue;
      }
      got_one = true;
      break;
    }
  }

  input.close ( );

  if ( !got_one )
  {
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Warning!\n";
    cerr << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( text );

  return column_num;
}
//****************************************************************************80

int file_row_count ( string filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 January 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  int record_num;
  int row_num;
  char text[255];

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the file: \"" << filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input.getline ( text, sizeof ( text ) );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( text[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( text ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
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
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
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

int *i4mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_DATA_READ reads data from an I4MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, int I4MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  int *table;
  int *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "I4MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new int[m*n];

  x = new int[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_i4vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80
 
void i4mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_HEADER_READ reads the header from an I4MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );
 
  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }
 
  *n = file_row_count ( input_filename );
 
  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "I4MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }
 
  return;
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Discussion:
//
//    An I4MAT is an MxN array of I4's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title for the matrix.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

void mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ZERO ensures that the element definition is zero-based.
//
//  Discussion:
//
//    The ELEMENT_NODE array contains nodes indices that form elements.
//    The convention for node indexing might start at 0 or at 1.
//    Since a C++ program will naturally assume a 0-based indexing, it is
//    necessary to check a given element definition and, if it is actually
//    1-based, to convert it.
//
//    This function attempts to detect 1-based node indexing and correct it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input/output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the element
//    definitions.
//
{
  int element;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = node_num + 1;
  node_max = -1;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      node_min = i4_min ( node_min, node );
      node_max = i4_max ( node_max, node );
    }
  }

  if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] - 1;
      }
    }
  }
  else if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

void q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[], double *q_min, double *q_max, double *q_ave, 
  double *q_area )

//****************************************************************************80
//
//  Purpose:
//
//    Q_MEASURE determines the triangulated pointset quality measure Q.
//
//  Discussion:
//
//    The Q measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    For a single triangle T, the value of Q(T) is defined as follows:
//
//      TAU_IN = radius of the inscribed circle,
//      TAU_OUT = radius of the circumscribed circle,
//
//      Q(T) = 2 * TAU_IN / TAU_OUT
//        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
//
//    where A, B and C are the lengths of the sides of the triangle T.
//
//    The Q measure computes the value of Q(T) for every triangle T in the
//    triangulation, and then computes the minimum of this
//    set of values:
//
//      Q_MEASURE = min ( all T in triangulation ) Q(T)
//
//    In an ideally regular mesh, all triangles would have the same
//    equilateral shape, for which Q = 1.  A good mesh would have
//    0.5 < Q.
//
//    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
//    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
//    triangles.  Generally, a maximal triangulation is expected, namely,
//    a triangulation whose image is a planar graph, but for which the
//    addition of any new triangle would mean the graph was no longer planar.
//    A Delaunay triangulation is a maximal triangulation which maximizes
//    the minimum angle that occurs in any triangle.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//    Per-Olof Persson and Gilbert Strang,
//    A Simple Mesh Generator in MATLAB,
//    SIAM Review,
//    Volume 46, Number 2, pages 329-345, June 2004.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double *Q_MIN, *Q_MAX, the minimum and maximum values
//    of Q over all triangles.
//
//    Output, double *Q_AVE, the average value of Q.
//
//    Output, double *Q_AREA, the average value of Q, weighted by
//    the area of each triangle.
//
{
  int a_index;
  double ab_length;
  double area;
  double area_total;
  int b_index;
  double bc_length;
  int c_index;
  double ca_length;
  double q;
  int triangle;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  *q_min = r8_huge ( );
  *q_max = - r8_huge ( );
  *q_ave = 0.0;
  *q_area = 0.0;
  area_total = 0.0;

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*triangle_order];
    b_index = triangle_node[1+triangle*triangle_order];
    c_index = triangle_node[2+triangle*triangle_order];

    ab_length = sqrt ( 
        pow ( z[0+(a_index-1)*2] - z[0+(b_index-1)*2], 2 ) 
      + pow ( z[1+(a_index-1)*2] - z[1+(b_index-1)*2], 2 ) );

    bc_length = sqrt ( 
        pow ( z[0+(b_index-1)*2] - z[0+(c_index-1)*2], 2 ) 
      + pow ( z[1+(b_index-1)*2] - z[1+(c_index-1)*2], 2 ) );

    ca_length = sqrt ( 
        pow ( z[0+(c_index-1)*2] - z[0+(a_index-1)*2], 2 ) 
      + pow ( z[1+(c_index-1)*2] - z[1+(a_index-1)*2], 2 ) );

    q = ( bc_length + ca_length - ab_length ) 
      * ( ca_length + ab_length - bc_length ) 
      * ( ab_length + bc_length - ca_length ) 
      / ( ab_length * bc_length * ca_length );

    x1 = z[0+(triangle_node[0+triangle*triangle_order]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*triangle_order]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*triangle_order]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*triangle_order]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*triangle_order]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*triangle_order]-1)*2];

    area = 0.5 * r8_abs ( x1 * ( y2 - y3 ) 
                        + x2 * ( y3 - y1 ) 
                        + x3 * ( y1 - y2 ) );

    *q_min = min ( *q_min, q );
    *q_max = max ( *q_max, q );
    *q_ave = *q_ave + q;
    *q_area = *q_area + q * area;

    area_total = area_total + area;
  }

  *q_ave = *q_ave / ( double ) ( triangle_num );

  if ( 0.0 < area_total )
  {
    *q_area = *q_area / area_total;
  }
  else
  {
    *q_area = 0.0;
  }
  return;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
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

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double *r8mat_data_read ( string input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DATA_READ reads the data from an R8MAT file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with '#' are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double R8MAT_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  string line;
  double *table;
  double *x;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "R8MAT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    getline ( input, line );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80
 
void r8mat_header_read ( string input_filename, int *m, int *n )
 
//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_HEADER_READ reads the header from an R8MAT file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points.
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_HEADER_READ - Fatal error!\n";
    cerr << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, an optional title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n ) 
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( string s, int *last, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int *LAST, the last character of S used to make IVAL.
//
//    Output, bool *ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
{
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  *error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for ( ; ; ) 
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if ( istate == 0 )
    {
      if ( c == ' ' )
      {
      }
      else if ( c == '-' )
      {
        istate = 1;
        isgn = -1;
      }
      else if ( c == '+' )
      {
        istate = 1;
        isgn = + 1;
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if ( istate == 1 )
    {
      if ( c == ' ' )
      {
      }
      else if ( '0' <= c && c <= '9' )
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        *error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if ( istate == 2 )
    {
      if ( '0' <= c && c <= '9' )
      {
        ival = 10 * (ival) + c - '0';
      }
      else
      {
        ival = isgn * ival;
        *last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if ( istate == 2 )
  {
    ival = isgn * ival;
    *last = s_len_trim ( s );
  }
  else
  {
    *error = true;
    *last = 0;
  }

  return ival;
}
//****************************************************************************80

bool s_to_i4vec ( string s, int n, int ivec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4VEC reads an I4VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, int IVEC[N], the values read from the string.
//
//    Output, bool S_TO_I4VEC, is TRUE if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    ivec[i] = s_to_i4 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( string s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  int begin;
  bool error;
  int i;
  int lchar;
  int length;

  begin = 0;
  length = s.length ( );
  error = 0;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s.substr(begin,length), &lchar, &error );

    if ( error )
    {
      return error;
    }
    begin = begin + lchar;
    length = length - lchar;
  }

  return error;
}
//****************************************************************************80

int s_word_count ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int char_count;
  int i;
  int word_count;

  word_count = 0;
  blank = true;

  char_count = s.length ( );

  for ( i = 0; i < char_count; i++ )
  {
    if ( isspace ( s[i] ) )
    {
      blank = true;
    }
    else if ( blank )
    {
      word_count = word_count + 1;
      blank = false;
    }
  }

  return word_count;
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
