# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void mesh_base_zero ( int node_num, int element_order, int element_num, 
  int element_node[] );
void timestamp ( );
void triangle_element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] );
void triangle_element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num );
void triangle_node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] );
void triangle_node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num );
void xml_mesh2d_write ( string xml_filename, int m, int node_num, double node_x[],
  int element_order, int element_num, int element_node[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_TO_XML.
//
//  Discussion:
//
//    The TRIANGLE program creates "node" and "element" files that define
//    a triangular mesh.  A typical pair of such files might have the names
//    "suv.node" and "suv.ele".
//
//    This program reads this pair of files and creates an equivalent XML file
//    whose name might be "suv.xml".
//
//  Usage:
//
//    triangle_to_xml prefix
//
//    where 'prefix' is the common filename prefix so that:
//
//    * prefix.node contains the coordinates of the nodes;
//    * prefix.ele contains the indices of nodes forming each element.
//    * prefix.xml will be the XML file created by the program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  int element_att_num;
  double *element_att;
  int *element_node;
  int element_num;
  int element_order;
  int i;
  int j;
  int m;
  double *node_att;
  int node_att_num;
  int node_marker_num;
  int *node_marker;
  double *node_x;
  int node_num;
  string prefix;
  string triangle_element_filename;
  string triangle_node_filename;
  string xml_filename;

  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_TO_XML:\n";
  cout << "  C++ version\n";
  cout << "  Read a pair of NODE and ELE files created by TRIANGLE.\n";
  cout << "  Write a corresponding XML mesh file.\n";
//
//  Get the filename prefix.
//
  if ( 1 <= argc )
  {
    prefix = argv[1];
  }
  else
  {
    cout << "\n";
    cout << "  Please enter the filename prefix:\n";

    cin >> prefix;
  }
//
//  Create the file names.
//
  triangle_node_filename = prefix + ".node";
  triangle_element_filename = prefix + ".ele";
  xml_filename = prefix + ".xml";

  cout << "\n";
  cout << "  Read:\n";
  cout << "  * TRIANGLE node file \"" << triangle_node_filename << "\"\n";
  cout << "  * TRIANGLE element file \"" << triangle_element_filename << "\".\n";
  cout << "  Create:\n";
  cout << "  * XML file \"" << xml_filename << "\".\n";
//
//  Read the TRIANGLE NODE data.
//
  triangle_node_size_read ( triangle_node_filename, &node_num, &m, &node_att_num, 
    &node_marker_num );

  node_x = new double[m*node_num];
  node_att = new double[node_att_num*node_num];
  node_marker = new int[node_marker_num*node_num];

  triangle_node_data_read ( triangle_node_filename, node_num, m, node_att_num, 
    node_marker_num, node_x, node_att, node_marker );
//
//  Read the TRIANGLE ELE data.
//
  triangle_element_size_read ( triangle_element_filename, &element_num, &element_order, 
    &element_att_num );

  element_node = new int[element_order*element_num];
  element_att = new double[element_att_num*element_num];

  triangle_element_data_read ( triangle_element_filename, element_num, element_order, 
    element_att_num, element_node, element_att );
//
//  Write the XML file.
//
  xml_mesh2d_write ( xml_filename, m, node_num, node_x, element_order,
    element_num, element_node );

  cout << "\n";
  cout << "  Created XML file \"" << xml_filename << "\".\n";
//
//  Free memory.
//
  delete [] element_att;
  delete [] element_node;
  delete [] node_att;
  delete [] node_x;
  delete [] node_marker;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_TO_XML:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 October 2014
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
  int i;
  const int i4_huge = 2147483647;
  int j;
  int node_max;
  int node_min;
  int t;

  node_min = + i4_huge;
  node_max = - i4_huge;

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < element_order; i++ )
    {
      t = element_node[i+j*element_order];
      if ( t < node_min )
      {
        node_min = t;
      }
      if ( node_max < t )
      {
        node_max = t;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  No conversion is necessary.\n";

  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ZERO:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  This will be converted to 0-based.\n";
    for ( j = 0; j < element_num; j++ )
    {
      for ( i = 0; i < element_order; i++ )
      {
        element_node[i+j*element_order] = element_node[i+j*element_order] - 1;
      }
    }
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

//****************************************************************************80

void triangle_element_data_read ( string element_file, int element_num, 
  int element_order, int element_att_num, int element_node[], double element_att[] )

//*****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ELEMENT_DATA_READ reads the header information from an element file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_ATT_NUM, number of element attributes listed on each 
//    node record.
//
//    Output, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the indices of the
//    nodes that make up each element.
//
//    Output, double ELEMENT_ATT[ELEMENT_ATT_NUM*ELEMENT_NUM], the attributes
//    of each element.
//
{
  int element;
  int i;
  int i1;
  int i2;
  int i3;
  ifstream input;
  int ival;
  double value;

  element = -1;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "TRIANGLE_ELEMENT_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
/*
  Read, but ignore, dimension line.
*/
    if ( element == -1 )
    {
      input >> i1 >> i2 >> i3;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < element_order; i++ )
      {
        input >> ival;
        element_node[i+element*element_order] = ival;
      }
      for ( i = 0; i < element_att_num; i++ )
      {
        input >> value;
        element_att[i+element*element_att_num] = value;
      }
    }

    element = element + 1;

    if ( element_num <= element )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void triangle_element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ELEMENT_SIZE_READ reads the header information from an element file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string ELEMENT_FILE, the name of the file to be read.
//
//    Output, int *ELEMENT_NUM, the number of elements.
//
//    Output, int *ELEMENT_ORDER, the order of the elements.
//
//    Output, int *ELEMENT_ATT_NUM, the number of element attributes.
//
{
  ifstream input;

  input.open ( element_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "TRIANGLE_ELEMENT_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  input >> *element_num >> *element_order >> *element_att_num;

  input.close ( );

  return;
}
//****************************************************************************80

void triangle_node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NODE_HEADER_READ reads the header information from a node file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Input, int NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
//    Output, double NODE_COORD[NODE_DIM*NODE_NUM], the nodal coordinates.
//
//    Output, double NODE_ATT[NODE_ATT_NUM*NODE_NUM], the nodal attributes.
//
//    Output, int NODE_MARKER[NODE_MARKER_NUM*NODE_NUM], the node markers.
//
{
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  ifstream input;
  int ival;
  int node;
  double value;

  node = -1;

  input.open ( node_file.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "TRIANGLE_NODE_DATA_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
//
//  Read, but ignore, dimension line.
//
    if ( node == -1 )
    {
      input >> i1 >> i2 >> i3 >> i4;
    }
    else
    {
      input >> ival;

      for ( i = 0; i < node_dim; i++ )
      {
        input >> value;
        node_coord[i+node*node_dim] = value;
      }
      for ( i = 0; i < node_att_num; i++ )
      {
        input >> value;
        node_att[i+node*node_att_num] = value;
      }
      for ( i = 0; i < node_marker_num; i++ )
      {
        input >> ival;
        node_marker[i+node*node_marker_num] = ival;
      }
    }

    node = node + 1;

    if ( node_num <= node )
    {
      break;
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void triangle_node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_NODE_SIZE_READ reads the header information from a node file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string NODE_FILE, the name of the node file to be read.
//
//    Output, int *NODE_NUM, the number of nodes.
//
//    Output, int *NODE_DIM, the spatial dimension.
//
//    Output, int *NODE_ATT_NUM, number of node attributes listed on each 
//    node record.
//
//    Output, int *NODE_MARKER_NUM, 1 if every node record includes a final
//    boundary marker value.
//
{
  ifstream input;

  input.open ( node_file.c_str ( ) );

  input >> *node_num
        >> *node_dim
        >> *node_att_num
        >> *node_marker_num;

  input.close ( );

  return;
}
//****************************************************************************80

void xml_mesh2d_write ( string xml_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    XML_MESH2D_WRITE writes a 2D mesh as a DOLFIN XML file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Anders Logg, Kent-Andre Mardal, Garth Wells,
//    Automated Solution of Differential Equations by the Finite Element
//    Method: The FEniCS Book,
//    Lecture Notes in Computational Science and Engineering,
//    Springer, 2011,
//    ISBN13: 978-364223098
//
//  Parameters:
//
//    Input, string XML_FILENAME, the name of the XML file 
//    to create.
//
//    Input, int M, the spatial dimension.
//
//    Input, inte NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[M*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int element;
  int element_type;
  ofstream xml;
  int i;
  int node;
//
//  Force 0-based indexing.
//
  mesh_base_zero ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  xml.open ( xml_filename.c_str ( ) );
//
//  Write the data.
//
  xml << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  xml << "\n";
  xml << "<dolfin xmlns:dolfin=\"http://www.fenics.org/dolfin/\">\n";
  xml << "  <mesh celltype=\"triangle\" dim=\"" << m << "\">\n";

  xml << "    <vertices size=\"" << node_num << "\">\n";
  for ( node = 0; node < node_num; node++ )
  {
    xml << "      <vertex index =\"" << node
        << "\" x =\"" << node_x[0+node*m] 
        << "\" y =\"" << node_x[1+node*m] << "/>\n";
  }
  xml << "    </vertices>\n";

  xml << "    <cells size=\"" << element_num << "\">\n";
  for ( element = 0; element < element_num; element++ )
  {
    xml << "      <triangle index =\"" << element
        << "\" v0 =\"" << element_node[0+element*element_order]
        << "\" v1 =\"" << element_node[1+element*element_order] 
        << "\" v2 =\"" << element_node[2+element*element_order] << "\"/>\n";
  }
  xml << "    </cells>\n";
  xml << "  </mesh>\n";
  xml << "</dolfin>\n";

  xml.close ( );

  return;
}
