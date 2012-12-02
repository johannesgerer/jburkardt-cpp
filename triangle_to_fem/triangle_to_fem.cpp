# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] );
void element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num );
void i4mat_write ( string output_filename, int m, int n, int table[] );
void node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] );
void node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num );
void r8mat_write ( string output_filename, int m, int n, double table[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_TO_FEM.
//
//  Discussion:
//
//    The TRIANGLE program creates "node" and "element" files that define
//    a triangular mesh.  A typical pair of such files might have the names
//    "suv.node" and "suv.ele".
//
//    This program reads this pair of files and creates a pair of FEM files, whose
//    names might be "suv_nodes.txt" and "suv_elements.txt".
//
//  Usage:
//
//    triangle_to_fem prefix
//
//    where 'prefix' is the common filename prefix so that:
//
//    * prefix.node contains the coordinates of the nodes;
//    * prefix.ele contains the indices of nodes forming each element.
//    * prefix_nodes.txt will be the FEM node file created by the program;
//    * prefix_elements.txt will be the FEM element file created by the program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  int dim_num;
  string element_filename;
  int element_att_num;
  double *element_att;
  int *element_node;
  int element_num;
  int element_order;
  string fem_element_filename;
  string fem_node_filename;
  int i;
  int j;
  double *node_att;
  int node_att_num;
  int node_marker_num;
  string node_filename;
  int *node_marker;
  double *node_coord;
  int node_dim;
  int node_num;
  string prefix;

  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_TO_FEM:\n";
  cout << "  C++ version\n";
  cout << "  Read a pair of NODE and ELE files created by TRIANGLE.\n";
  cout << "  Write a corresponding pair of FEM node and element files.\n";
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
  node_filename = prefix + ".node";
  element_filename = prefix + ".ele";
  fem_node_filename = prefix + "_nodes.txt";
  fem_element_filename = prefix + "_elements.txt";

  cout << "\n";
  cout << "  Read Node file \"" << node_filename << "\"\n";
  cout << "    and Element file \"" << element_filename << "\".\n";
  cout << "  Create FEM node file \"" << fem_node_filename << "\"\n";
  cout << "    and FEM element file \"" << fem_element_filename << "\".\n";
//
//  Read the TRIANGLE NODE data.
//
  node_size_read ( node_filename, &node_num, &node_dim, &node_att_num, 
    &node_marker_num );

  node_coord = new double[2*node_num];
  node_att = new double[node_att_num*node_num];
  node_marker = new int[node_num];

  node_data_read ( node_filename, node_num, node_dim, node_att_num, 
    node_marker_num, node_coord, node_att, node_marker );
//
//  Read the TRIANGLE ELE data.
//
  element_size_read ( element_filename, &element_num, &element_order, 
    &element_att_num );

  element_node = new int[element_order*element_num];
  element_att = new double[element_att_num*element_num];

  element_data_read ( element_filename, element_num, element_order, 
    element_att_num, element_node, element_att );
//
//  Write the FEM NODE and ELEMENT data.
//
  dim_num = 2;

  r8mat_write ( fem_node_filename, dim_num, node_num, node_coord );

  i4mat_write ( fem_element_filename, element_order, element_num, element_node );
//
//  Free memory.
//
  delete [] element_att;
  delete [] element_node;
  delete [] node_att;
  delete [] node_coord;
  delete [] node_marker;
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_TO_FEM:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );
}
//****************************************************************************80

void element_data_read ( string element_file, int element_num, int element_order, 
  int element_att_num, int element_node[], double element_att[] )

//*****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_DATA_READ reads the header information from an element file.
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
    cerr << "ELEMENT_DATA_READ - Fatal error!\n";
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

void element_size_read ( string element_file, int *element_num, 
  int *element_order, int *element_att_num )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_SIZE_READ reads the header information from an element file.
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
    cerr << "ELEMENT_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open file.\n";
    exit ( 1 );
  }

  input >> *element_num >> *element_order >> *element_att_num;

  input.close ( );

  return;
}
//****************************************************************************80

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
//
//  Discussion:
//
//    An I4MAT is an array of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "I4MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
//****************************************************************************80

void node_data_read ( string node_file, int node_num, int node_dim, 
  int node_att_num, int node_marker_num, double node_coord[], double node_att[],
  int node_marker[] )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_HEADER_READ reads the header information from a node file.
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
    cerr << "NODE_DATA_READ - Fatal error!\n";
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

void node_size_read ( string node_file, int *node_num, int *node_dim, 
  int *node_att_num, int *node_marker_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_SIZE_READ reads the header information from a node file.
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

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
