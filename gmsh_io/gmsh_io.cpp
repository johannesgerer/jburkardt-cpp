# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

# include "gmsh_io.hpp"

//****************************************************************************80

char ch_cap ( char ch )

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
//    Input, char CH, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= ch && ch <= 122 ) 
  {
    ch = ch - 32;
  }   

  return ch;
}
//****************************************************************************80

bool ch_eqi ( char ch1, char ch2 )

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
//    Input, char CH1, CH2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= ch1 && ch1 <= 122 ) 
  {
    ch1 = ch1 - 32;
  } 
  if ( 97 <= ch2 && ch2 <= 122 ) 
  {
    ch2 = ch2 - 32;
  }     

  return ( ch1 == ch2 );
}
//****************************************************************************80

int ch_to_digit ( char ch )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     CH  DIGIT
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
//    Input, char CH, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the 
//    character was 'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= ch && ch <= '9' )
  {
    digit = ch - '0';
  }
  else if ( ch == ' ' )
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

void gmsh_data_read ( string gmsh_filename, int node_dim, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_DATA_READ reads data from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, character *GMSH_FILENAME, the GMSH filename.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_X[NODE_DIM*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
{
  int i;
  int i4_dummy;
  bool ierror;
  int indx;
  ifstream input;
  int j;
  int k;
  int length;
  int level;
  string text;
  double x;
  double y;
  double z;

  input.open ( gmsh_filename.c_str ( ) );

  if ( ! input )
  {
    cerr << "\n";
    cerr << "GMSH_DATA_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      i4_dummy = s_to_i4 ( text, length, ierror );
      level = 2;
      j = 0;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        i4_dummy = s_to_i4 ( text, length, ierror );
        text.erase ( 0, length );
        for ( i = 0; i < node_dim; i++ )
        {
          x = s_to_r8 ( text, length, ierror );
          text.erase ( 0, length );
          node_x[i+j*node_dim] = x;
        }
        j = j + 1;
      }
    }
  }
//
//  Now read element information.
//
  level = 0;
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      i4_dummy = s_to_i4 ( text, length, ierror );
      level = 2;
      j = 0;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndElements" ) )
      {
        break;
      }
      else
      {
        for ( k = 1; k <= 5; k++ )
        {
          i4_dummy = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
        }
        for ( i = 0; i < element_order; i++ )
        {
          k = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
          element_node[i+j*element_order] = k;
        }
        j = j + 1;
      }
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

void gmsh_size_read ( string gmsh_filename, int &node_num, int &node_dim,
  int &element_num, int &element_order )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_SIZE_READ reads sizes from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the GMSH filename.
//
//    Output, int &NODE_NUM, the number of nodes.
//
//    Output, int &NODE_DIM, the spatial dimension.
//
//    Output, int &ELEMENT_NUM, the number of elements.
//
//    Output, int &ELEMENT_ORDER, the order of the elements.
//
{
  bool ierror;
  int indx;
  ifstream input;
  int k;
  int length;
  int level;
  const double r8_big = 1.0E+30;
  string text;
  double x;
  double x_max;
  double x_min;
  double y;
  double y_max;
  double y_min;
  double z;
  double z_max;
  double z_min;

  node_num = 0;
  node_dim = 0;

  x_max = - r8_big;
  x_min = + r8_big;
  y_max = - r8_big;
  y_min = + r8_big;
  z_max = - r8_big;
  z_min = + r8_big;

  input.open ( gmsh_filename.c_str ( ) );

  if ( ! input )
  {
    cerr << "\n";
    cerr << "GMSH_SIZE_READ - Fatal error!\n";
    cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
    exit ( 1 );
  }

  level = 0;
 
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Nodes" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      node_num = s_to_i4 ( text, length, ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndNodes" ) )
      {
        break;
      }
      else
      {
        indx = s_to_i4 ( text, length, ierror );
        text.erase ( 0, length );

        x = s_to_r8 ( text, length, ierror );
        x_min = r8_min ( x_min, x );
        x_max = r8_max ( x_max, x );
        text.erase ( 0, length );

        y = s_to_r8 ( text, length, ierror );
        y_min = r8_min ( y_min, y );
        y_max = r8_max ( y_max, y );
        text.erase ( 0, length );

        z = s_to_r8 ( text, length, ierror);
        z_min = r8_min ( z_min, z );
        z_max = r8_max ( z_max, z );
        text.erase ( 0, length );
      }
    }
  }
//
//  Make a very simple guess as to the dimensionality of the data.
//
  node_dim = 3;
  if ( z_max == z_min )
  {
    node_dim = 2;
    if ( y_max == y_min )
    {
      node_dim = 1;
    }
  }
//
//  Now read element information.
//
  level = 0;

  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( level == 0 )
    {
      if ( s_begin ( text, "$Elements" ) )
      {
        level = 1;
      }
    }
    else if ( level == 1 )
    {
      element_num = s_to_i4 ( text, length, ierror );
      level = 2;
    }
    else if ( level == 2 )
    {
      if ( s_begin ( text, "$EndElements" ) )
      {
        break;
      }
      else
      {
        k = 0;
        for ( ; ; )
        {
          indx = s_to_i4 ( text, length, ierror );
          text.erase ( 0, length );
          if ( ierror != 0 )
          {
            break;
          }
          k = k + 1;
        }
        element_order = k - 5;
        break;
      }
    }
  }

  input.close ( );

  return;
}
//****************************************************************************80

int *gmsh_mesh2d_element_data_example ( int element_num, int element_order )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_ELEMENT_DATA_EXAMPLE returns element information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Output, int GMSH_MESH2D_ELEMENT_DATA_EXAMPLE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the indices of the nodes that make up each element.
//
{
  int *element_node;
  int element_node_save[3*24] = {
    1,  2,  6,
    7,  6,  2,
    2,  3,  7,
    8,  7,  3,
    3,  4,  8,
    9,  8,  4,
    4,  5,  9,
   10,  9,  5,
    6,  7, 11,
   12, 11,  7,
    7,  8, 12,
   13, 12,  8,
    8,  9, 13,
   14, 13,  9,
    9, 10, 14,
   15, 14, 10,
   11, 12, 16,
   17, 16, 12,
   12, 13, 17,
   18, 17, 13,
   16, 17, 19,
   20, 19, 17,
   17, 18, 20,
   21, 20, 18 };

  element_node = i4mat_copy_new ( element_order, element_num, 
    element_node_save );

  return element_node;
}
//****************************************************************************80

void gmsh_mesh2d_element_size_example ( int &element_num, int &element_order )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_ELEMENT_SIZE_EXAMPLE returns element sizes for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int &ELEMENT_NUM, the number of elements.
//
//    Output, int &ELEMENT_ORDER, the order of the elements.
//
{
  element_num = 24;
  element_order = 3;

  return;
}
//****************************************************************************80

double *gmsh_mesh2d_node_data_example ( int node_num, int node_dim )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_NODE_DATA_EXAMPLE returns node information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Output, double GMSH_MESH2D_NODE_DATA_EXAMPLE[NODE_DIM*NODE_NUM], 
//    the nodal coordinates.
//
{
  double *node_coord;
  double node_coord_save[2*21] = {
    0.0, 0.0,
    1.0, 0.0,
    2.0, 0.0,
    3.0, 0.0,
    4.0, 0.0,
    0.0, 1.0,
    1.0, 1.0,
    2.0, 1.0,
    3.0, 1.0,
    4.0, 1.0,
    0.0, 2.0,
    1.0, 2.0,
    2.0, 2.0,
    3.0, 2.0,
    4.0, 2.0,
    0.0, 3.0,
    1.0, 3.0,
    2.0, 3.0,
    0.0, 4.0,
    1.0, 4.0,
    2.0, 4.0 };

  node_coord = r8mat_copy_new ( 2, 21, node_coord_save );

  return node_coord;
}
//****************************************************************************80

void gmsh_mesh2d_node_size_example ( int &node_num, int &node_dim )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_NODE_SIZE_EXAMPLE returns the sizes of node information for the example.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int &NODE_NUM, the number of nodes.
//
//    Output, int &NODE_DIM, the spatial dimension.
//
{
  node_num = 21;
  node_dim = 2;

  return;
}
//****************************************************************************80

void gmsh_mesh1d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH1D_WRITE writes 1d mesh data as a Gmsh mesh file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
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
  ofstream gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  0.0  0.0\n";
  }
  gmsh << "$EndNodes\n";

  element_type = 1;

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    for ( i = 0; i < element_order; i++ )
    {
      gmsh << "  " << element_node[i+element*element_order];
    }
    gmsh << "\n";
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

void gmsh_mesh2d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH2D_WRITE writes 2d mesh data as a Gmsh mesh file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
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
  ofstream gmsh;
  int i;
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  " << node_x[1+node*m]
         << "  0.0\n";
  }
  gmsh << "$EndNodes\n";

  if ( element_order == 3 )
  {
    element_type = 2;
  }
  else if ( element_order == 6 )
  {
    element_type = 9;
  }

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    for ( i = 0; i < element_order; i++ )
    {
      gmsh << "  " << element_node[i+element*element_order];
    }
    gmsh << "\n";
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

void gmsh_mesh3d_write ( string gmsh_filename, int m, int node_num, 
  double node_x[], int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_MESH3D_WRITE writes 3D mesh data as a Gmsh mesh file.
//
//  Discussion:
//
//    The node ordering for the 20 node element is not standard.
//
//    Assuming the vertices are A, B, C and D, Gmsh uses the following ordering:
//
//    1:    a
//    2:        b
//    3:            c
//    4:                d
//    5: (2*a  +b        )/3
//    6: (  a+2*b        )/3
//    7: (    2*b+  c    )/3
//    8: (      b+2*c    )/3
//    9: (  a    +2*c    )/3
//   10: (2*a    +  c    )/3
//   11: (2*a        +  d)/3
//   12: (  a        +2*d)/3
//   13: (     b     +2*d)/3
//   14: (   2*b     +  d)/3
//   15: (       +  c+2*d)/3
//   16: (       +2*c+  d)/3
//   17: (  a+  b+  c    )/3
//   18: (  a+  b    +  d)/3
//   19: (      b+  c+  d)/3
//   20: (  a+      c+  d)/3
//
//    Leo Rebholz used the following ordering:
//
//    1:    a
//    2:        b
//    3:            c
//    4:                d
//    5: (2*a  +b        )/3
//    6: (2*a    +  c    )/3
//    7: (  a+2*b        )/3
//    8: (  a    +2*c    )/3
//    9: (  a+  b+  c    )/3
//   10: (    2*b+  c    )/3
//   11: (      b+2*c    )/3
//   12: (2*a        +  d)/3
//   13: (   2*b     +  d)/3
//   14: (       +2*c+  d)/3
//   15: (  a+  b    +  d)/3
//   16: (      b+  c+  d)/3
//   17: (  a+      c+  d)/3
//   18: (  a        +2*d)/3
//   19: (     b     +2*d)/3
//   20: (       +  c+2*d)/3
//
//    Since the only 20 node data we have is from Leo, we will assume that
//    all 20 node input data is in Leo's format, and needs to be converted
//    to the Gmsh convention.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Christophe Geuzaine, Jean-Francois Remacle,
//    Gmsh: a three-dimensional finite element mesh generator with
//    built-in pre- and post-processing facilities,
//    International Journal for Numerical Methods in Engineering,
//    Volume 79, Number 11, pages 1309-1331, 2009.
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the name of the Gmsh file.
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
  ofstream gmsh;
  int i;
  int i2;
  int leo_to_gmsh[20] = {
     0,  1,  2,  3,  4, 
     6,  9, 10,  7,  5, 
    11, 17, 18, 12, 19, 
    13,  8, 14, 15, 16 };
  int node;
  int tag_num;
  int tag1;
//
//  Detect and correct 0-based node indexing.
//
  mesh_base_one ( node_num, element_order, element_num, element_node );
//
//  Open the file.
//
  gmsh.open ( gmsh_filename.c_str ( ) );
//
//  Write the data.
//
  gmsh << "$MeshFormat\n";
  gmsh << "2.2 0 8\n";
  gmsh << "$EndMeshFormat\n";

  gmsh << "$Nodes\n";
  gmsh << node_num << "\n";
  for ( node = 0; node < node_num; node++ )
  {
    gmsh << node + 1 
         << "  " << node_x[0+node*m]
         << "  " << node_x[1+node*m]
         << "  " << node_x[2+node*m] << "\n";
  }
  gmsh << "$EndNodes\n";

  if ( element_order == 4 )
  {
    element_type = 4;
  }
  else if ( element_order == 10 )
  {
    element_type = 11;
  }
  else if ( element_order == 20 )
  {
    element_type = 29;
  }

  tag_num = 2;
  tag1 = 0;
  gmsh << "$Elements\n";
  gmsh << element_num << "\n";
  for ( element = 0; element < element_num; element++ )
  {
    gmsh << element + 1
         << "  " << element_type
         << "  " << tag_num
         << "  " << tag1
         << "  " << element + 1;
    if ( element_order == 20 )
    {
      for ( i = 0; i < element_order; i++ )
      {
        i2 = leo_to_gmsh[i];
        gmsh << "  " << element_node[i2+element*element_order];
      }
      gmsh << "\n";
    }
    else
    {
      for ( i = 0; i < element_order; i++ )
      {
        gmsh << "  " << element_node[i+element*element_order];
      }
      gmsh << "\n";
    }
  }
  gmsh << "$EndElements\n";

  gmsh.close ( );

  return;
}
//****************************************************************************80

int *i4mat_copy_new ( int m, int n, int a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_COPY_NEW copies an I4MAT to a "new" I4MAT.
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
//    27 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A1[M*N], the matrix to be copied.
//
//    Output, int I4MAT_COPY_NEW[M*N], the copy of A1.
//
{
  int *a2;
  int i;
  int j;

  a2 = new int[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
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
//    15 October 2014
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
//    Input, string TITLE, a title.
//
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }
    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(6) << i - 1 << "  ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = jlo;
    if ( j2lo < 1 )
    {
      j2lo = 1;
    }
    j2hi = jhi;
    if ( n < j2hi )
    {
      j2hi = n;
    }
    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j - 1 << ":";
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

void mesh_base_one ( int node_num, int element_order, int element_num, 
  int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    MESH_BASE_ONE ensures that the element definition is one-based.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2014
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
  const int i4_huge = 2147483647;
  int node;
  int node_max;
  int node_min;
  int order;

  node_min = + i4_huge;
  node_max = - i4_huge;
  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < element_order; order++ )
    {
      node = element_node[order+element*element_order];
      if ( node < node_min )
      {
        node_min = node;
      }
      if ( node_max < node )
      {
        node_max = node;
      }
    }
  }

  if ( node_min == 0 && node_max == node_num - 1 )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 0-based!\n";
    cout << "  This will be converted to 1-based.\n";
    for ( element = 0; element < element_num; element++ )
    {
      for ( order = 0; order < element_order; order++ )
      {
        element_node[order+element*element_order] =
          element_node[order+element*element_order] + 1;
      }
    }
  }
  else if ( node_min == 1 && node_max == node_num )
  {
    cout << "\n";
    cout << "MESH_BASE_ONE:\n";
    cout << "  The element indexing appears to be 1-based!\n";
    cout << "  No conversion is necessary.\n";
  }
  else
  {
    cout << "\n";
    cout << "MESH_BASE_ONE - Warning!\n";
    cout << "  The element indexing is not of a recognized type.\n";
    cout << "  NODE_MIN = " << node_min << "\n";
    cout << "  NODE_MAX = " << node_max << "\n";
    cout << "  NODE_NUM = " << node_num << "\n";
  }
  return;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Discussion:
//
//    The C++ math library provides the function fmax() which is preferred.
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
//  Discussion:
//
//    The C++ math library provides the function fmin() which is preferred.
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

double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
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
//    Input, string TITLE, a title.
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
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2014
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
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int i2lo_hi;
  int i2lo_lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

  if ( ilo < 1 )
  {
    i2lo_lo = 1;
  }
  else
  {
    i2lo_lo = ilo;
  }

  if ( ihi < m )
  {
    i2lo_hi = m;
  }
  else
  {
    i2lo_hi = ihi;
  }

  for ( i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;

    if ( m < i2hi )
    {
      i2hi = m;
    }
    if ( ihi < i2hi )
    {
      i2hi = ihi;
    }

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    if ( jlo < 1 )
    {
      j2lo = 1;
    }
    else
    {
      j2lo = jlo;
    }
    if ( n < jhi )
    {
      j2hi = n;
    }
    else
    {
      j2hi = jhi;
    }

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
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

bool s_begin ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
{
  int i;
  int n1;
  int n2;

  n1 = s1.length ( );
  n2 = s2.length ( );

  if ( n1 < n2 )
  {
    return false;
  }

  for ( i = 0; i < n2; i++ )
  {
    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
  return true;
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
    if ( s[n-1] != ' ' && s[n-1] != '\n' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

int s_to_i4 ( string s, int &last, bool &error )

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
//    Output, int &LAST, the last character of S used to make IVAL.
//
//    Output, bool &ERROR is TRUE if an error occurred.
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

  error = false;
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
        error = true;
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
        error = true;
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
        ival = 10 * ival + c - '0';
      }
      else
      {
        ival = isgn * ival;
        last = i - 1;
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
    last = s_len_trim ( s );
  }
  else
  {
    error = true;
    last = 0;
  }

  return ival;
}
//****************************************************************************80

double s_to_r8 ( string s, int &lchar, bool &error )

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
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
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
//    Output, int &LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool &ERROR, is true if an error occurred.
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
  error = false;
  r = 0.0;
  lchar = -1;
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
    c = s[lchar+1];
    lchar = lchar + 1;
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
        lchar = lchar + 1;
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
    if ( iterm == 1 || nchar <= lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && lchar + 1 == nchar )
  {
    lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    error = true;
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
//**********************************************************************

void timestamp ( )

//**********************************************************************
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
//    24 September 2003
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
