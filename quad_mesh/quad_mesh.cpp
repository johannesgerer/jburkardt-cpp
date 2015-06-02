# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>

using namespace std;

# include "quad_mesh.hpp"

//****************************************************************************80

int *adj_set_q4_mesh ( int node_num, int element_num,
  int element_node[], int element_neighbor[], int adj_num, int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_SET_Q4_MESH sets adjacencies in a triangulation.
//
//  Discussion:
//
//    This routine is called to set the adjacencies, after the
//    appropriate amount of memory has been set aside for storage.
//
//    The mesh is assumed to involve 4-node quadrilaterals.
//
//    Two nodes are "adjacent" if they are both nodes in some element.
//    Also, a node is considered to be adjacent to itself.
//
//    This routine can be used to create the compressed column storage
//    for a linear element finite element discretization of
//    Poisson's equation in two dimensions.
//
//  Diagram:
//
//         side 3
//       4-------3
//    s  |       |  s
//    i  |       |  i
//    d  |       |  d
//    e  |       |  e
//       |       |
//    4  |       |  2
//       |       |
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   20-21-22-23-24
//    |  |  |  |  |
//    |  |  |  |  |
//   15-16-17-18-19
//    |  |  |  |  |
//    |  |  |  |  |
//   10-11-12-13-14
//    |  |  |  |  |
//    |  |  |  |  |
//    5--6--7--8--9
//    |  |  |  |  |
//    |  |  |  |  |
//    0--1--2--3--4
//
//    A sample grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the nodes that
//    make up each element in counterclockwise order.
//
//    Input, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], for each side of
//    an element, lists the neighboring element, or -1 if there is
//    no neighbor.
//
//    Input, int ADJ_NUM, the number of adjacencies.
//
//    Input, int ADJ_ROW[NODE_NUM+1].  Information about column J is stored
//    in entries ADJ_ROW(J) through ADJ_ROW(J+1)-1 of ADJ.
//
//    Output, int ADJ_SET_Q4_MESH[ADJ_NUM], the adjacency information.
//
{
  int *adj;
  int *adj_copy;
  int k;
  int k1;
  int k2;
  int n1;
  int n2;
  int n3;
  int n4;
  int node;
  int element;
  int element2;
  int element_order = 4;

  adj = new int[adj_num];
  for ( k = 0; k < adj_num; k++ )
  {
    adj[k] = -1;
  }

  adj_copy = new int[node_num];
  for ( node = 0; node < node_num; node++ )
  {
    adj_copy[node] = adj_row[node];
  }
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj[adj_copy[node]] = node;
    adj_copy[node] = adj_copy[node] + 1;
  }
//
//  Examine each element.
//
  for ( element = 0; element < element_num; element++ )
  {
    n1 = element_node[0+element*element_order];
    n2 = element_node[1+element*element_order];
    n3 = element_node[2+element*element_order];
    n4 = element_node[3+element*element_order];
//
//  Add edges (1,3) and (2,4).  There is no need to check for redundancy,
//  since this is the only case when these nodes can share an element.
//
    adj[adj_copy[n1]] = n3;
    adj_copy[n1] = adj_copy[n1] + 1;
    adj[adj_copy[n3]] = n1;
    adj_copy[n3] = adj_copy[n3] + 1;

    adj[adj_copy[n2]] = n4;
    adj_copy[n2] = adj_copy[n2] + 1;
    adj[adj_copy[n4]] = n2;
    adj_copy[n4] = adj_copy[n4] + 1;
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
//  or if this element is the first of the pair in which the edge
//  occurs (ELEMENT < ELEMENT2).
//
    element2 = element_neighbor[0+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj[adj_copy[n1]] = n2;
      adj_copy[n1] = adj_copy[n1] + 1;
      adj[adj_copy[n2]] = n1;
      adj_copy[n2] = adj_copy[n2] + 1;
    }
//
//  Add edge (2,3).
//
    element2 = element_neighbor[1+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj[adj_copy[n2]] = n3;
      adj_copy[n2] = adj_copy[n2] + 1;
      adj[adj_copy[n3]] = n2;
      adj_copy[n3] = adj_copy[n3] + 1;
    }
//
//  Add edge (3,4).
//
    element2 = element_neighbor[2+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj[adj_copy[n4]] = n3;
      adj_copy[n4] = adj_copy[n4] + 1;
      adj[adj_copy[n3]] = n4;
      adj_copy[n3] = adj_copy[n3] + 1;
    }
//
//  Add edge (4,1).
//
    element2 = element_neighbor[3+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj[adj_copy[n1]] = n4;
      adj_copy[n1] = adj_copy[n1] + 1;
      adj[adj_copy[n4]] = n1;
      adj_copy[n4] = adj_copy[n4] + 1;
    }
  }
//
//  Ascending sort the entries for each node.
//
  for ( node = 0; node < node_num; node++ )
  {
    k1 = adj_row[node];
    k2 = adj_row[node+1]-1;
    i4vec_sort_heap_a ( k2+1-k1, adj+k1 );
  }

  delete [] adj_copy;

  return adj;
}
//****************************************************************************80

int adj_size_q4_mesh ( int node_num, int element_num, int element_node[], 
  int element_neighbor[], int adj_row[] )

//****************************************************************************80
//
//  Purpose:
//
//    ADJ_SIZE_Q4_MESH counts adjacencies in a Q4 mesh.
//
//  Discussion:
//
//    This routine is called to count the adjacencies, so that the
//    appropriate amount of memory can be set aside for storage when
//    the adjacency structure is created.
//
//    The mesh is assumed to involve 4-node quadrilaterals.
//
//    Two nodes are "adjacent" if they are both nodes in some quadrilateral.
//    Also, a node is considered to be adjacent to itself.
//
//  Diagram:
//
//         side 3
//       4-------3
//    s  |       |  s
//    i  |       |  i
//    d  |       |  d
//    e  |       |  e
//       |       |
//    4  |       |  2
//       |       |
//       1-------2
//
//         side 1
//
//    The local node numbering
//
//
//   20-21-22-23-24
//    |  |  |  |  |
//    |  |  |  |  |
//   15-16-17-18-19
//    |  |  |  |  |
//    |  |  |  |  |
//   10-11-12-13-14
//    |  |  |  |  |
//    |  |  |  |  |
//    5--6--7--8--9
//    |  |  |  |  |
//    |  |  |  |  |
//    0--1--2--3--4
//
//    A sample grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the nodes that
//    make up each element, in counterclockwise order.
//
//    Input, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], for each side of
//    a element, lists the neighboring elment, or -1 if there is
//    no neighbor.
//
//    Output, int ADJ_ROW[NODE_NUM+1], Information about column J is stored
//    in entries ADJ_ROW[J] through ADJ_ROW[J+1]-1 of ADJ.
//
//    Output, int ADJ_SIZE_Q4_MESH, the number of adjacencies.
//
{
  int adj_num;
  int element;
  int element_order = 4;
  int element2;
  int i;
  int n1;
  int n2;
  int n3;
  int n4;
  int node;

  adj_num = 0;
//
//  Set every node to be adjacent to itself.
//
  for ( node = 0; node < node_num; node++ )
  {
    adj_row[node] = 1;
  }
//
//  Examine each element.
//
  for ( element = 0; element < element_num; element++ )
  {
    n1 = element_node[0+element*element_order];
    n2 = element_node[1+element*element_order];
    n3 = element_node[2+element*element_order];
    n4 = element_node[3+element*element_order];
//
//  Add edge (1,3).
//
    adj_row[n1] = adj_row[n1] + 1;
    adj_row[n3] = adj_row[n3] + 1;
//
//  Add edge (2,4).
//
    adj_row[n2] = adj_row[n2] + 1;
    adj_row[n4] = adj_row[n4] + 1;
//
//  Add edge (1,2) if this is the first occurrence,
//  that is, if the edge (1,2) is on a boundary (ELEMENT2 <= 0)
//  or if this element is the first of the pair in which the edge
//  occurs (ELEMENT < ELEMENT2).
//
    element2 = element_neighbor[0+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj_row[n1] = adj_row[n1] + 1;
      adj_row[n2] = adj_row[n2] + 1;
    }
//
//  Add edge (2,3).
//
    element2 = element_neighbor[1+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj_row[n2] = adj_row[n2] + 1;
      adj_row[n3] = adj_row[n3] + 1;
    }
//
//  Add edge (3,4).
//
    element2 = element_neighbor[2+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj_row[n3] = adj_row[n3] + 1;
      adj_row[n4] = adj_row[n4] + 1;
    }
//
//  Add edge (4,1).
//
    element2 = element_neighbor[3+element*4];

    if ( element2 < 0 || element < element2 )
    {
      adj_row[n4] = adj_row[n4] + 1;
      adj_row[n1] = adj_row[n1] + 1;
    }
  }
//
//  We used ADJ_ROW to count the number of entries in each column.
//  Convert it to pointers into the ADJ array.
//
  for ( node = node_num; 1 <= node; node-- )
  {
    adj_row[node] = adj_row[node-1];
  }
  adj_row[0] = 0;
  for ( i = 1; i <= node_num; i++ )
  {
    adj_row[i] = adj_row[i] + adj_row[i-1];
  }
//
//  Finally, record the total number of adjacencies.
//
  adj_num = adj_row[node_num];

  return adj_num;
}
//****************************************************************************80

void area_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], double element_area[], double *mesh_area )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_Q4_MESH computes areas of elements in a Q4 mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, double NODE_XY[2*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the 
//    nodes that make up each element, in counterclockwise order.
//
//    Output, double ELEMENT_AREA[ELEMENT_NUM], the element areas.
//
//    Output, double *MESH_AREA, the mesh area.
//
{
  int dim;
  int element;
  int node;
  double q4[2*4];

  for ( element = 0; element < element_num; element++ )
  {
    for ( node = 0; node < 4; node++ )
    {
      for ( dim = 0; dim < 2; dim++ )
      {
        q4[dim+2*node] = node_xy[dim+2*element_node[node+4*element]];
      }
    }
    element_area[element] = area_quad ( q4 );
  }

  *mesh_area = r8vec_sum ( element_num, element_area );

  return;
}
//****************************************************************************80

double area_quad ( double quad_xy[2*4] )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_QUAD returns the area of a quadrilateral.
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
//    Input, double QUAD_XY[2*4], the coordinates of the nodes.
//
//    Output, double AREA_QUAD, the area.
//
{
  double area;
  double area1;
  double area2;
  double t1[2*3];
  double t2[2*3];

  t1[0+0*2] = quad_xy[0+0*2];
  t1[1+0*2] = quad_xy[1+0*2];
  t1[0+1*2] = quad_xy[0+1*2];
  t1[1+1*2] = quad_xy[1+1*2];
  t1[0+2*2] = quad_xy[0+2*2];
  t1[1+2*2] = quad_xy[1+2*2];

  area1 = triangle_area ( t1 );

  t2[0+0*2] = quad_xy[0+2*2];
  t2[1+0*2] = quad_xy[1+2*2];
  t2[0+1*2] = quad_xy[0+3*2];
  t2[1+1*2] = quad_xy[1+3*2];
  t2[0+2*2] = quad_xy[0+0*2];
  t2[1+2*2] = quad_xy[1+0*2];

  area2 = triangle_area ( t2 );

  if ( area1 < 0.0 || area2 < 0.0 )
  {
    t1[0+0*2] = quad_xy[0+1*2];
    t1[1+0*2] = quad_xy[1+1*2];
    t1[0+1*2] = quad_xy[0+2*2];
    t1[1+1*2] = quad_xy[1+2*2];
    t1[0+2*2] = quad_xy[0+3*2];
    t1[1+2*2] = quad_xy[1+3*2];

    area1 = triangle_area ( t1 );

    t2[0+0*2] = quad_xy[0+3*2];
    t2[1+0*2] = quad_xy[1+3*2];
    t2[0+1*2] = quad_xy[0+0*2];
    t2[1+1*2] = quad_xy[1+0*2];
    t2[0+2*2] = quad_xy[0+1*2];
    t2[1+2*2] = quad_xy[1+1*2];

    area2 = triangle_area ( t2 );

    if ( area1 < 0.0 || area2 < 0.0 )
    {
      cerr << "\n";
      cerr << "AREA_QUAD - Fatal error!\n";
      cerr << "  The quadrilateral nodes seem to be listed in\n";
      cerr << "  the wrong order, or the quadrilateral is\n";
      cerr << "  degenerate.\n";
      exit ( 1 );
    }
  }
  area = area1 + area2;

  return area;
}
//****************************************************************************80

void bandwidth ( int element_order, int element_num, int element_node[], 
  int *ml, int *mu, int *m )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH determines the bandwidth associated with the finite element mesh.
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
//    23 September 2006
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
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
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

int boundary_edge_count_q4_mesh ( int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOUNDARY_EDGE_COUNT_Q4_MESH counts the boundary edges.
//
//  Discussion:
//
//    This routine is given a Q4 mesh, an abstract list of sets of 4 nodes.
//    It is assumed that the nodes in each Q4 are listed
//    in a counterclockwise order, although the routine should work 
//    if the nodes are consistently listed in a clockwise order as well.
//
//    It is assumed that each edge of the mesh is either 
//    * an INTERIOR edge, which is listed twice, once with positive
//      orientation and once with negative orientation, or;
//    * a BOUNDARY edge, which will occur only once.
//
//    This routine should work even if the region has holes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes 
//    that make up the elements.  These should be listed in counterclockwise 
//    order.
//
//    Output, int BOUNDARY_EDGE_COUNT_Q4_MESH, the number of boundary 
//    edges.
//
{
  int boundary_edge_num;
  int e1;
  int e2;
  int *edge;
  int element;
  int interior_edge_num;
  int j;
  int m;
  int n;
  int unique_num;

  m = 2;
  n = 4 * element_num;
//
//  Set up the edge array.
//
  edge = new int[2*4*element_num];

  for ( element = 0; element < element_num; element++ )
  {
    edge[0+element*2+element_num*2*0] = element_node[0+element*4];
    edge[1+element*2+element_num*2*0] = element_node[1+element*4];

    edge[0+element*2+element_num*2*1] = element_node[1+element*4];
    edge[1+element*2+element_num*2*1] = element_node[2+element*4];

    edge[0+element*2+element_num*2*2] = element_node[2+element*4];
    edge[1+element*2+element_num*2*2] = element_node[3+element*4];

    edge[0+element*2+element_num*2*3] = element_node[3+element*4];
    edge[1+element*2+element_num*2*3] = element_node[0+element*4];
  }
//
//  In each column, force the smaller entry to appear first.
//
  for ( j = 0; j < n; j++ )
  {
    e1 = i4_min ( edge[0+2*j], edge[1+2*j] );
    e2 = i4_max ( edge[0+2*j], edge[1+2*j] );
    edge[0+2*j] = e1;
    edge[1+2*j] = e2;
  }
//
//  Ascending sort the column array.
//
  i4col_sort_a ( m, n, edge );
//
//  Get the number of unique columns in EDGE.
//
  unique_num = i4col_sorted_unique_count ( m, n, edge );

  interior_edge_num = 4 * element_num - unique_num;

  boundary_edge_num = 4 * element_num - 2 * interior_edge_num;

  delete [] edge;

  return boundary_edge_num;
}
//****************************************************************************80

int boundary_edge_count_euler_q4_mesh ( int node_num, int element_num, 
  int hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts boundary edges.
//
//  Discussion:
//
//    We assume we are given information about a quadrilateral mesh
//    of a set of nodes in the plane.
//
//    Given the number of nodes, elements and holes, we are going to apply
//    Euler's formula to determine the number of edges that lie on the
//    boundary of the set of nodes.
//
//    The number of faces, including the infinite face and internal holes, 
//    is ELEMENT_NUM + HOLE_NUM + 1.
//
//    Let BOUNDARY_NUM denote the number of edges on the boundary.
//    Each of the ELEMENT_NUM quadrilaterals uses four edges.  Every edge
//    occurs in two different elements, so the number of edges must be
//    ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2.
//
//    The number of nodes used in the mesh is NODE_NUM.
//
//    Euler's formula asserts that, for a simple connected figure in the
//    plane with no edge crossings, NODE_NUM nodes, EDGE_NUM edges and
//    FACE_NUM faces:
//
//      NODE_NUM - EDGE_NUM + FACE_NUM = 2
//
//    In our context, this becomes
//
//      NODE_NUM - ( 4 * ELEMENT_NUM + BOUNDARY_NUM ) / 2 
//      + ELEMENT_NUM + HOLE_NUM + 1 = 2
//
//    or
//
//      BOUNDARY_NUM = 2 * NODE_NUM + 2 * HOLE_NUM - 2 * ELEMENT_NUM - 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Marc de Berg, Marc Krevald, Mark Overmars, Otfried Schwarzkopf,
//    Computational Geometry, Section 9.1,
//    Springer, 2000.
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int HOLE_NUM, the number of internal holes.
//
//    Output, int BOUNDARY_EDGE_COUNT_EULER_Q4_MESH, the number of edges that 
//    lie on the boundary of the mesh.
//
{
  int boundary_num;

  boundary_num = 2 * node_num + 2 * hole_num - 2 * element_num - 2;

  return boundary_num;
}
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
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If the character was
//    'illegal', then DIGIT is -1.
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

void example1_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int element_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.
//
//  Discussion:
//
//    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
//    calling EXAMPLE1_Q4_MESH_SIZE first.
//
//   24---25---26---27---28
//    | 14 | 15 | 16 | 17 |
//   18---19---20---21---22---23
//    | 10 | -2 | 11 | 12 | 13 |
//   12---13---14---15---16---17
//    |  5 |  6 |  7 |  8 |  9 |
//    6----7----8----9---10---11
//    |  1 |  2 |  3 |  4 |
//    1----2----3----4----5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.  
//
//    Input, int ELEMENT_NUM, the number of elements.  
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the
//    nodes.
//
//    Output, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes
//    that make up the elements.
//
//    Output, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], the 
//    element neighbors on each side.  Negative values indicate edges that 
//    lie on the exterior.
//
{
# define ELEMENT_NUM_DATA 17
# define NODE_NUM_DATA 28

  int element_neighbor_data[4*ELEMENT_NUM_DATA] = {
       -1,  1,  4, -1, 
       -1,  2,  5,  0, 
       -1,  3,  6,  1, 
       -1, -1,  7,  2, 
        0,  5,  9, -1, 
        1,  6, -2,  4, 
        2,  7, 10,  5, 
        3,  8, 11,  6, 
       -1, -1, 12,  7, 
        4, -2, 13, -1, 
        6, 11, 15, -2, 
        7, 12, 16, 10, 
        8, -1, -1, 11, 
        9, 14, -1, -1, 
       -2, 15, -1, 13, 
       10, 16, -1, 14, 
       11, -1, -1, 15 };

  int element_node_data[4*ELEMENT_NUM_DATA] = {
     0,  1,  6,  5, 
     1,  2,  7,  6, 
     2,  3,  8,  7, 
     3,  4,  9,  8, 
     5,  6, 12, 11, 
     6,  7, 13, 12, 
     7,  8, 14, 13, 
     8,  9, 15, 14, 
     9, 10, 16, 15, 
    11, 12, 18, 17, 
    13, 14, 20, 19, 
    14, 15, 21, 20, 
    15, 16, 22, 21, 
    17, 18, 24, 23, 
    18, 19, 25, 24, 
    19, 20, 26, 25, 
    20, 21, 27, 26 };

  double node_xy_data[2*NODE_NUM_DATA] = {
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
       5.0, 1.0, 
       0.0, 2.0, 
       1.0, 2.0, 
       2.0, 2.0, 
       3.0, 2.0, 
       4.0, 2.0, 
       5.0, 2.0, 
       0.0, 3.0, 
       1.0, 3.0, 
       2.0, 3.0, 
       3.0, 3.0, 
       4.0, 3.0, 
       5.0, 3.0, 
       0.0, 4.0, 
       1.0, 4.0, 
       2.0, 4.0, 
       3.0, 4.0, 
       4.0, 4.0 };

  i4mat_copy ( 4, element_num, element_neighbor_data, element_neighbor );

  i4mat_copy ( 4, element_num, element_node_data, element_node );

  r8mat_copy ( 2, node_num, node_xy_data, node_xy );

  return;
# undef ELEMENT_NUM_DATA
# undef NODE_NUM_DATA
}
//****************************************************************************80

void example1_q4_mesh_size ( int *node_num, int *element_num, int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    EXAMPLE1_Q4_MESH_SIZE sets sizes for example #1 Q4 mesh
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.  
//
//    Output, int *ELEMENT_NUM, the number of elements. 
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *element_num = 17;
  *hole_num = 1;
  *node_num = 28;

  return;
}
//****************************************************************************80

void example2_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int element_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.
//
//  Discussion:
//
//    The region is a semicircle.  This example includes degenerate elements
//    (the first layer of elements is touching the origin, and so has a side
//    of length zero).  The elements are not parallelograms.  And the elements
//    vary in size.
//
//    Because of the treatment of node 1, algorithms for counting boundary 
//    edges may become "confused".
//
//    The appropriate values of NODE_NUM and ELEMENT_NUM can be found by
//    calling EXAMPLE1_Q4_MESH_SIZE first.
//
//   29---30---31---32---33---34---35---36---37
//    | 25 | 26 | 27 | 28 | 29 | 30 | 31 | 32 |
//   20---21---22---23---24---25---26---27---28
//    | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
//   11---12---13---14---15---16---17---18---19
//    |  9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
//    2----3----4----5----6----7----8----9---10
//    |  1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |
//    1----1----1----1----1----1----1----1----1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.  
//
//    Input, int ELEMENT_NUM, the number of elements.  
//
//    Output, double NODE_XY[2*NODE_NUM], the coordinates of the
//    nodes.
//
//    Output, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes
//    that make up the elements.
//
//    Output, int ELEMENT_NEIGHBOR[4*ELEMENT_NUM], the 
//    element neighbors on each side.  Negative values indicate edges that 
//    lie on the exterior.
//
{
  double a;
  int col;
  int element;
  int k;
  double pi = 3.141592653589793;
  double r;
  int row;

  k = 0;
  node_xy[0+k*2] = 0.0;
  node_xy[1+k*2] = 0.0;

  for ( row = 1; row <= 4; row++ )
  {
    r = ( double ) ( row );
    for ( col = 0; col <= 8; col++ )
    {
      a = ( double ) ( 8 - col ) * pi / 8.0;
      k = k + 1;
      node_xy[0+k*2] = r * cos ( a );
      node_xy[1+k*2] = r * sin ( a );
    }
  }

  element = 0;
  for ( row = 0; row <= 3; row++ )
  {
    for ( col = 0; col <= 7; col++ )
    {
      if ( row == 0 )
      {
        element_node[0+element*4] = 1;
        element_node[1+element*4] = 1;
        element_node[2+element*4] = col + 3;
        element_node[3+element*4] = col + 2;
      }
      else
      {
        element_node[0+element*4] = element_node[3+(element-8)*4];
        element_node[1+element*4] = element_node[2+(element-8)*4];
        element_node[2+element*4] = element_node[1+element*4] + 9;
        element_node[3+element*4] = element_node[0+element*4] + 9;
      }
      element = element + 1;
    }
  }

  element = 0;
  for ( row = 0; row <= 3; row++ )
  {
    for ( col = 0; col <= 7; col++ )
    {
      if ( row == 0 )
      {
        element_neighbor[0+element*4] = -1;
      }
      else
      {
        element_neighbor[0+element*4] = element - 8;
      }
      if ( col == 7 )
      {
        element_neighbor[1+element*4] = -1;
      }
      else
      {
        element_neighbor[1+element*4] = element + 1;
      }
      if ( row == 3 )
      {
        element_neighbor[2+element*4] = - 1;
      }
      else
      {
        element_neighbor[2+element*4] = element + 8;
      }
      if ( col == 0 )
      {
        element_neighbor[3+element*4] = - 1;
      }
      else
      {
        element_neighbor[3+element*4] = element - 1;
      }
      element = element + 1;
    }
  }
  return;
}
//****************************************************************************80

void example2_q4_mesh_size ( int *node_num, int *element_num, int *hole_num )

//****************************************************************************80
//
//  Purpose:
//
//    EXAMPLE2_Q4_MESH_SIZE sets sizes for example #2 Q4 mesh
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *NODE_NUM, the number of nodes.  
//
//    Output, int *ELEMENT_NUM, the number of elements. 
//
//    Output, int *HOLE_NUM, the number of holes.
//
{
  *element_num = 32;
  *hole_num = 0;
  *node_num = 37;

  return;
}
//****************************************************************************80

int file_column_count ( string input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed 
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char line[255];
//
//  Open the file.
//
  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    column_num = -1;
    cerr << "\n";
    cerr << "FILE_COLUMN_COUNT - Fatal error!\n";
    cerr << "  Could not open the file:\n";
    cerr << "  \"" << input_filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = true;
    break;

  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( input_filename.c_str ( ) );

    for ( ; ; )
    {
      input.getline ( line, sizeof ( line ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
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

  column_num = s_word_count ( line );

  return column_num;
}
//****************************************************************************80

int file_row_count ( string input_filename )

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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  char line[255];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename.c_str ( ) );

  if ( !input )
  {
    cerr << "\n";
    cerr << "FILE_ROW_COUNT - Fatal error!\n";
    cerr << "  Could not open the input file: \"" << input_filename << "\"\n";
    return (-1);
  }

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
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

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If 
//      NREM = I4_MODP ( I, J ) 
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
// 
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is 
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

int i4col_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_COMPARE compares columns I and J of an I4COL.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 4
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4COL_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], an array of N columns of vectors of length M.
//
//    Input, int I, J, the columns to be compared.
//    I and J must be between 1 and N.
//
//    Output, int I4COL_COMPARE, the results of the comparison:
//    -1, column I < column J,
//     0, column I = column J,
//    +1, column J < column I.
//
{
  int k;
//
//  Check.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index I = " << i << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < i )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index I = " << i << ".\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  Column index J = " << j << " is less than 1.\n";
    exit ( 1 );
  }

  if ( n < j )
  {
    cout << "\n";
    cout << "I4COL_COMPARE - Fatal error!\n";
    cout << "  N = " << n << " is less than column index J = " << j << ".\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  k = 1;

  while ( k <= m )
  {
    if ( a[k-1+(i-1)*m] < a[k-1+(j-1)*m] )
    {
      return (-1);
    }
    else if ( a[k-1+(j-1)*m] < a[k-1+(i-1)*m] )
    {
      return 1;
    }
    k = k + 1;
  }

  return 0;
}
//****************************************************************************80

void i4col_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORT_A ascending sorts the columns of an I4COL.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of N columns of M vectors;
//    On output, the columns of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4col_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4col_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

int i4col_sorted_unique_count ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SORTED_UNIQUE_COUNT counts unique elements in an I4COL.
//
//  Discussion:
//
//    The columns of the array may be ascending or descending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int A[M*N], a sorted array, containing
//    N columns of data.
//
//    Output, int I4COL_SORTED_UNIQUE_COUNT, the number of unique columns.
//
{
  int i;
  int j1;
  int j2;
  int unique_num;

  if ( n <= 0 )
  {
    unique_num = 0;
    return unique_num;
  }

  unique_num = 1;
  j1 = 0;

  for ( j2 = 1; j2 < n; j2++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j1*m] != a[i+j2*m] )
      {
        unique_num = unique_num + 1;
        j1 = j2;
        break;
      }
    }
  }

  return unique_num;
}
//****************************************************************************80

void i4col_swap ( int m, int n, int a[], int icol1, int icol2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4COL_SWAP swaps two columns of an I4COL.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based//  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int ICOL1, ICOL2, the two columns to swap.
//    These indices should be between 1 and N.
//
{
# define OFFSET 1

  int i;
  int t;
//
//  Check.
//
  if ( icol1 - OFFSET < 0 || n-1 < icol1 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL1 is out of range.\n";
    exit ( 1 );
  }

  if ( icol2 - OFFSET < 0 || n-1 < icol2 - OFFSET )
  {
    cout << "\n";
    cout << "I4COL_SWAP - Fatal error!\n";
    cout << "  ICOL2 is out of range.\n";
    exit ( 1 );
  }

  if ( icol1 == icol2 )
  {
    return;
  }
  for ( i = 0; i < m; i++ )
  {
    t                     = a[i+(icol1-OFFSET)*m];
    a[i+(icol1-OFFSET)*m] = a[i+(icol2-OFFSET)*m];
    a[i+(icol2-OFFSET)*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4mat_copy ( int m, int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_COPY copies one I4MAT to another.
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
//    Output, int A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
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
//    Input, string TITLE, a title.
//
{
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

void i4mat_write ( string output_filename, int m, int n, int table[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_WRITE writes an I4MAT file with no header.
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
//    Input, int TABLE[M*N], the table data.
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
    return;
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(10) << table[i+j*m];
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

int i4row_compare ( int m, int n, int a[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_COMPARE compares two rows of a integer array.
//
//  Discussion:
//
//    The two dimensional information is stored in a one dimensional array,
//    by columns.  The entry A(I,J) is stored in A[I+J*M].
//
//    The input arguments I and J are row indices.  They DO NOT use the
//    C convention of starting at 0, but rather start at 1.
//
//  Example:
//
//    Input:
//
//      M = 3, N = 4, I = 2, J = 3
//
//      A = (
//        1  2  3  4
//        5  6  7  8
//        9 10 11 12 )
//
//    Output:
//
//      I4ROW_COMPARE = -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int  A[M*N], the array of data.
//
//    Input, int I, J, the rows to be compared.
//    I and J must be between 1 and M.
//
//    Output, int I4ROW_COMPARE, the results of the comparison:
//    -1, row I < row J,
//     0, row I = row J,
//    +1, row J < row I.
//
{
  int k;
//
//  Check that I and J are legal.
//
  if ( i < 1 )
  {
    cout << "\n";
    cout << "I4ROW_COMPARE - Fatal error!\n";
    cout << "  Row index I is less than 1.\n";
    cout << "  I = " << i << "\n";
    exit ( 1 );
  }
  else if ( m < i )
  {
    cout << "\n";
    cout << "I4ROW_COMPARE - Fatal error!\n";
    cout << "  Row index I is out of bounds.\n";
    cout << "  I = " << i << "\n";
    cout << "  Maximum legal value is M = " << m << "\n";
    exit ( 1 );
  }

  if ( j < 1 )
  {
    cout << "\n";
    cout << "I4ROW_COMPARE - Fatal error!\n";
    cout << "  Row index J is less than 1.\n";
    cout << "  J = " << j << "\n";
    exit ( 1 );
  }
  else if ( m < j )
  {
    cout << "\n";
    cout << "I4ROW_COMPARE - Fatal error!\n";
    cout << "  Row index J is out of bounds.\n";
    cout << "  J = " << j << "\n";
    cout << "  Maximum legal value is M = " << m << "\n";
    exit ( 1 );
  }

  if ( i == j )
  {
    return 0;
  }

  for ( k = 0; k < n; k++ )
  {
    if ( a[(i-1)+k*m] < a[(j-1)+k*m] )
    {
      return -1;
    }
    else if ( a[(j-1)+k*m] < a[(i-1)+k*m] )
    {
      return +1;
    }
  }

  return 0;
}
//****************************************************************************80

void i4row_sort_a ( int m, int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_SORT_A ascending sorts the rows of an I4ROW.
//
//  Discussion:
//
//    In lexicographic order, the statement "X < Y", applied to two
//    vectors X and Y of length M, means that there is some index I, with
//    1 <= I <= M, with the property that
//
//      X(J) = Y(J) for J < I,
//    and
//      X(I) < Y(I).
//
//    In other words, X is less than Y if, at the first index where they
//    differ, the X value is less than the Y value.
//
//  Example:
//
//    Input:
//
//      M = 5, N = 3
//
//      A =
//        3  2  1
//        2  4  3
//        3  1  8
//        2  4  2
//        1  9  9
//
//    Output:
//
//      A =
//        1  9  9
//        2  4  2
//        2  4  3
//        3  1  8
//        3  2  1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of A.
//
//    Input, int N, the number of columns of A.
//
//    Input/output, int A[M*N].
//    On input, the array of M rows of N-vectors.
//    On output, the rows of A have been sorted in ascending
//    lexicographic order.
//
{
  int i;
  int indx;
  int isgn;
  int j;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( m, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      i4row_swap ( m, n, a, i, j );
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4row_compare ( m, n, a, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4row_swap ( int m, int n, int a[], int irow1, int irow2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4ROW_SWAP swaps two rows of an I4ROW.
//
//  Discussion:
//
//    The two dimensional information is stored as a one dimensional
//    array, by columns.
//
//    The row indices are 1 based, NOT 0 based//  However, a preprocessor
//    variable, called OFFSET, can be reset from 1 to 0 if you wish to
//    use 0-based indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int A[M*N], an array of data.
//
//    Input, int IROW1, IROW2, the two rows to swap.
//    These indices should be between 1 and M.
//
{
# define OFFSET 1

  int j;
  int t;
//
//  Check.
//
  if ( irow1 < 0+OFFSET || m-1+OFFSET < irow1 )
  {
    cout << "\n";
    cout << "I4ROW_SWAP - Fatal error!\n";
    cout << "  IROW1 is out of range.\n";
    exit ( 1 );
  }

  if ( irow2 < 0+OFFSET || m-1+OFFSET < irow2 )
  {
    cout << "\n";
    cout << "I4ROW_SWAP - Fatal error!\n";
    cout << "  IROW2 is out of range.\n";
    exit ( 1 );
  }

  if ( irow1 == irow2 )
  {
    return;
  }
  for ( j = 0; j < n; j++ )
  {
    t                   = a[irow1-OFFSET+j*m];
    a[irow1-OFFSET+j*m] = a[irow2-OFFSET+j*m];
    a[irow2-OFFSET+j*m] = t;
  }

  return;
# undef OFFSET
}
//****************************************************************************80

void i4vec_heap_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                        
//            A(1)         A(2)
//                             
//      A(3)       A(4)  A(5) A(6)
//                      
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  { 
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ; ; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }
      }
    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;
  }

  return;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8) << i 
         << "  " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

void i4vec_sort_heap_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;
  }
  return;
}
//****************************************************************************80

int *i4vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ZERO_NEW creates and zeroes an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, int I4VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0;
  }
  return a;
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

int *neighbor_elements_q4_mesh ( int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    NEIGHBOR_ELEMENTS_Q4_MESH determines element neighbors in a Q4 mesh.
//
//  Discussion:
//
//    A quadrilateral mesh of a set of nodes can be completely described by
//    the coordinates of the nodes, and the list of nodes that make up
//    each element.  However, in some cases, it is necessary to know
//    element adjacency information, that is, which element, if any,
//    is adjacent to a given element on a particular side.
//
//    This routine creates a data structure recording this information.
//
//    The primary amount of work occurs in sorting a list of 4 * ELEMENT_NUM
//    data items.
//
//    Note that COL is a work array allocated dynamically inside this
//    routine.  It is possible, for very large values of ELEMENT_NUM,
//    that the necessary amount of memory will not be accessible, and the
//    routine will fail.  This is a limitation of the implementation of
//    dynamic arrays in FORTRAN90.  One way to get around this would be
//    to require the user to declare ROW in the calling routine
//    as an allocatable array, get the necessary memory explicitly with
//    an ALLOCATE statement, and then pass ROW into this routine.
//
//    Of course, the point of dynamic arrays was to make it easy to
//    hide these sorts of temporary work arrays from the poor user!
//
//    This routine was revised to store the edge data in a column
//    array rather than a row array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes that make up each element.
//
//    Output, int NEIGHBOR_ELEMENTS_Q4_MESH[4*ELEMENT_NUM], the elements that are direct 
//    neighbors of a given element, or -1 if there is no neighbor on that side.
//
{
  int *col;
  int element;
  int *element_neighbor;
  int element_order = 4;
  int element1;
  int element2;
  int i;
  int icol;
  int j;
  int k;
  int l;
  int side1;
  int side2;

  element_neighbor = new int[4*element_num];
  col = new int[4*4*element_num];
//
//  Step 1.
//  From the list of nodes for element E, of the form: (I,J,K,L)
//  construct the four neighbor relations:
//
//    (I,J,0,E) or (J,I,0,E),
//    (J,K,1,E) or (K,J,1,E),
//    (K,L,2,E) or (L,K,2,E)
//    (L,I,3,E) or (I,L,3,E)
//
//  where we choose (I,J,0,E) if I < J, or else (J,I,0,E)
//
  for ( element = 0; element < element_num; element++ )
  {
    i = element_node[0+element*element_order];
    j = element_node[1+element*element_order];
    k = element_node[2+element*element_order];
    l = element_node[3+element*element_order];

    col[0+0*4+16*element] = i4_min ( i, j );
    col[1+0*4+16*element] = i4_max ( i, j );
    col[2+0*4+16*element] = 0;
    col[3+0*4+16*element] = element;

    col[0+1*4+16*element] = i4_min ( j, k );
    col[1+1*4+16*element] = i4_max ( j, k );
    col[2+1*4+16*element] = 1;
    col[3+1*4+16*element] = element;

    col[0+2*4+16*element] = i4_min ( k, l );
    col[1+2*4+16*element] = i4_max ( k, l );
    col[2+2*4+16*element] = 2;
    col[3+2*4+16*element] = element;

    col[0+3*4+16*element] = i4_min ( l, i );
    col[1+3*4+16*element] = i4_max ( l, i );
    col[2+3*4+16*element] = 3;
    col[3+3*4+16*element] = element;
  }
//
//  Step 2. Perform an ascending dictionary sort on the neighbor relations.
//  We only intend to sort on rows 1 and 2; the routine we call here
//  sorts on rows 1 through 4 but that won't hurt us.
//
//  What we need is to find cases where two elements share an edge.
//  Say they share an edge defined by the nodes I and J.  Then there are
// two columns of COL that start out ( I, J, ?, ? ).  By sorting COL,
//  we make sure that these two columns occur consecutively.  That will
//  make it easy to notice that the elements are neighbors.
//
  i4col_sort_a ( 4, 4*element_num, col );
//
//  Step 3. Neighboring elements show up as consecutive columns with
//  identical first two entries.  Whenever you spot this happening,
//  make the appropriate entries in ELEMENT_NEIGHBOR.
//
  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      element_neighbor[i+j*4] = -1;
    }
  }

  icol = 0;

  for ( ; ; )
  {
    if ( 4 * element_num <= icol )
    {
      break;
    }

    if ( col[0+icol*4] != col[0+(icol+1)*4] || 
         col[1+icol*4] != col[1+(icol+1)*4] )
    {
      icol = icol + 1;
      continue;
    }

    side1    = col[2+icol*4];
    element1 = col[3+icol*4];

    side2    = col[2+(icol+1)*4];
    element2 = col[3+(icol+1)*4];

    element_neighbor[side1+element1*4] = element2;
    element_neighbor[side2+element2*4] = element1;

    icol = icol + 2;
  }
 
  delete [] col;

  return element_neighbor;
}
//****************************************************************************80

int *node_order_q4_mesh ( int element_num, int element_node[], int node_num )

//****************************************************************************80
//
//  Purpose:
//
//    NODE_ORDER_Q4_MESH determines the order of nodes in a Q4 mesh.
//
//  Discussion:
//
//    The order of a node is the number of elements that use that node
//    as a vertex.
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
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], 
//    the nodes that make up the elements.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int NODE_ORDER_Q4_MESH[NODE_NUM], the order of each node.
//
{
  int element;
  int i;
  int node;
  int *node_order;

  node_order = i4vec_zero_new ( node_num );

  for ( element = 0; element < element_num; element++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      node = element_node[i+element*4];
      if ( node < 0 || node_num <= node )
      {
        cerr << "\n";
        cerr << "NODE_ORDER_Q4_MESH - Fatal error!\n";
        cerr << "  Illegal entry in ELEMENT_NODE.\n";
        exit ( 1 );
      }
      else
      {
        node_order[node] = node_order[node] + 1;
      }
    }
  }
  return node_order;
}
//****************************************************************************80

void plot_q4_mesh ( int node_num, int element_num, double node_xy[], 
  int element_node[], int node_show, int element_show, string output_filename )

//****************************************************************************80
//
//  Purpose:
//
//    PLOT_Q4_MESH plots a Q4 mesh.
//
//  Discussion:
//
//    The triangulation is most usually a Delaunay triangulation,
//    but this is not necessary.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int ELEMENT_NODE[4*ELEMENT_NUM], the nodes that form the elements.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int ELEMENT_SHOW:
//    0, do not show elements;
//    1, show elements;
//    2, show elements and label them.
//
//    Input, string OUTPUT_FILENAME, the name of the output file.
//
{
  double ave_x;
  double ave_y;
  int circle_size;
  int delta;
  int e;
  int element;
  int element_order = 4;
  int i;
  int node;
  ofstream output_unit;
  double x_max;
  double x_min;
  int x_ps;
  int x_ps_max = 576;
  int x_ps_max_clip = 594;
  int x_ps_min = 36;
  int x_ps_min_clip = 18;
  double x_scale;
  double y_max;
  double y_min;
  int y_ps;
  int y_ps_max = 666;
  int y_ps_max_clip = 684;
  int y_ps_min = 126;
  int y_ps_min_clip = 108;
  double y_scale;
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( x_max < node_xy[0+node*2] )
     {
       x_max = node_xy[0+node*2];
     }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[0+node*2] < x_min )
     {
       x_min = node_xy[0+node*2];
     }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( y_max < node_xy[1+node*2] )
     {
       y_max = node_xy[1+node*2];
     }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_xy[1+node*2] < y_min )
     {
       y_min = node_xy[1+node*2];
     }
  }
  y_scale = y_max - y_min;

  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  if ( x_scale < y_scale )
  {
    delta = r8_nint ( ( double ) ( x_ps_max - x_ps_min )
      * ( y_scale - x_scale ) / ( 2.0 * y_scale ) );

    x_ps_max = x_ps_max - delta;
    x_ps_min = x_ps_min + delta;

    x_ps_max_clip = x_ps_max_clip - delta;
    x_ps_min_clip = x_ps_min_clip + delta;

    x_scale = y_scale;
  }
  else if ( y_scale < x_scale )
  {
    delta = r8_nint ( ( double ) ( y_ps_max - y_ps_min ) 
      * ( x_scale - y_scale ) / ( 2.0 * x_scale ) );

    y_ps_max = y_ps_max - delta;
    y_ps_min = y_ps_min + delta;

    y_ps_max_clip = y_ps_max_clip - delta;
    y_ps_min_clip = y_ps_min_clip + delta;

    y_scale = x_scale;
  }

  output_unit.open ( output_filename.c_str() );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "PLOT_Q4_MESH - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  output_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  output_unit << "%%Creator: plot_q4_mesh.C\n";
  output_unit << "%%Title: " << output_filename << "\n";

  output_unit << "%%Pages: 1\n";
  output_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  output_unit << "%%Document-Fonts: Times-Roman\n";
  output_unit << "%%LanguageLevel: 1\n";
  output_unit << "%%EndComments\n";
  output_unit << "%%BeginProlog\n";
  output_unit << "/inch {72 mul} def\n";
  output_unit << "%%EndProlog\n";
  output_unit << "%%Page:      1     1\n";
  output_unit << "save\n";
  output_unit << "%\n";
  output_unit << "% Set the RGB line color to very light gray.\n";
  output_unit << "%\n";
  output_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  output_unit << "%\n";
  output_unit << "% Draw a gray border around the page.\n";
  output_unit << "%\n";
  output_unit << "newpath\n";
  output_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  output_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  output_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  output_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  output_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  output_unit << "stroke\n";
  output_unit << "%\n";
  output_unit << "% Set RGB line color to black.\n";
  output_unit << "%\n";
  output_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  output_unit << "%\n";
  output_unit << "%  Set the font and its size:\n";
  output_unit << "%\n";
  output_unit << "/Times-Roman findfont\n";
  output_unit << "0.50 inch scalefont\n";
  output_unit << "setfont\n";
  output_unit << "%\n";
  output_unit << "%  Print a title:\n";
  output_unit << "%\n";
  output_unit << "%  210  702 moveto\n";
  output_unit << "%(Pointset) show\n";
  output_unit << "%\n";
  output_unit << "% Define a clipping polygon\n";
  output_unit << "%\n";
  output_unit << "newpath\n";
  output_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  output_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  output_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  output_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  output_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  output_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( node_num <= 200 )
  {
    circle_size = 5;
  }
  else if ( node_num <= 500 )
  {
    circle_size = 4;
  }
  else if ( node_num <= 1000 )
  {
    circle_size = 3;
  }
  else if ( node_num <= 5000 )
  {
    circle_size = 2;
  }
  else
  {
    circle_size = 1;
  }

  if ( 1 <= node_show )
  {
    output_unit << "%\n";
    output_unit << "%  Draw filled dots at each node:\n";
    output_unit << "%\n";
    output_unit << "%  Set the color to blue:\n";
    output_unit << "%\n";
    output_unit << "0.000  0.150  0.750  setrgbcolor\n";
    output_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max                     - y_min ) );

      output_unit << "newpath  " 
        << x_ps << "  " 
        << y_ps << "  "
        << circle_size << " 0 360 arc closepath fill\n";
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    output_unit << "%\n";
    output_unit << "%  Label the nodes:\n";
    output_unit << "%\n";
    output_unit << "%  Set the color to darker blue:\n";
    output_unit << "%\n";
    output_unit << "0.000  0.250  0.850  setrgbcolor\n";
    output_unit << "/Times-Roman findfont\n";
    output_unit << "0.20 inch scalefont\n";
    output_unit << "setfont\n";

    output_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    { 
      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max                     - y_min ) );

      output_unit << "newpath  " 
        << x_ps     << "  " 
        << y_ps + 5 << "  moveto ("
        << node+1   << ") show\n";
    }
  }
//
//  Draw the elements.
//
  if ( 1 <= element_show )
  {
    output_unit << "%\n";
    output_unit << "%  Set the RGB color to red.\n";
    output_unit << "%\n";
    output_unit << "0.900  0.200  0.100 setrgbcolor\n";
    output_unit << "%\n";
    output_unit << "%  Draw the elements.\n";
    output_unit << "%\n";

    for ( element = 0; element < element_num; element++ )
    {
      output_unit << "newpath\n";

      for ( i = 1; i <= element_order+1; i++ )
      {
        e = i4_wrap ( i, 1, element_order );

        node = element_node[e-1+element*element_order] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        if ( i == 1 )
        {
          output_unit << x_ps << "  " << y_ps << "  moveto\n";
        } 
        else
        {
          output_unit << x_ps << "  " << y_ps << "  lineto\n";
        }
      }
      output_unit << "stroke\n";
    }
  }
//
//  Label the elements.
//
  if ( 2 <= element_show )
  {
    output_unit << "%\n";
    output_unit << "%  Label the elements.\n";
    output_unit << "%\n";
    output_unit << "%  Set the RGB color to darker red.\n";
    output_unit << "%\n";
    output_unit << "0.950  0.250  0.150 setrgbcolor\n";
    output_unit << "/Times-Roman findfont\n";
    output_unit << "0.20 inch scalefont\n";
    output_unit << "setfont\n";
    output_unit << "%\n";

    for ( element = 0; element < element_num; element++ )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( i = 1; i <= element_order; i++ )
      {
        node = element_node[i-1+element*element_order] - 1;
        ave_x = ave_x + node_xy[0+node*2];
        ave_y = ave_y + node_xy[1+node*2];
      }
      ave_x = ave_x / ( double ) ( element_order );
      ave_y = ave_y / ( double ) ( element_order );

      x_ps = ( int ) (
        ( ( x_max - ave_x         ) * ( double ) ( x_ps_min )  
        + (       + ave_x - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max         - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - ave_y         ) * ( double ) ( y_ps_min )  
        + (         ave_y - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max         - y_min ) );

      output_unit << x_ps << "  " 
                << y_ps << "  moveto ("
                << element+1 << ") show\n";
    }
  }

  output_unit << "%\n";
  output_unit << "restore  showpage\n";
  output_unit << "%\n";
  output_unit << "%  End of page.\n";
  output_unit << "%\n";
  output_unit << "%%Trailer\n";
  output_unit << "%%EOF\n";

  output_unit.close ( );

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

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

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

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to an R8.
//
//  Example:
//
//        X         Value
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the value.
//
//    Output, int R8_NINT, the nearest integer to X.
//
{
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r8_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r8_abs ( x ) + 0.5 );
  }

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
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
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_copy ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY copies one R8MAT to another.
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
//    16 October 2005
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
//    Output, double A2[M*N], the copy of A1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return;
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

void r8mat_mm ( int n1, int n2, int n3, double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM multiplies two matrices.
//
//  Discussion: 							    
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double C[N1*N3], the product matrix C = A * B.
//
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n1; i ++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
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
//    10 September 2009
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
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
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
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

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

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
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
//    Input, double TABLE[M*N], the table data.
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
    return;
  }
//
//  Write the data.
//  For greater precision, try
//
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(10) << table[i+j*m];
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

void r8vec_bracket ( int n, double x[], double xval, int *left, 
  int *right )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1]; 
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ ) 
  {
    if ( xval < x[i-1] ) 
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
}
//****************************************************************************80

void reference_to_physical_q4 ( double q4[2*4], int n, double rs[], 
  double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_Q4 maps Q4 reference points to physical points.
//
//  Discussion:
//
//    XY(R,S) = XY(0,0) * (1-R) * (1-S)
//            + XY(1,0) *    R  * (1-S)
//            + XY(1,1) *    R  *    S
//            + XY(0,1) * (1-R) *    S
//
//  Reference Element Q4:
//
//    |
//    1  4-----3
//    |  |     |
//    |  |     |
//    S  |     |
//    |  |     |
//    |  |     |
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double Q4[2*4], the coordinates of the vertices.
//    The vertices are assumed to be the images of the reference vertices
//    (0,0), (1,0), (1,1) and (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double RS[2*N], (R,S) points in the reference element.
//
//    Output, double XY[2*N], (X,Y) points in the physical element.
//
{
  int j;
  double *psi;

  psi = new double[4*n];

  for ( j = 0; j < n; j++ )
  {
    psi[0+j*2] = ( 1.0 - rs[0+j*2] ) * ( 1.0 - rs[1+j*2] );
    psi[1+j*2] =         rs[0+j*2]   * ( 1.0 - rs[1+j*2] );
    psi[2+j*2] =         rs[0+j*2]   *         rs[1+j*2];
    psi[3+j*2] = ( 1.0 - rs[0+j*2] ) *         rs[1+j*2];
  }

  r8mat_mm ( 2, 4, n, q4, psi, xy );

  delete [] psi;

  return;
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

void sample_q4_mesh ( int node_num, double node_xy[], int element_num, 
  int element_node[], int sample_num, int *seed, double sample_xy[], 
  int sample_element[] )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_Q4_MESH returns random points in a Q4 mesh.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY(2,NODE_NUM), the coordinates of the nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE(4,ELEMENT_NUM), the nodes
//    that form the elements.
//
//    Input, int SAMPLE_NUM, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random
//     number generator.
//
//    Output, double SAMPLE_XY(2,SAMPLE_NUM), the sample points.
//
//    Output, int SAMPLE_ELEMENT(SAMPLE_NUM), the elements from
//    which each point was drawn.
//
{
  double area;
  double *area_cum;
  double area_total;
  int element;
  int i1;
  int i2;
  int i3;
  int i4;
  int left;
  double quad_xy[2*4];
  double r;
  int right;
  int sample;
//
//  Compute the areas of the quadrilaterals.
//
  area_cum = new double[element_num+1];

  area_cum[0] = 0.0;

  for ( element = 1; element <= element_num; element++ )
  {
    i1 = element_node[0+(element-1)*4];
    i2 = element_node[1+(element-1)*4];
    i3 = element_node[2+(element-1)*4];
    i4 = element_node[3+(element-1)*4];

    quad_xy[0+0*2] = node_xy[0+i1*2];
    quad_xy[1+0*2] = node_xy[1+i1*2];
    quad_xy[0+1*2] = node_xy[0+i2*2];
    quad_xy[1+1*2] = node_xy[1+i2*2];
    quad_xy[0+2*2] = node_xy[0+i3*2];
    quad_xy[1+2*2] = node_xy[1+i3*2];
    quad_xy[0+3*2] = node_xy[0+i4*2];
    quad_xy[1+3*2] = node_xy[1+i4*2];

    area = area_quad ( quad_xy );

    area_cum[element] = area_cum[element-1] + area;
  }

  area_total = area_cum[element_num];

  for ( element = 0; element <= element_num; element++ )
  {
    area_cum[element] = area_cum[element] / area_total;
  }
//
//  A random value R indicates the corresponding quadrilateral whose
//  cumulative relative area first includes the number R.
//
  for ( sample = 0; sample < sample_num; sample++ )
  {
    r = r8_uniform_01 ( seed );

    r8vec_bracket ( element_num + 1, area_cum, r, &left, &right );

    element = right - 1;

    i1 = element_node[0+(element-1)*4];
    i2 = element_node[1+(element-1)*4];
    i3 = element_node[2+(element-1)*4];
    i4 = element_node[3+(element-1)*4];

    quad_xy[0+0*2] = node_xy[0+i1*2];
    quad_xy[1+0*2] = node_xy[1+i1*2];
    quad_xy[0+1*2] = node_xy[0+i2*2];
    quad_xy[1+1*2] = node_xy[1+i2*2];
    quad_xy[0+2*2] = node_xy[0+i3*2];
    quad_xy[1+2*2] = node_xy[1+i3*2];
    quad_xy[0+3*2] = node_xy[0+i4*2];
    quad_xy[1+3*2] = node_xy[1+i4*2];

    sample_quad ( quad_xy, 1, seed, sample_xy+sample*2 );

    sample_element[sample] = element;
  }

  delete [] area_cum;

  return;
}
//****************************************************************************80

void sample_quad ( double quad_xy[2*4], int n, int *seed, double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_QUAD returns random points in a quadrilateral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double QUAD_XY[2*4], the coordinates of the nodes.
//
//    Input, int N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random 
//     number generator.
//
//    Output, double XY[2*N], the sample points.
//
{
  double area1;
  double area2;
  double area_total;
  int i;
  double r;
  double t1[2*3];
  double t2[2*3];

  t1[0+0*2] = quad_xy[0+0*2];
  t1[1+0*2] = quad_xy[1+0*2];
  t1[0+1*2] = quad_xy[0+1*2];
  t1[1+1*2] = quad_xy[1+1*2];
  t1[0+2*2] = quad_xy[0+2*2];
  t1[1+2*2] = quad_xy[1+2*2];

  area1 = triangle_area ( t1 );

  t2[0+0*2] = quad_xy[0+2*2];
  t2[1+0*2] = quad_xy[1+2*2];
  t2[0+1*2] = quad_xy[0+3*2];
  t2[1+1*2] = quad_xy[1+3*2];
  t2[0+2*2] = quad_xy[0+0*2];
  t2[1+2*2] = quad_xy[1+0*2];

  area2 = triangle_area ( t2 );

  if ( area1 < 0.0 || area2 < 0.0 )
  {
    t1[0+0*2] = quad_xy[0+1*2];
    t1[1+0*2] = quad_xy[1+1*2];
    t1[0+1*2] = quad_xy[0+2*2];
    t1[1+1*2] = quad_xy[1+2*2];
    t1[0+2*2] = quad_xy[0+3*2];
    t1[1+2*2] = quad_xy[1+3*2];

    area1 = triangle_area ( t1 );

    t2[0+0*2] = quad_xy[0+3*2];
    t2[1+0*2] = quad_xy[1+3*2];
    t2[0+1*2] = quad_xy[0+0*2];
    t2[1+1*2] = quad_xy[1+0*2];
    t2[0+2*2] = quad_xy[0+1*2];
    t2[1+2*2] = quad_xy[1+1*2];

    area2 = triangle_area ( t2 );

    if ( area1 < 0.0 || area2 < 0.0 )
    {
      cerr << "\n";
      cerr << "SAMPLE_QUAD - Fatal error!\n";
      cerr << "  The quadrilateral nodes seem to be listed in\n";
      cerr << "  the wrong order, or the quadrilateral is\n";
      cerr << "  degenerate.\n";
      exit ( 1 );
    }
  }
  area_total = area1 + area2;
//
//  Choose a triangle at random, weighted by the areas.
//  Then choose a point in that triangle.
//
  for ( i = 0; i < n; i++ )
  {
    r = r8_uniform_01 ( seed );

    if ( r * area_total < area1 )
    {
      triangle_sample ( t1, 1, seed, xy+i*2 );
    }
    else
    {
      triangle_sample ( t2, 1, seed, xy+i*2 );
    }
  }

  return;
}
//****************************************************************************80

double *sample_quad_new ( double quad_xy[2*4], int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_QUAD_NEW returns random points in a quadrilateral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 February 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double QUAD_XY[2*4], the coordinates of the nodes.
//
//    Input, int N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random 
//     number generator.
//
//    Output, double SAMPLE_QUAD[2*N], the sample points.
//
{
  double area1;
  double area2;
  double area_total;
  int i;
  double r;
  double t1[2*3];
  double t2[2*3];
  double *xy;

  t1[0+0*2] = quad_xy[0+0*2];
  t1[1+0*2] = quad_xy[1+0*2];
  t1[0+1*2] = quad_xy[0+1*2];
  t1[1+1*2] = quad_xy[1+1*2];
  t1[0+2*2] = quad_xy[0+2*2];
  t1[1+2*2] = quad_xy[1+2*2];

  area1 = triangle_area ( t1 );

  t2[0+0*2] = quad_xy[0+2*2];
  t2[1+0*2] = quad_xy[1+2*2];
  t2[0+1*2] = quad_xy[0+3*2];
  t2[1+1*2] = quad_xy[1+3*2];
  t2[0+2*2] = quad_xy[0+0*2];
  t2[1+2*2] = quad_xy[1+0*2];

  area2 = triangle_area ( t2 );

  if ( area1 < 0.0 || area2 < 0.0 )
  {
    t1[0+0*2] = quad_xy[0+1*2];
    t1[1+0*2] = quad_xy[1+1*2];
    t1[0+1*2] = quad_xy[0+2*2];
    t1[1+1*2] = quad_xy[1+2*2];
    t1[0+2*2] = quad_xy[0+3*2];
    t1[1+2*2] = quad_xy[1+3*2];

    area1 = triangle_area ( t1 );

    t2[0+0*2] = quad_xy[0+3*2];
    t2[1+0*2] = quad_xy[1+3*2];
    t2[0+1*2] = quad_xy[0+0*2];
    t2[1+1*2] = quad_xy[1+0*2];
    t2[0+2*2] = quad_xy[0+1*2];
    t2[1+2*2] = quad_xy[1+1*2];

    area2 = triangle_area ( t2 );

    if ( area1 < 0.0 || area2 < 0.0 )
    {
      cerr << "\n";
      cerr << "SAMPLE_QUAD - Fatal error!\n";
      cerr << "  The quadrilateral nodes seem to be listed in\n";
      cerr << "  the wrong order, or the quadrilateral is\n";
      cerr << "  degenerate.\n";
      exit ( 1 );
    }
  }
  area_total = area1 + area2;
//
//  Choose a triangle at random, weighted by the areas.
//  Then choose a point in that triangle.
//
  xy = new double[2*n];

  for ( i = 0; i < n; i++ )
  {
    r = r8_uniform_01 ( seed );

    if ( r * area_total < area1 )
    {
      triangle_sample ( t1, 1, seed, xy+i*2 );
    }
    else
    {
      triangle_sample ( t2, 1, seed, xy+i*2 );
    }
  }

  return xy;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 ) 
    {
      if ( isgn < 0 ) 
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn ) 
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      if ( n1 == 1 ) 
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else 
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }
    k = k - 1;
    k1 = k;
  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 ) 
  {
    k1 = k;
  }

  for ( ; ; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 ) 
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 ) 
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 ) 
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 ) 
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else 
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
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
//    03 October 2003
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

double triangle_area ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed 
//    area of a triangle, it is possible to easily compute the area 
//    of a nonconvex polygon as the sum of the (possibly negative) 
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA, the area of the triangle.
//
{
  double area;

  area = 0.5 * ( 
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) + 
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) + 
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );
 
  return area;
}
//****************************************************************************80

void triangle_sample ( double t[2*3], int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE returns random points in a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
# define DIM_NUM 2

  double alpha;
  double beta;
  int j;
  double r;
  double p12[DIM_NUM];
  double p13[DIM_NUM];

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    alpha = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    p12[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+1*2];
    p12[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+1*2];

    p13[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+2*2];;
    p13[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+2*2];;
//
//  Now choose, uniformly at random, a point on the line L.
//
    beta = r8_uniform_01 ( seed );

    p[0+j*2] = ( 1.0 - beta ) * p12[0] + beta * p13[0];
    p[1+j*2] = ( 1.0 - beta ) * p12[1] + beta * p13[1];
  }

  return;
# undef DIM_NUM
}
