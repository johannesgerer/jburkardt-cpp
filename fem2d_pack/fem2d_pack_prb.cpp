# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem2d_pack.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test105 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM2D_PACK_PRB calls the various FEM2D_PACK tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  timestamp ( );

  cout << "\n";
  cout << "FEM2D_PACK_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the routines in the FEM2D_PACK library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test105 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );

  cout << "\n";
  cout << "FEM2D_PACK_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BANDWIDTH_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int *element_node;
  int element_num;
  int element_order;
  int m;
  int ml;
  int mu;
  int nelemx;
  int nelemy;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  BANDWIDTH_MESH computes the geometric bandwidth:\n";
  cout << "  of a finite element mesh.\n";

  nelemx = 2;
  nelemy = 6;

  cout << "\n";
  cout << "  NELEMX = " << nelemx << "\n";
  cout << "  NELEMY = " << nelemy << "\n";

  element_order = 6;
  element_num = grid_element_num ( "T6", nelemx, nelemy );

  cout << "\n";
  cout << "  ELEMENT_ORDER = " << element_order << "\n";
  cout << "  ELEMENT_NUM   = " << element_num << "\n";

  element_node = grid_t6_element ( nelemx, nelemy );

  grid_print ( element_order, element_num, element_node );

  bandwidth_mesh ( element_order, element_num, element_node, &ml, &mu, &m );

  cout << "\n";
  cout << "  Lower bandwidth ML = " << ml << "\n";
  cout << "  Upper bandwidth MU = " << mu << "\n";
  cout << "  Total bandwidth M  = " << m  << "\n";

  delete [] element_node;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests BANDWIDTH_VAR, NS_T6_VAR_COUNT, NS_T6_VAR_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int *element_node;
  int element_num;
  int element_order;
  int i;
  int ihi;
  int ilo;
  int m;
  int ml;
  int mu;
  int nelemx;
  int nelemy;
  int node;
  int node_num;
  int *var;
  int *var_node;
  int var_num;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  For the Navier Stokes variables associated with\n";
  cout << "  a mesh of T6 elements,\n";
  cout << "  NS_T6_VAR_COUNT counts variables,\n";
  cout << "  NS_T6_VAR_SET sets them,\n";
  cout << "  BANDWIDTH_VAR computes the variable bandwidth.\n";

  nelemx = 2;
  nelemy = 6;

  cout << "\n";
  cout << "  NELEMX = " << nelemx << "\n";
  cout << "  NELEMY = " << nelemy << "\n";

  element_order = 6;
  element_num = grid_element_num ( "T6", nelemx, nelemy );
  node_num = grid_node_num ( "T6", nelemx, nelemy );

  cout << "\n";
  cout << "  ELEMENT_ORDER = " << element_order << "\n";
  cout << "  ELEMENT_NUM   = " << element_num << "\n";
  cout << "  NODE_NUM      = " << node_num << "\n";

  element_node = grid_t6_element ( nelemx, nelemy );

  grid_print ( element_order, element_num, element_node );

  var_node = new int[node_num+1];

  var_num = ns_t6_var_count ( element_num, element_node, node_num, var_node );

  cout << "\n";
  cout << "  Number of variables VAR_NUM = " << var_num << "\n";

  i4vec_print ( node_num+1, var_node, "  VAR_NODE pointer vector:" );

  var = ns_t6_var_set ( element_num, element_node, node_num, var_node, 
    var_num );

  cout << "\n";
  cout << "  Node    U_Var    V_Var    P_Var\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    ilo = var_node[node];
    ihi = var_node[node+1]-1;

    cout << "  " << setw(8) << node;
    cout << "  ";

    for ( i = ilo; i <= ihi; i++ )
    {
      cout << "  " << setw(8) << var[i-1];
    }
    cout << "\n";
  }

  bandwidth_var ( element_order, element_num, element_node, 
    node_num, var_node, var_num, var, &ml, &mu, &m );

  cout << "\n";
  cout << "  Lower bandwidth ML = " << ml << "\n";
  cout << "  Upper bandwidth MU = " << mu << "\n";
  cout << "  Total bandwidth M  = " << m  << "\n";

  delete [] element_node;
  delete [] var;
  delete [] var_node;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests BASIS_11_**_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST03\n";
  cout << "  BASIS_11_T3_TEST - Test the T3 basis functions.\n";
  cout << "  BASIS_11_T4_TEST - Test the T4 basis functions.\n";
  cout << "  BASIS_11_T6_TEST - Test the T6 basis functions.\n";

  basis_11_t3_test ( );

  basis_11_t4_test ( );

  basis_11_t6_test ( );

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests BASIS_MN_**_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the computation of basis functions by evaluating them\n";
  cout << "  at the nodes that define the basis functions.\n";
  cout << "\n";
  cout << "  BASIS_MN_Q4_TEST - for the Q4 element.\n";
  cout << "  BASIS_MN_T3_TEST - for the T3 element.\n";
  cout << "  BASIS_MN_T4_TEST - for the T4 element.\n";
  cout << "  BASIS_MN_T6_TEST - for the T6 element.\n";

  basis_mn_q4_test ( );

  basis_mn_t3_test ( );

  basis_mn_t4_test ( );

  basis_mn_t6_test ( );

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 demonstrates DERIVATIVE_AVERAGE_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 June 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_ORDER 3
# define NELEMX 7
# define NELEMY 5

  double angle;
  double *c;
  int col;
  double *dcdx;
  double dcdx_exact;
  double *dcdy;
  double dcdy_exact;
  int element;
  int *element_node;
  int element_num;
  int node;
  int node_num;
  double *node_xy;
  double r;
  int row;
  double x;
  double y;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  DERIVATIVE_AVERAGE_T3 averages the spatial derivatives\n";
  cout << "  of a finite element function at the nodes.\n";
//
//  How many elements are there?
//
  element_num = grid_t3_element_num ( NELEMX, NELEMY );
//
//  How many nodes are there?
//
  node_num = grid_t3_node_num ( NELEMX, NELEMY );

  c = new double[node_num];
  dcdx = new double[node_num];
  dcdy = new double[node_num];
  node_xy = new double[2*node_num];
//
//  Get the nodes that make up each element.
//
  element_node = grid_t3_element ( NELEMX, NELEMY );
//
//  Generate the coordinates of the nodes.
//
  node = 0;

  for ( row = 0; row <= NELEMY; row++ )
  {
    r = ( ( double ) ( NELEMY - row ) *  1.0   
        + ( double ) (        + row ) *  3.0 ) 
        / ( double ) ( NELEMY       );

    for ( col = 0; col <= NELEMX; col++ )
    {
      angle = ( ( double ) ( NELEMX - col ) * 135.0   
              + ( double ) (        + col ) *  45.0 ) 
              / ( double ) ( NELEMX       );

      angle = degrees_to_radians ( angle );

      node_xy[0+node*2] = r * cos ( angle );
      node_xy[1+node*2] = r * sin ( angle );

      node = node + 1;
    }
  }
//
//  Set the finite element function.
//
  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];
    c[node] = sin ( x ) * ( 1.0 + y * y );
  }

  derivative_average_t3 ( node_num, node_xy, element_num, 
    element_node, c, dcdx, dcdy );

  cout << "\n";
  cout << "  C         X               Y\n";
  cout << "         dCdX(computed)  dCdY(computed)\n";
  cout << "         dCdX(exact)     dCdY(exact)\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    x = node_xy[0+node*2];
    y = node_xy[1+node*2];

    dcdx_exact = cos ( x ) * ( 1.0 * y * y );
    dcdy_exact = sin ( x ) * 2.0 * y;

    cout << "\n";
    cout << "  " << setw(14) << c[node]
         << "  " << setw(14) << node_xy[0+node*2] 
         << "  " << setw(14) << node_xy[1+node*2] << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << dcdx[node] 
         << "  " << setw(14) << dcdy[node] << "\n";
    cout << "  " << "              "
         << "  " << setw(14) << dcdx_exact 
         << "  " << setw(14) << dcdy_exact << "\n";
  }

  delete [] c;
  delete [] dcdx;
  delete [] dcdy;
  delete [] element_node;
  delete [] node_xy;

  return;
# undef ELEMENT_ORDER
# undef NELEMX
# undef NELEMY
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests ELEMENT_EPS using Q4 elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_ORDER 4
# define NELEMX 7
# define NELEMY 5

# define ELEMENT_NUM ( NELEMX * NELEMY )
# define NODE_NUM ( NELEMX + 1 ) * ( NELEMY + 1 )

  double angle;
  char code[3] = "Q4";
  int col;
  int element;
  bool element_mask[ELEMENT_NUM];
  int *element_node;
  int element_show = 2;
  char file_name[] = "fem2d_pack_prb_q4.eps";
  int node;
  int node_show = 2;
  double node_xy[2*NODE_NUM];
  double r;
  int row;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  ELEMENTS_EPS creates an Encapsulated PostScript\n";
  cout << "  file containing an image of a mesh.\n";

  element_node = grid_q4_element ( NELEMX, NELEMY );

  node = 0;

  for ( row = 0; row <= NELEMY; row++ )
  {
    r = ( ( double ) ( NELEMY - row ) *  1.0   
        + ( double ) (        + row ) *  3.0 ) 
        / ( double ) ( NELEMY       );

    for ( col = 0; col <= NELEMX; col++ )
    {
      angle = ( ( double ) ( NELEMX - col ) * 135.0   
              + ( double ) (        + col ) *  45.0 ) 
              / ( double ) ( NELEMX       );

      angle = degrees_to_radians ( angle );

      node_xy[0+node*2] = r * cos ( angle );
      node_xy[1+node*2] = r * sin ( angle );

      node = node + 1;
    }

  }

  for ( element = 0; element < ELEMENT_NUM; element++ )
  {
    element_mask[element] = true;
  }

  elements_eps ( file_name, NODE_NUM, node_xy, code, 
    ELEMENT_NUM, element_mask, element_node, node_show, element_show );

  delete [] element_node;

  return;
# undef ELEMENT_ORDER
# undef NELEMX
# undef NELEMY
# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests ELEMENT_EPS, using T3 elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 September 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_ORDER 3
# define NELEMX 7
# define NELEMY 5

# define ELEMENT_NUM ( 2 * NELEMX * NELEMY )
# define NODE_NUM ( NELEMX + 1 ) * ( NELEMY + 1 )

  double angle;
  char code[3] = "T3";
  int col;
  int element;
  bool element_mask[ELEMENT_NUM];
  int *element_node;
  int element_show = 2;
  char file_name[] = "fem2d_pack_prb_t3.eps";
  int node;
  int node_show = 2;
  double node_xy[2*NODE_NUM];
  double r;
  int row;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  ELEMENTS_EPS creates an Encapsulated PostScript\n";
  cout << "  file containing an image of a mesh.\n";

  element_node = grid_t3_element ( NELEMX, NELEMY );

  node = 0;

  for ( row = 0; row <= NELEMY; row++ )
  {
    r = ( ( double ) ( NELEMY - row ) *  1.0   
        + ( double ) (        + row ) *  3.0 ) 
        / ( double ) ( NELEMY       );

    for ( col = 0; col <= NELEMX; col++ )
    {
      angle = ( ( double ) ( NELEMX - col ) * 135.0   
              + ( double ) (        + col ) *  45.0 ) 
              / ( double ) ( NELEMX       );

      angle = degrees_to_radians ( angle );

      node_xy[0+node*2] = r * cos ( angle );
      node_xy[1+node*2] = r * sin ( angle );

      node = node + 1;
    }
  }

  for ( element = 0; element < ELEMENT_NUM; element++ )
  {
    element_mask[element] = true;
  }

  elements_eps ( file_name, NODE_NUM, node_xy, code, 
    ELEMENT_NUM, element_mask, element_node, node_show, element_show );

  delete [] element_node;

  return;
# undef ELEMENT_ORDER
# undef NELEMX
# undef NELEMY
# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests ELEMENT_EPS, using T4 elements.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 March 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_ORDER 4
# define NELEMX 7
# define NELEMY 5

  char code[3] = "T4";
  int col;
  int element;
  bool *element_mask;
  int *element_node;
  int element_num;
  int element_show = 1;
  char file_name[] = "fem2d_pack_prb_t4.eps";
  int i;
  int node;
  int node_num;
  int node_show = 2;
  double *node_xy;
  int row;
  double x;
  double y;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  ELEMENTS_EPS creates an Encapsulated PostScript\n";
  cout << "  file containing an image of a T4 mesh.\n";
//
//  How many elements are there?
//
  element_num = grid_t4_element_num ( NELEMX, NELEMY );

  element_mask = new bool[element_num];
//
//  How many nodes are there?
//
  node_num = grid_t4_node_num ( NELEMX, NELEMY );

  node_xy = new double[2*node_num];
//
//  Get the nodes that make up each element.
//
  element_node = grid_t4_element ( NELEMX, NELEMY );

  node = 0;

  for ( row = 0; row <= NELEMY; row++ )
  {
    y = ( ( double ) ( 3*NELEMY - row ) *  0.0   
        + ( double ) (          + row ) *  6.0 ) 
        / ( double ) ( 3*NELEMY       );

    for ( col = 0; col <= NELEMX; col++ )
    {
      x = ( ( double ) ( 2*NELEMX - col ) * 0.0   
          + ( double ) (          + col ) * 6.0 ) 
          / ( double ) ( 2*NELEMX       );

      node_xy[0+node*2] = x;
      node_xy[1+node*2] = y;

      node = node + 1;
    }
//
//  Skip over the two rows of interior nodes.
//
    node = node + NELEMX;
    node = node + NELEMX;
  }
//
//  The coordinates of interior nodes are the average of the vertices.
//
  for ( element = 0; element < element_num; element++ )
  {
    x = 0.0;
    y = 0.0;
    for ( i = 0; i < 3; i++ )
    {
      node = element_node[i+element*4];
      x = x + node_xy[0+(node-1)*2];
      y = y + node_xy[1+(node-1)*2];
    }
    node = element_node[3+element*4];
    node_xy[0+(node-1)*2] = x / 3.0;
    node_xy[1+(node-1)*2] = y / 3.0;
  }

  for ( element = 0; element < element_num; element++ )
  {
    element_mask[element] = true;
  }

  elements_eps ( file_name, node_num, node_xy, code, 
    element_num, element_mask, element_node, node_show, element_show );

  delete [] element_mask;
  delete [] element_node;
  delete [] node_xy;

  return;
# undef ELEMENT_ORDER
# undef NELEMX
# undef NELEMY
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests ELEMENT_EPS, using T6 elements.
//
//  Discussion:
//
//    We generate a big grid of T6 elements, but we only want to
//    look at the six elements shared by node 85.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_ORDER 6
# define NELEMX 6
# define NELEMY 6

# define ELEMENT_NUM ( 2 * NELEMX * NELEMY )
# define NODE_NUM ( 2 * NELEMX + 1 ) * ( 2 * NELEMY + 1 )

  char code[3] = "T6";
  int col;
  int element;
  bool element_mask[ELEMENT_NUM];
  int *element_node;
  int element_show = 2;
  char file_name[] = "fem2d_pack_prb_t6.eps";
  int node;
  int node_show = 2;
  double node_xy[2*NODE_NUM];
  double x;
  double y;
  int row;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  ELEMENTS_EPS creates an Encapsulated PostScript\n";
  cout << "  file containing an image of a mesh.\n";

  element_node = grid_t6_element ( NELEMX, NELEMY );

  node = 0;

  for ( row = 0; row <= 2 * NELEMY; row++ )
  {
    y = ( ( double ) ( 2 * NELEMY - row ) *  0.0   
        + ( double ) (            + row ) *  6.0 ) 
        / ( double ) ( 2 * NELEMY       );

    for ( col = 0; col <= 2 * NELEMX; col++ )
    {
      x = ( ( double ) ( 2 * NELEMX - col ) * 0.0   
          + ( double ) (            + col ) * 6.0 ) 
          / ( double ) ( 2 * NELEMX       );

      node_xy[0+node*2] = x;
      node_xy[1+node*2] = y;

      node = node + 1;
    }

  }

  for ( element = 0; element < ELEMENT_NUM; element++ )
  {
    element_mask[element] = false;
  }
  element_mask[30-1] = true;
  element_mask[31-1] = true;
  element_mask[32-1] = true;
  element_mask[41-1] = true;
  element_mask[42-1] = true;
  element_mask[43-1] = true;

  elements_eps ( file_name, NODE_NUM, node_xy, code, 
    ELEMENT_NUM, element_mask, element_node, node_show, element_show );

  delete [] element_node;

  return;
# undef ELEMENT_ORDER
# undef NELEMX
# undef NELEMY

# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests GRID_NODES_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 May 2008
//
//  Author:
//
//    John Burkardt
//
{
  int node;
  int node_num;
  double *node_xy;
  int num_x = 5;
  int num_y = 3;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  GRID_NODES_01 computes a regular grid in the unit square.\n";
  cout << "\n";
  cout << "  NUM_X =    " << num_x << "\n";
  cout << "  NUM_Y =    " << num_y << "\n";
  node_num = num_x * num_y;
  cout << "  NODE_NUM = " << node_num << "\n";
  cout << "\n";

  node_xy = grid_nodes_01 ( num_x, num_y );

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(8)  << node
         << "  " << setw(14) << node_xy[0+node*2]
         << "  " << setw(14) << node_xy[1+node*2] << "\n";
  }

  delete [] node_xy;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests GRID_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST11\n";
  cout << "  GRID_TEST tests the grid routines.\n";

  grid_test ( "Q4" );

  grid_test ( "Q8" );

  grid_test ( "Q9" );

  grid_test ( "Q12" );

  grid_test ( "Q16" );

  grid_test ( "QL" );

  grid_test ( "T3" );

  grid_test ( "T6" );

  grid_test ( "T10" );

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests INTERP_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST12\n";
  cout << "  INTERP_TEST tests the interpolating power\n";
  cout << "  of the element.\n";

  interp_test ( "Q4" );

  interp_test ( "Q8" );

  interp_test ( "Q9" );

  interp_test ( "Q12" );

  interp_test ( "Q16" );

  interp_test ( "QL" );

  interp_test ( "T3" );

  interp_test ( "T4" );

  interp_test ( "T6" );

  interp_test ( "T10" );

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests MAP_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST13\n";
  cout << "  MAP_TEST tests the map routines.\n";

  map_test ( "Q4" );

  map_test ( "Q8" );

  map_test ( "Q9" );

  map_test ( "Q12" );

  map_test ( "Q16" );

  map_test ( "QL" );

  map_test ( "T3" );

  map_test ( "T6" );

  map_test ( "T10" );

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests MASS_MATRIX_T6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ELEMENT_NUM 2
# define NODE_NUM 9

  double *a;
  int element_node[6*ELEMENT_NUM] = {
    1, 3, 7, 2, 5, 4,
    9, 7, 3, 8, 5, 6 };
  double node_xy[2*NODE_NUM] = {
    0.0, 0.0,
    0.0, 0.5,
    0.0, 1.0,
    0.5, 0.0,
    0.5, 0.5,
    0.5, 1.0,
    1.0, 0.0,
    1.0, 0.5,
    1.0, 1.0 };

  cout << "\n";
  cout << "TEST14\n";
  cout << "  MASS_MATRIX_T6 computes the mass matrix for\n";
  cout << "  a finite element system using T6 elements\n";
  cout << "  (quadratic triangles).\n";

  a = mass_matrix_t6 ( NODE_NUM, ELEMENT_NUM, element_node, node_xy );

  r8mat_print ( NODE_NUM, NODE_NUM, a, "  The T6 mass matrix:" );

  delete [] a;

  return;
# undef ELEMENT_NUM
# undef NODE_NUM
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests PHYSICAL_TO_REFERENCE_T3 and REFERENCE_TO_PHYSICAL_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int j;
  double phy[2*N];
  double ref[2*N];
  double ref2[2*N];
  int seed;
  double t[2*3] = {
    1.0, 1.0, 
    3.0, 1.0, 
    2.0, 5.0 };

  seed = 123456789;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  For an order 3 triangle,\n";
  cout << "  PHYSICAL_TO_REFERENCE_T3 maps a physical point to\n";
  cout << "    a reference point.\n";
  cout << "  REFERENCE_TO_PHYSICAL_T3 maps a reference point to\n";
  cout << "    a physical point.\n";
  cout << "\n";
  cout << "      XSI     ETA  ==>  X       Y    ==>  XSI2    ETA2\n";
  cout << "\n";

  for ( j = 0; j < N; j++ )
  {
    for ( i = 0; i < 2; i++ )
    {
      ref[i+j*2] = r8_uniform_01 ( &seed );
    }

    if ( 1.0 < ref[0+j*2] + ref[1+j*2] )
    {
      ref[0+j*2] = 1.0 - ref[0+j*2];
      ref[1+j*2] = 1.0 - ref[1+j*2];
    }
  }

  reference_to_physical_t3 ( t, N, ref, phy );

  physical_to_reference_t3 ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(10) << ref[0+j*2]
         << "  " << setw(10) << ref[1+j*2]
         << "  " << setw(10) << phy[0+j*2]
         << "  " << setw(10) << phy[1+j*2]
         << "  " << setw(10) << ref2[0+j*2]
         << "  " << setw(10) << ref2[1+j*2] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests REFERENCE_TO_PHYSICAL_T6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 16

  int i;
  int j;
  double phy[2*N];
  double ref[2*N] = {
    0.00, 0.00, 
    1.00, 0.00, 
    0.00, 1.00, 
    0.50, 0.00, 
    0.50, 0.50, 
    0.00, 0.50, 
    0.25, 0.75, 
    0.75, 0.25, 
    0.40, 0.10, 
    0.30, 0.20, 
    0.20, 0.30, 
    0.10, 0.40, 
    0.10, 0.10, 
    0.20, 0.20, 
    0.30, 0.30, 
    0.40, 0.40 };
  double t[2*6] = {
    0.0, 0.0, 
    2.0, 0.0, 
    0.0, 4.0, 
    1.0, 0.0, 
    1.0, 1.0, 
    0.0, 2.0 };

  cout << "\n";
  cout << "TEST16\n";
  cout << "  For an order 6 triangle,\n";
  cout << "  REFERENCE_TO_PHYSICAL_T6 maps a reference point to\n";
  cout << "    a physical point.\n";
  cout << "\n";
  cout << "      XSI     ETA  ==>  X       Y\n";
  cout << "\n";

  reference_to_physical_t6 ( t, N, ref, phy );

  for ( j = 0; j < N; j++ )
  {
    cout << "  " << setw(8) << ref[0+j*2] 
         << "  " << setw(8) << ref[1+j*2] 
         << "  " << setw(8) << phy[0+j*2] 
         << "  " << setw(8) << phy[1+j*2] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests the shape routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST18\n";
  cout << "  SHAPE_TEST tests the shape routines.\n";

  shape_test ( "Q4" );

  shape_test ( "Q8" );

  shape_test ( "Q9" );

  shape_test ( "Q12" );

  shape_test ( "Q16" );

  shape_test ( "QL" );

  shape_test ( "T3" );

  shape_test ( "T6" );

  shape_test ( "T10" );

  return;
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests SPHERE_GRID_Q4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  int *element_node;
  int element_num;
  int element_order = 4;
  int nelemx = 8;
  int nelemy = 8;
  int node;
  int node_num;
  double *node_xyz;
  int order;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  SPHERE_GRID_Q4_ELEMENT sets up a grid of\n";
  cout << "    Q4 quadrilaterals on a sphere.\n";
  cout << "  SPHERE_GRID_Q4_ELEMENT_NUM returns the number\n";
  cout << "    of elements in the grid\n";
  cout << "  SPHERE_GRID_Q4_NODE_NUM returns the number\n";
  cout << "    of nodes in the grid.\n";
  cout << "  SPHERE_GRID_Q4_NODE_XYZ returns the coordinates\n";
  cout << "    of nodes in the grid.\n";

  element_num = sphere_grid_q4_element_num ( nelemx, nelemy );
  node_num = sphere_grid_q4_node_num ( nelemx, nelemy );

  cout << "\n";
  cout << "  Expected number of nodes =    " << node_num    << "\n";
  cout << "  Expected number of elements = " << element_num << "\n";

  element_node = sphere_grid_q4_element ( nelemx, nelemy );

  cout << "\n";
  cout << "  The elements and their nodes:\n";
  cout << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    cout << "  " << setw(4) << element + 1 << "  ";
    for ( order = 0; order < element_order; order++ )
    {
      cout << "  " << setw(4) << element_node[order+element*element_order];
    }
    cout << "\n";
  }

  node_xyz = sphere_grid_q4_node_xyz ( nelemx, nelemy );

  cout << "\n";
  cout << "  The node coordinates:\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node + 1
         << "  " << setw(12) << node_xyz[0+node*3]
         << "  " << setw(12) << node_xyz[1+node*3]
         << "  " << setw(12) << node_xyz[2+node*3] << "\n";
  }
//
//  Write the elements and nodes to files.
//
  r8mat_write ( "sphere_q4_nodes.txt", 3, node_num, node_xyz );

  i4mat_write ( "sphere_q4_elements.txt", element_order, element_num, 
    element_node );

  delete [] element_node;
  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests SPHERE_GRID_Q9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  int *element_node;
  int element_num;
  int element_order = 9;
  int nelemx = 3;
  int nelemy = 4;
  int node;
  int node_num;
  double *node_xyz;
  int order;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  SPHERE_GRID_Q9_ELEMENT sets up a grid of\n";
  cout << "    Q9 quadrilaterals on a sphere.\n";
  cout << "  SPHERE_GRID_Q9_ELEMENT_NUM returns the number\n";
  cout << "    of elements in the grid\n";
  cout << "  SPHERE_GRID_Q9_NODE_NUM returns the number\n";
  cout << "    of nodes in the grid.\n";
  cout << "  SPHERE_GRID_Q9_NODE_XYZ returns the coordinates\n";
  cout << "    of nodes in the grid.\n";

  element_num = sphere_grid_q9_element_num ( nelemx, nelemy );
  node_num = sphere_grid_q9_node_num ( nelemx, nelemy );

  cout << "\n";
  cout << "  Expected number of nodes =    " << node_num    << "\n";
  cout << "  Expected number of elements = " << element_num << "\n";

  element_node = sphere_grid_q9_element ( nelemx, nelemy );

  cout << "\n";
  cout << "  The elements and their nodes:\n";
  cout << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    cout << "  " << setw(4) << element + 1 << "  ";
    for ( order = 0; order < element_order; order++ )
    {
      cout << "  " << setw(4) << element_node[order+element*element_order];
    }
    cout << "\n";
  }

  node_xyz = sphere_grid_q9_node_xyz ( nelemx, nelemy );

  cout << "\n";
  cout << "  The node coordinates:\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node + 1
         << "  " << setw(12) << node_xyz[0+node*3]
         << "  " << setw(12) << node_xyz[1+node*3]
         << "  " << setw(12) << node_xyz[2+node*3] << "\n";
  }
//
//  Write the elements and nodes to files.
//
  r8mat_write ( "sphere_q9_nodes.txt", 3, node_num, node_xyz );

  i4mat_write ( "sphere_q9_elements.txt", element_order, element_num, 
    element_node );

  delete [] element_node;
  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests SPHERE_GRID_Q16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  int *element_node;
  int element_num;
  int element_order = 16;
  int i;
  int ilo;
  int j;
  int k;
  int nelemx = 2;
  int nelemy = 2;
  int node;
  int node_num;
  double *node_xyz;
  int order;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  SPHERE_GRID_Q16_ELEMENT sets up a grid of\n";
  cout << "    Q16 quadrilaterals on a sphere.\n";
  cout << "  SPHERE_GRID_Q16_ELEMENT_NUM returns the number\n";
  cout << "    of elements in the grid\n";
  cout << "  SPHERE_GRID_Q16_NODE_NUM returns the number\n";
  cout << "    of nodes in the grid.\n";
  cout << "  SPHERE_GRID_Q16_NODE_XYZ returns the coordinates\n";
  cout << "    of nodes in the grid.\n";

  element_num = sphere_grid_q16_element_num ( nelemx, nelemy );
  node_num = sphere_grid_q16_node_num ( nelemx, nelemy );

  cout << "\n";
  cout << "  Expected number of nodes =    " << node_num    << "\n";
  cout << "  Expected number of elements = " << element_num << "\n";

  element_node = sphere_grid_q16_element ( nelemx, nelemy );

  cout << "\n";
  cout << "  The elements and their nodes, listed in a way\n";
  cout << "  that suggests their geometry:\n";
  cout << "\n";

  element = element_num;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      element = element - 1;
      cout << "\n";

      for ( ilo = 12; 0 <= ilo; ilo = ilo - 4 )
      {
        if ( ilo == 12 )
        {
          cout << "  " << setw(4) << element + 1;
        }
        else
        {
          cout << "      ";
        }

        for ( k = ilo; k <= ilo + 3; k++ )
        {
           cout << "  " << setw(4) << element_node[k+element*element_order];
        }
        cout << "\n";
      }
    }
  }

  node_xyz = sphere_grid_q16_node_xyz ( nelemx, nelemy );

  cout << "\n";
  cout << "  The node coordinates:\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node + 1
         << "  " << setw(12) << node_xyz[0+node*3]
         << "  " << setw(12) << node_xyz[1+node*3]
         << "  " << setw(12) << node_xyz[2+node*3] << "\n";
  }
//
//  Write the elements and nodes to files.
//
  r8mat_write ( "sphere_q16_nodes.txt", 3, node_num, node_xyz );

  i4mat_write ( "sphere_q16_elements.txt", element_order, element_num, 
    element_node );

  delete [] element_node;
  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests SPHERE_GRID_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int dim;
  int dim_num = 3;
  int element;
  char element_file_name[80] = "sphere_t3_elements.txt";
  int *element_node;
  int element_num;
  int element_order = 3;
  int nelemx = 8;
  int nelemy = 8;
  int node;
  int node_num;
  char node_file_name[80] = "sphere_t3_nodes.txt";
  double *node_xyz;
  int order;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  SPHERE_GRID_T3_ELEMENT sets up a grid of T3 triangles\n";
  cout << "    on a sphere.\n";
  cout << "  SPHERE_GRID_T3_ELEMENT_NUM returns the number\n";
  cout << "    of elements in the grid\n";
  cout << "  SPHERE_GRID_T3_NODE_NUM returns the number\n";
  cout << "    of nodes in the grid.\n";
  cout << "  SPHERE_GRID_T3_NODE_XYZ returns the coordinates\n";
  cout << "    of nodes in the grid.\n";

  element_num = sphere_grid_t3_element_num ( nelemx, nelemy );
  node_num = sphere_grid_t3_node_num ( nelemx, nelemy );

  cout << "\n";
  cout << "  Expected number of nodes =    " << node_num    << "\n";
  cout << "  Expected number of elements = " << element_num << "\n";
//
//  Generate the ELEMENT_NODE array, print it, and write it to a file.
//
  element_node = sphere_grid_t3_element ( nelemx, nelemy );

  cout << "\n";
  cout << "  The elements and their nodes:\n";
  cout << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    cout << "  " << setw(4) << element + 1 << "  ";
    for ( order = 0; order < element_order; order++ )
    {
      cout << "  " << setw(4) << element_node[order+element*element_order];
    }
    cout << "\n";
  }

  i4mat_write ( element_file_name, element_order, element_num, 
    element_node );
//
//  Generate the NODE_XYZ array, print it, and write it to a file.
//
  node_xyz = sphere_grid_t3_node_xyz ( nelemx, nelemy );

  cout << "\n";
  cout << "  The node coordinates:\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node + 1;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(12) << node_xyz[dim+node*3];
    }
    cout << "\n";
  }
//
//  Write the elements and nodes to files.
//
  r8mat_write ( node_file_name, dim_num, node_num, node_xyz );

  delete [] element_node;
  delete [] node_xyz;

  return;
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests SPHERE_GRID_T6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 September 2006
//
//  Author:
//
//    John Burkardt
//
{
  int element;
  int *element_node;
  int element_num;
  int element_order = 6;
  int nelemx = 3;
  int nelemy = 4;
  int node;
  int node_num;
  double *node_xyz;
  int order;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  SPHERE_GRID_T6_ELEMENT sets up a grid of T6 triangles\n";
  cout << "    on a sphere.\n";
  cout << "  SPHERE_GRID_T6_ELEMENT_NUM returns the number\n";
  cout << "    of elements in the grid\n";
  cout << "  SPHERE_GRID_T6_NODE_NUM returns the number\n";
  cout << "    of nodes in the grid\n";
  cout << "  SPHERE_GRID_T6_NODE_XYZ returns the coordinates\n";
  cout << "    of nodes in the grid.\n";

  element_num = sphere_grid_t6_element_num ( nelemx, nelemy );
  node_num = sphere_grid_t6_node_num ( nelemx, nelemy );

  cout << "\n";
  cout << "  Expected number of nodes =    " << node_num    << "\n";
  cout << "  Expected number of elements = " << element_num << "\n";

  element_node = sphere_grid_t6_element ( nelemx, nelemy );

  cout << "\n";
  cout << "  The elements and their nodes:\n";
  cout << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    cout << "  " << setw(4) << element + 1 << "  ";
    for ( order = 0; order < element_order; order++ )
    {
      cout << "  " << setw(4) << element_node[order+element*element_order];
    }
    cout << "\n";
  }

  node_xyz = sphere_grid_t6_node_xyz ( nelemx, nelemy );

  cout << "\n";
  cout << "  The node coordinates:\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(4) << node + 1
         << "  " << setw(12) << node_xyz[0+node*3]
         << "  " << setw(12) << node_xyz[1+node*3]
         << "  " << setw(12) << node_xyz[2+node*3] << "\n";
  }
//
//  Write the elements and nodes to files.
//
  r8mat_write ( "sphere_t6_nodes.txt", 3, node_num, node_xyz );

  i4mat_write ( "sphere_t6_elements.txt", element_order, element_num, 
    element_node );

  delete [] element_node;
  delete [] node_xyz;

  return;
}
//*****************************************************************************

void test24 ( )

//*****************************************************************************
//
//  Purpose:
//
//    TEST24 tests TRIANGLE_UNIT_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2006
//
//  Author:
//
//    John Burkardt
//
{
# define ORDER_MAX 64

  int a;
  int b;
  double coef;
  double err;
  double exact;
  int i;
  int order;
  double quad;
  int rule;
  int rule_max = 20;
  double value;
  double weight[ORDER_MAX];
  double x;
  double xtab[ORDER_MAX];
  double y;
  double ytab[ORDER_MAX];

  cout << "\n";
  cout << "TEST24\n";
  cout << "  TRIANGLE_UNIT_SET sets up a quadrature\n";
  cout << "    in the unit triangle,\n";
  cout << "\n";

  for ( a = 0; a <= 10; a++ )
  {
    for ( b = 0; b <= 10 - a; b++ )
    {
      coef = ( double ) ( a + b + 2 ) * ( double ) ( a + b + 1 );
      for ( i = 1; i <= b; i++ )
      {
        coef = coef * ( double ) ( a + i ) / ( double ) ( i );
      }

      cout << "\n";
      cout << "  A = " << a << "  B = " << b << "\n";
      cout << "\n";
      cout << "  Rule       QUAD           ERROR\n";
      cout << "\n";

      for ( rule = 1; rule <= rule_max; rule++ )
      {
        order = triangle_unit_size ( rule );

        triangle_unit_set ( rule, xtab, ytab, weight );

        quad = 0.0;

        for ( i = 0; i < order; i++ )
        {
          x = xtab[i];
          y = ytab[i];

          if ( a == 0 && b == 0 )
          {
            value = coef;
          }
          else if ( a == 0 && b != 0 )
          {
            value = coef * pow ( ytab[i], b );
          }
          else if ( a != 0 && b == 0 )
          {
            value = coef * pow ( xtab[i], a );
          }
          else if ( a != 0 && b != 0 )
          {
            value = coef * pow ( xtab[i], a ) * pow ( ytab[i], b );
          }

          quad = quad + 0.5 * weight[i] * value;

        }

        exact = 1.0;
        err = fabs ( exact - quad );

        cout << "  " << setw(4) << rule
             << "  " << setw(14) << quad 
             << "  " << setw(14) << err << "\n";
      }
    }
  }

  return;
# undef ORDER_MAX
}
