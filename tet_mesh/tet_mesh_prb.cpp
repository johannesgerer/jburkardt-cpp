# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "tet_mesh.hpp"

int main ( );
void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TET_MESH_PRB
//
//  Discussion:
//
//    TET_MESH_PRB tests the routines in TET_MESH.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TET_MESH_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TET_MESH library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TET_MESH_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST001 tests R8MAT_SOLVE.
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
{
# define N 3
# define RHS_NUM 2

  double a[N*(N+RHS_NUM)] = {
     1.0,  4.0,  7.0,
     2.0,  5.0,  8.0,
     3.0,  6.0,  0.0,
    14.0, 32.0, 23.0,
     7.0, 16.0,  7.0 };
  int i;
  int info;
  int j;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  R8MAT_SOLVE solves linear systems.\n";
//
//  Print out the matrix to be inverted.
//
  r8mat_print ( N, N+RHS_NUM, a, "  The linear system:" );
//
//  Solve the systems.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "  The input matrix was singular.\n";
    cout << "  The solutions could not be computed.\n";
    return;
  }

  cout << "\n";
  cout << "  The computed solutions:\n";
  cout << "\n";
  for ( i = 0; i < N; i++ )
  {
    for ( j = N; j < N+RHS_NUM; j++ )
    {
      cout << setw(10) << a[i+j*N] << "  ";
    }
    cout << "\n";
  }

  return;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST002 tests TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE, TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 10

  int i;
  int j;
  double phy[3*N];
  double ref[3*N];
  double ref2[3*N];
  int seed;
  double t[3*4] = {
    5.0, 0.0, 0.0,
    8.0, 0.0, 0.0,
    5.0, 2.0, 0.0,
    6.0, 1.0, 2.0 };

  seed = 123456789;

  cout << "\n";
  cout << "TEST002\n";
  cout << "  For an order 4 tetrahedron,\n";
  cout << "  TETRAHEDRON_ORDER4_PHYSICAL_TO_REFERENCE\n";
  cout << "    maps a physical point to a reference point.\n";
  cout << "  TETRAHEDRON_ORDER4_REFERENCE_TO_PHYSICAL \n";
  cout << "    maps a reference point to a physical point.\n";
  cout << "\n";
  cout << "     ( R, S, T )          ==>  ( X, Y, Z )           ==> ( R2, S2, T2 )\n";
  cout << "\n";

  tetrahedron_reference_sample ( N, &seed, ref );

  tetrahedron_order4_reference_to_physical ( t, N, ref, phy );
  tetrahedron_order4_physical_to_reference ( t, N, phy, ref2 );

  for ( j = 0; j < N; j++ )
  {
   cout << "  " << setprecision(4) << setw(8) << ref[0+j*3]
        << "  " << setprecision(4) << setw(8) << ref[1+j*3]
        << "  " << setprecision(4) << setw(8) << ref[2+j*3]
        << "  "
        << "  " << setprecision(4) << setw(8) << phy[0+j*3]
        << "  " << setprecision(4) << setw(8) << phy[1+j*3]
        << "  " << setprecision(4) << setw(8) << phy[2+j*3]
        << "  "
        << "  " << setprecision(4) << setw(8) << ref2[0+j*3]
        << "  " << setprecision(4) << setw(8) << ref2[1+j*3]
        << "  " << setprecision(4) << setw(8) << ref2[2+j*3] << "\n";
  }

  return;
# undef N
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests TETRAHEDRON_ORDER10_TO_ORDER4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int node_num1;
  int node_num2;
  double *node_xyz;
  int *tet_node1;
  int *tet_node2;
  int tet_num1;
  int tet_num2;
  int tet_order1 = 10;
  int tet_order2 = 4;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  For an order 10 tet mesh,\n";
  cout << "  TETRAHEDRON_ORDER10_TO_ORDER4\n";
  cout << "    makes a linear (order 4) tet mesh by using\n";
  cout << "    the existing nodes, and replacing each\n";
  cout << "    quadratic tetrahedron by 8 linear tetrahedrons.\n";

  tet_mesh_order10_example_size ( &node_num1, &tet_num1 );

  node_xyz = new double[3*node_num1];
  tet_node1 = new int[tet_order1*tet_num1];

  tet_mesh_order10_example_set ( node_num1, tet_num1,
    node_xyz, tet_node1 );

  i4mat_transpose_print_some ( tet_order1, tet_num1, tet_node1,
    1, 1, tet_order1, 5, "  First 5 quadratic tetrahedrons:" );

  tet_mesh_order10_to_order4_size ( node_num1, tet_num1,
    &node_num2, &tet_num2 );

  cout << "\n";
  cout << "  Quadratic mesh size is       " << tet_num1 << "\n";
  cout << "  Linearized mesh size will be " << tet_num2 << "\n";

  tet_node2 = new int[tet_order2*tet_num2];

  tet_mesh_order10_to_order4_compute ( tet_num1, tet_node1,
    tet_num2, tet_node2 );

  i4mat_transpose_print_some ( tet_order2, tet_num2, tet_node2,
    1, 1, tet_order2, 5, "  First 5 linear tetrahedrons:" );

  delete [] node_xyz;
  delete [] tet_node1;
  delete [] tet_node2;

  return;
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests TETRAHEDRON_ORDER10_TO_ORDER4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
{
  int node_num;
  int *node_order;
  double *node_xyz;
  int *tet_node;
  int tet_num;
  int tet_order = 10;

  cout << "\n";
  cout << "TEST004\n";
  cout << "  TET_MESH_NODE_ORDER determines the order of \n";
  cout << "  each node in a tet mesh.\n";
  cout << "\n";
  cout << "  The order of a node is the number of tetrahedrons\n";
  cout << "  that use the node as part of their definition.\n";

  tet_mesh_order10_example_size ( &node_num, &tet_num );

  cout << "\n";
  cout << "  This mesh has tetrahedron order " << tet_order << "\n";
  cout << "  The number of tetrahedrons is   " << tet_num << "\n";

  node_xyz = new double[3*node_num];
  tet_node = new int[tet_order*tet_num];

  tet_mesh_order10_example_set ( node_num, tet_num,
    node_xyz, tet_node );

  i4mat_transpose_print ( tet_order, tet_num, tet_node,
    "  The tet mesh:" );

  node_order = tet_mesh_node_order ( tet_order, tet_num, tet_node, node_num );

  i4vec_print ( node_num, node_order, "  Node orders:" );

  cout << "\n";
  cout << "  Check that the following are equal:\n";
  cout << "\n";
  cout << "  Number of tetrahedrons * order = " << tet_num * tet_order << "\n";
  cout << "  Sum of node orders             = " << i4vec_sum ( node_num, node_order ) << "\n";

  delete [] node_order;
  delete [] node_xyz;
  delete [] tet_node;

  return;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests TETRAHEDRON_BARYCENTRIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  double *c1;
  double c1_sum;
  double *c2;
  int i;
  double *p;
  int seed;
  int test1;
  int test1_num = 3;
  int test2;
  int test2_num = 5;
  double *tet_xyz;

  seed = 123456789;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  TETRAHEDRON_BARYCENTRIC computes the barycentric\n";
  cout << "  coordinates of a point.\n";
//
//  Choose a random tetrahedron.
//
  for ( test1 = 1; test1 <= test1_num; test1++ )
  {
    tet_xyz = r8mat_uniform_01 ( 3, 4, &seed );

    r8mat_transpose_print ( 3, 4, tet_xyz, "  Random tetrahedron:" );
//
//  Choose barycentric coordinates C1 at random.
//
//  Define a point P.
//
//  Have TETRAHEDRON_BARYCENTRIC compute C2, the barycentric coordinates of P.
//
    for ( test2 = 1; test2 <= test2_num; test2++ )
    {
      c1 = r8vec_uniform_01 ( 4, &seed );
      c1_sum = r8vec_sum ( 4, c1 );
      for ( i = 0; i < 4; i++ )
      {
        c1[i] = c1[i] / c1_sum;
      }

      p = r8mat_mv ( 3, 4, tet_xyz, c1 );

      c2 = tetrahedron_barycentric ( tet_xyz, p );

      cout << "\n";
      cout << "  C1 = ";
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << setprecision(6) << setw(14) << c1[i];
      }
      cout << "\n";
      cout << "  C2 = ";
      for ( i = 0; i < 4; i++ )
      {
        cout << "  " << setprecision(6) << setw(14) << c2[i];
      }
      cout << "\n";

      delete [] c1;
      delete [] c2;
      delete [] p;
    }
    delete tet_xyz;
  }

  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests TET_MESH_TET_NEIGHBORS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int node_num;
  double *node_xyz;
  int *tet_neighbor;
  int *tet_node;
  int tet_num;
  int tet_order = 4;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  TET_MESH_TET_NEIGHBORS computes the 4 neighboring\n";
  cout << "  tetrahedrons of each tetrahedron in a tet mesh.\n";
  cout << "  containing a point.\n";
//
//  Set up the example tetrahedron mesh.
//
  tet_mesh_order4_example_size ( &node_num, &tet_num );

  cout << "\n";
  cout << "  This mesh has tetrahedron order " << tet_order << "\n";
  cout << "  The number of tetrahedrons is   " << tet_num << "\n";

  node_xyz = new double[3*node_num];
  tet_node = new int[tet_order*tet_num];

  tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node );
//
//  Print the tets.
//
  i4mat_transpose_print_some ( tet_order, tet_num, tet_node,
    1, 1, tet_order, 10, "  First 10 Tets:" );
//
//  The TET_NEIGHBOR array is needed by TET_MESH_DELAUNAY_SEARCH.
//
  tet_neighbor = tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node );

  i4mat_transpose_print_some ( 4, tet_num, tet_neighbor,
    1, 1, 4, 10, "  First 10 Tet Neighbors:" );

  delete [] node_xyz;
  delete [] tet_neighbor;
  delete [] tet_node;

  return;
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests TET_MESH_SEARCH_NAIVE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2009
//
//  Author:
//
//    John Burkardt
//
{
  int face;
  int i;
  int j;
  int k;
  int node_num;
  double *node_xyz;
  double p[3];
  int seed;
  int step_num;
  int test;
  int test_num = 5;
  int *tet_neighbor;
  int *tet_node;
  int tet_num;
  int tet_order = 4;
  double tet_xyz[3*4];
  int tet1;
  int tet2;
  int tet3;

  seed = 123456789;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  TET_MESH_SEARCH_NAIVE uses a naive algorithm\n";
  cout << "  to search a tetrahedral mesh for the tetrahedron\n";
  cout << "  containing a point.\n";
//
//  Set up the example tetrahedron mesh.
//
  tet_mesh_order4_example_size ( &node_num, &tet_num );

  cout << "\n";
  cout << "  This mesh has tetrahedron order " << tet_order << "\n";
  cout << "  The number of tetrahedrons is   " << tet_num << "\n";

  node_xyz = new double[3*node_num];
  tet_node = new int[tet_order*tet_num];

  tet_mesh_order4_example_set ( node_num, tet_num, node_xyz, tet_node );
//
//  The TET_NEIGHBOR array is needed for the Delaunay search.
//
  tet_neighbor = tet_mesh_neighbor_tets ( tet_order, tet_num, tet_node );

  for ( test = 1; test <= test_num; test++ )
  {
//
//  Choose a tetrahedron at random.
//
    tet1 = i4_uniform ( 0, tet_num - 1, &seed );

    cout << "\n";
    cout << "  Point was chosen from tetrahedron    " << setw(8) << tet1 << "\n";

    for ( j = 0; j < 4; j++ )
    {
      k = tet_node[j+tet1*4];
      for ( i = 0; i < 3; i++ )
      {
        tet_xyz[i+j*3] = node_xyz[i+k*3];
      }
    }
//
//  Choose a point in the tetrahedron at random.
//
    tetrahedron_sample ( tet_xyz, 1, &seed, p );
//
//  Naive search.
//
    tet2 = tet_mesh_search_naive ( node_num, node_xyz, tet_order, tet_num,
      tet_node, p, &step_num );

    cout << "  Naive search ended in tetrahedron    " << setw(8) << tet2
         << ", number of steps = " << step_num << "\n";
//
//  Delaunay search.
//
    tet3 = tet_mesh_search_delaunay ( node_num, node_xyz, tet_order,
      tet_num, tet_node, tet_neighbor, p, &face, &step_num );

    cout << "  Delaunay search ended in tetrahedron " << setw(8) << tet3
         << ", number of steps = " << step_num << "\n";
  }

  delete [] node_xyz;
  delete [] tet_neighbor;
  delete [] tet_node;

  return;
}
