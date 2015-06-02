# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "treepack.hpp"

int main ( );
void test005 ( );
void test006 ( );
void test01 ( );
void test02 ( );
void test025 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TREEPACK_PRB.
//
//  Discussion:
//
//    TREEPACK_PRB tests the TREEPACK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TREEPACK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TREEPACK library.\n";

  test005 ( );
  test006 ( );
  test01 ( );
  test02 ( );
  test025 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TREEPACK_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests CATALAN and CATALAN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int c;
  int *c2;
  int n;
  int n_data;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  CATALAN computes Catalan numbers.\n";
  cout << "  CATALAN_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  N  exact C(I)  computed C(I)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    catalan_values ( n_data, n, c );

    if ( n_data == 0 )
    {
      break;
    }

    c2 = catalan ( n );

    cout << "  " << setw(4) << n
         << "  " << setw(6) << c
         << "  " << setw(6) << c2[n] << "\n";
    delete [] c2;

  }

  return;
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests CBT_TRAVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int depth = 4;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  CBT_TRAVERSE traverses a complete binary tree.\n";
  cout << "\n";
  cout << "  For this demonstration, we simply print our path.\n";
  cout << "  The tree depth is " << depth << "\n";
  cout << "\n";

  cbt_traverse ( depth );

  return;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests PRUEFER_TO_TREE_ARC.
//
//  Discussion:
//
//    The tree is
//
//          5
//          |
//    2-3-6-8-1-9
//      |   |
//      7   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int code[7] = { 1, 3, 8, 8, 3, 6, 8 };
  int *inode;
  int *jnode;
  int nnode = 9;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  PRUEFER_TO_TREE_ARC computes a tree from its Pruefer code.\n";
  cout << "\n";
  cout << "          5\n";
  cout << "          |\n";
  cout << "    2-3-6-8-1-9\n";
  cout << "      |   |\n";
  cout << "      7   4\n";

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );

  inode = new int[nnode-1];
  jnode = new int[nnode-1];

  pruefer_to_tree_arc ( nnode, code, inode, jnode );
 
  graph_arc_print ( nnode-1, inode, jnode, "  The graph:" );

  delete [] inode;
  delete [] jnode;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PRUEFER_TO_TREE_2_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 9

  int code[NNODE] = { 1, 3, 8, 8, 3, 6, 8, 0, 0 };
  int *itree;
  int nnode = NNODE;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PRUEFER_TO_TREE_2_NEW produces a tree from its Pruefer code\n";

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );

  itree = pruefer_to_tree_2_new ( nnode, code );
 
  i4vec_print ( nnode-1, itree, "  The edge list of the tree:" );
 
  delete [] itree;

  return;
# undef NNODE
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests PRUEFER_TO_TREE_2.
//
//  Discussion:
//
//    This example is used to illustrate the Nijenhuis and Wilf algorithm
//    LBLTRE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 4

  int code[NNODE];
  int i;
  int itree[NNODE];
  int j;
  int nnode = NNODE;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  PRUEFER_TO_TREE_2 produces a tree from its Pruefer code\n";
  cout << "\n";
  cout << "   Code      Tree\n";
  cout << "\n";
  for ( j = 1; j <= nnode; j++ )
  {
    code[1] = j;
    for ( i = 1; i <= nnode; i++ )
    {
      code[0] = i;
      pruefer_to_tree_2 ( nnode, code, itree );
      cout << "  " << setw(2) << code[0] 
           << "  " << setw(2) << code[1]
           << "  " 
           << "  " << setw(2) << itree[0]
           << "  " << setw(2) << itree[1]
           << "  " << setw(2) << itree[2] << "\n";
    }
  }

  return;
# undef NNODE
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests TREE_ARC_TO_PRUEFER.
//
//  Discussion:
//
//    The tree is
//
//          5
//          |
//    2-3-6-8-1-9
//      |   |
//      7   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int *code;
  int inode[8] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[8] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int nnode = 9;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  TREE_ARC_TO_PRUEFER: Tree => Pruefer code\n";
  cout << "\n";
  cout << "          5\n";
  cout << "          |\n";
  cout << "    2-3-6-8-1-9\n";
  cout << "      |   |\n";
  cout << "      7   4\n";

  graph_arc_print ( nnode-1, inode, jnode, "  The graph:" );
 
  code = tree_arc_to_pruefer ( nnode, inode, jnode );

  i4vec_print ( nnode-2, code, "  The Pruefer code:" );
 
  delete [] code;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests TREE_ARC_CENTER.
//
//  Discussion:
//
//    The tree is
//
//    2---3---6---8---1---9
//       /       / \
//      7       5   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 9

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[NNODE-1] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int nnode = NNODE;
  int parity;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  TREE_ARC_CENTER computes the center of a tree.\n";

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, eccent, parity );

  cout << "\n";
  cout << "  Parity = " << parity << "\n";
  cout << "  Eccentricity is " << eccent << "\n";

  if ( parity == 0 )
  {
    cout << "  No center node (degenerate case).\n";
  }
  else if ( parity == 1 )
  {
    cout << "  Center node: " << center[0] << "\n";
  }
  else
  {
    cout << "  Center nodes: " << center[0] << "  " << center[1] << "\n";
  }

  return;
# undef NNODE
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests TREE_ARC_CENTER.
//
//  Discussion:
//
//    Compare:
//
//    2--1--3
//
//    1--2--3
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 3

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1];
  int jnode[NNODE-1];
  int nnode = NNODE;
  int parity;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  TREE_ARC_CENTER computes the center of a tree.\n";

  inode[0] = 1;
  inode[1] = 1;
  jnode[0] = 2;
  jnode[1] = 3;

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, eccent, parity );

  cout << "\n";
  cout << "  Parity = " << parity << "\n";
  cout << "  Eccentricity is " << eccent << "\n";

  if ( parity == 0 )
  {
    cout << "  No center node (degenerate case).\n";
  }
  else if ( parity == 1 )
  {
    cout << "  Center node: " << center[0] << "\n";
  }
  else
  {
    cout << "  Center nodes: " << center[0] << "  " << center[1] << "\n";
  }

  inode[0] = 2;
  inode[1] = 2;
  jnode[0] = 1;
  jnode[1] = 3;

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, eccent, parity );

  cout << "\n";
  cout << "  Parity = " << parity << "\n";
  cout << "  Eccentricity is " << eccent << "\n";

  if ( parity == 0 )
  {
    cout << "  No center node (degenerate case).\n";
  }
  else if ( parity == 1 )
  {
    cout << "  Center node: " << center[0] << "\n";
  }
  else
  {
    cout << "  Center nodes: " << center[0] << "  " << center[1] << "\n";
  }

  return;
# undef NNODE
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests TREE_ARC_CENTER.
//
//  Discussion:
//
//    The tree is
//
//     1-----2-----3
//    /|\   / \   /|\
//   4 5 6 8  10 7 9 11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 11

  int center[2];
  int eccent;
  int i;
  int inode[NNODE-1] = { 1, 1, 1, 2,  2, 3, 3,  3, 1, 2 };
  int jnode[NNODE-1] = { 4, 5, 6, 8, 10, 7, 9, 11, 2, 3 };
  int nnode = NNODE;
  int parity;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  TREE_ARC_CENTER computes the center of a tree.\n";

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_center ( nnode, inode, jnode, center, eccent, parity );

  cout << "\n";
  cout << "  Parity = " << parity << "\n";
  cout << "  Eccentricity is " << eccent << "\n";

  if ( parity == 0 )
  {
    cout << "  No center node (degenerate case).\n";
  }
  else if ( parity == 1 )
  {
    cout << "  Center node: " << center[0] << "\n";
  }
  else
  {
    cout << "  Center nodes: " << center[0] << "  " << center[1] << "\n";
  }

  return;
# undef NNODE
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests TREE_ARC_DIAM.
//
//  Discussion:
//
//    The tree is:
//
//    2---3---6---8---1---9
//       /       / \
//      7       5   4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int diam;
  int inode[8] = { 2, 3, 3, 6, 8, 8, 8, 1 };
  int jnode[8] = { 3, 7, 6, 8, 4, 5, 1, 9 };
  int label[9];
  int nnode = 9;
  int nnode1;
  int nnode2;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  TREE_ARC_DIAM computes the diameter of a tree.\n";

  graph_arc_print ( nnode-1, inode, jnode, "  The edge list of the tree:" );

  tree_arc_diam ( nnode, inode, jnode, diam, label, nnode1, nnode2 );

  cout << "\n";
  cout << "  This tree has a diameter of " << diam << "\n";
  cout << "  between nodes " << nnode1 << " and " << nnode2 << "\n";

  i4vec_print ( nnode, label, "  Nodes and labels:" );

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests TREE_ARC_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int icode[2];
  int inode[3];
  int jnode[3];
  int nnode = 4;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  TREE_ARC_RANDOM produces a random labeled\n";
  cout << "  tree and its Pruefer code.\n";
  cout << "\n";
 
  for ( i = 1; i <= 5; i++ )
  {
    tree_arc_random ( nnode, seed, icode, inode, jnode );

    graph_arc_print ( nnode-1, inode, jnode, "  The random tree:" );

    i4vec_print ( nnode-2, icode, "  The Pruefer code:" );
  }
  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests TREE_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int nnode;
  int num;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  TREE_ENUM enumerates the labeled trees on a given\n";
  cout << "  number of nodes.\n";
  cout << "\n";

  for ( nnode = 0; nnode <= 10; nnode++ )
  {
    num = tree_enum ( nnode );
    cout << "  " << setw(8) << nnode
         << "  " << setw(10) << num << "\n";
  } 
  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests TREE_PARENT_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
# define NNODE 4

  int icode[NNODE];
  int itree[NNODE];
  int more;
  int nnode = NNODE;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  TREE_PARENT_NEXT finds all labeled trees of a given \n";
  cout << "  order, and their Pruefer codes.\n";
  cout << "\n";
  cout << "  Pruefer code     Tree\n";
  cout << "\n";
 
  more = 0;
 
  for ( ; ; )
  {
    tree_parent_next ( nnode, icode, itree, more );
 
    cout << "  " << icode[0]
         << "  " << icode[1]
         << "            "
         << "  " << setw(2) << itree[0]
         << "  " << setw(2) << itree[1]
         << "  " << setw(2) << itree[2] << "\n";

    if ( ! more )
    {
      break;
    }
  }
  return;
# undef NNODE
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests TREE_RB_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int nnode;
  int num;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  TREE_RB_ENUM enumerates the rooted binary trees on a \n";
  cout << "  given number of nodes.\n";
  cout << "\n";

  for ( nnode = 0; nnode <= 11; nnode++ )
  {
    num = tree_rb_enum ( nnode );

    cout << "  " << setw(8) << nnode
         << "  " << setw(8) << num << "\n";
  }
  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests TREE_RB_LEX_NEXT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int a[11];
  int i;
  int j;
  int more;
  int n = 11;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  TREE_RB_LEX_NEXT produces all rooted binary trees with\n";
  cout << "  a given number of nodes, in lexicographic order, using\n";
  cout << "  the preorder traversal representation.\n";
  cout << "\n";
  cout << "  The number of nodes N = " << n << "\n";
  cout << "\n";

  more = 0;
  i = 0;

  for ( ; ; )
  {
    tree_rb_lex_next ( n, a, more );

    if ( ! more )
    {
      break;
    }

    i = i + 1;
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j < n; j++ )
    {
      cout << a[j];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests TREE_RB_LEX_NEXT, TREE_RB_TO_PARENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int a[11];
  int i;
  int j;
  int more;
  int n = 11;
  int *parent;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  TREE_RB_LEX_NEXT produces all rooted binary trees with\n";
  cout << "  a given number of nodes, in lexicographic order,\n";
  cout << "  using the preorder traversal representation.\n";
  cout << "  TREE_RB_TO_PARENT converts the preorder traversal form\n";
  cout << "  to the more comprehensible parent node representation.\n";
  cout << "\n";
  cout << "  The number of nodes N = " << n << "\n";
  cout << "\n";

  more = 0;
  i = 0;

  for ( ; ; )
  {
    tree_rb_lex_next ( n, a, more );

    if ( ! more )
    {
      break;
    }

    parent = tree_rb_to_parent ( n, a );

    i = i + 1;
    cout << "  " << setw(2) << i << "  ";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(3) << parent[j];
    }
    cout << "\n";

    delete [] parent;
  }
  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests TREE_RB_YULE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int a[11];
  int i;
  int j;
  int n;
  int n_max = 11;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  TREE_RB_YULE carries out one step of the Yule model\n";
  cout << "  on a rooted binary tree stored in preorder traversal form.\n";
  cout << "\n";
  cout << "  Each call adds two children to an arbitary leaf.\n";

  for ( i = 1; i <= 5; i++ )
  {
    cout << "\n";
    cout << "  Simulation " << i << "\n";
    cout << "\n";
    cout << "  Nodes  Preorder code\n";
    cout << "\n";

    n = 0;

    for ( ; ; )
    {
      tree_rb_yule ( n, seed, a );

      cout << "  " << setw(2) << n << "  ";
      for ( j = 0; j < n; j++ )
      {
        cout << a[j];
      }
      cout << "\n";

      if ( n_max < n + 2 )
      {
        break;
      }
    }
  }
  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests TREE_ROOTED_CODE.
//
//  Discussion:
//
//      1
//      |\
//      | \
//      |  \
//      2   3
//     /|\  |\
//    4 5 6 7 8
//     /|  \
//    9 10  11
//      |
//      12
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int *code;
  int nnode = 12;
  int parent[12] = { 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 };

  cout << "\n";
  cout << "TEST15\n";
  cout << "  TREE_ROOTED_CODE: code of a rooted tree.\n";

  i4vec_print ( nnode, parent, "  Parent vector for tree:" );

  code = tree_rooted_code ( nnode, parent );

  i4vec_print ( 2*nnode, code, "  The tree code:" );

  delete [] code;

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests TREE_ROOTED_DEPTH.
//
//  Discussion:
//
//      1
//      |\
//      | \
//      |  \
//      2   3
//     /|\  |\
//    4 5 6 7 8
//     /|  \
//    9 10  11
//      |
//      12
//
//    Depths
//
//    1  2  3  4  5  6  7  8  9 10 11 12
//    0  1  1  2  2  2  2  2  3  3  3  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int depth;
  int *depth_node;
  int nnode = 12;
  int parent[12] = { 0, 1, 1, 2, 2, 2, 3, 3, 5, 5, 6, 10 };

  cout << "\n";
  cout << "TEST16\n";
  cout << "  TREE_ROOTED_DEPTH: depth of a rooted tree.\n";

  i4vec_print ( nnode, parent, "  Parent vector for tree:" );

  depth_node = new int[nnode];

  tree_rooted_depth ( nnode, parent, depth, depth_node );

  i4vec_print ( nnode, depth_node, "  Individual node depths:" );

  cout << "\n";
  cout << "  Overall rooted tree depth: " << depth << "\n";

  delete [] depth_node;

  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests TREE_ROOTED_ENUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int nnode = 10;
  int *ntree;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  TREE_ROOTED_ENUM counts unlabeled rooted trees.\n";

  ntree = tree_rooted_enum ( nnode );

  i4vec_print ( nnode, ntree, 
    "  Number of trees with given number of nodes:" );

  delete [] ntree;

  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests TREE_ROOTED_RANDOM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int *itree;
  int j;
  int nnode = 5;
  int seed;

  seed = 123456789;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  TREE_ROOTED_RANDOM: random unlabeled rooted trees.\n";
  cout << "\n";
  cout << "  Selecting random trees, rooted at 1\n";
  cout << "  Number of nodes is NNODE = " << nnode << "\n";
  cout << "\n";
  cout << "  Each tree is described by the nodes that\n";
  cout << "  connect nodes 2 through NNODE.\n";
  cout << "\n";
  for ( i = 1; i <= 5; i++ )
  {
    itree = tree_rooted_random ( nnode, seed );

    cout << "  ";
    for ( j = 1; j < nnode; j++ )
    {
      cout << setw(4) << itree[j];
    }
    cout << "\n";

    delete [] itree;
  }
  return;
}
