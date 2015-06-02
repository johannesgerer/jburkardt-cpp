# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "treepack.hpp"

//****************************************************************************80

int *catalan ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN computes the Catalan numbers, from C(0) to C(N).
//
//  Discussion:
//
//    The Catalan number C(N) counts:
//
//    1) the number of binary trees on N vertices;
//    2) the number of ordered trees on N+1 vertices;
//    3) the number of full binary trees on 2N+1 vertices;
//    4) the number of well formed sequences of 2N parentheses;
//    5) the number of ways 2N ballots can be counted, in order,
//       with N positive and N negative, so that the running sum
//       is never negative;
//    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
//    7) the number of monotone functions from [1..N} to [1..N} which
//       satisfy f(i) <= i for all i;
//    8) the number of ways to triangulate a polygon with N+2 vertices.
//
//    The formula is:
//
//      C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
//           = 1 / (N+1) * COMB ( 2N, N )
//           = 1 / (2N+1) * COMB ( 2N+1, N+1).
//
//  First values:
//
//     C(0)     1
//     C(1)     1
//     C(2)     2
//     C(3)     5
//     C(4)    14
//     C(5)    42
//     C(6)   132
//     C(7)   429
//     C(8)  1430
//     C(9)  4862
//    C(10) 16796
//
//  Recursion:
//
//    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
//    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
//
//  Example:
//
//    N = 3
//
//    ()()()
//    ()(())
//    (()())
//    (())()
//    ((()))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer, 1986,
//    ISBN: 0387963472,
//    LC: QA164.S79.
//
//  Parameters:
//
//    Input, int N, the number of Catalan numbers desired.
//
//    Output, int CATALAN[N+1], the Catalan numbers from C(0) to C(N).
//
{
  int *c;
  int i;

  if ( n < 0 )
  {
    return NULL;
  }

  c = new int[n+1];

  c[0] = 1;
//
//  The extra parentheses ensure that the integer division is
//  done AFTER the integer multiplication.
//
  for ( i = 1; i <= n; i++ )
  {
    c[i] = ( c[i-1] * 2 * ( 2 * i - 1 ) ) / ( i + 1 );
  }

  return c;
}
//****************************************************************************80

void catalan_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_VALUES returns some values of the Catalan numbers for testing.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the order of the Catalan number.
//
//    Output, int &C, the value of the Catalan number.
//
{
# define N_MAX 11

  int c_vec[N_MAX] = { 1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 };
  int n_vec[N_MAX] = { 0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data];
    c = c_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cbt_traverse ( int depth )

//****************************************************************************80
//
//  Purpose:
//
//    CBT_TRAVERSE traverses a complete binary tree of given depth.
//
//  Discussion:
//
//    There will be 2^DEPTH terminal nodes of the complete binary tree.
//
//    This function traverses the tree, and prints out a binary code of 0's
//    and 1's each time it encounters a terminal node.  This results in a 
//    printout of the binary digits from 0 to 2^DEPTH - 1.
//
//    The function is intended as a framework to be used to traverse a binary
//    tree.  Thus, in practice, a user would insert some action when a terminal
//    node is encountered.
//
//    Another use would occur when a combinatorial search is being made, for
//    example in a knapsack problem.  Each binary string then represents which
//    objects are to be included in the knapsack.  In that case, the traversal
//    could be speeded up by noticing cases where a nonterminal node has been
//    reached, but the knapsack is already full, in which case the only solution
//    uses none of the succeeding items, or overfull, in which case no solutions
//    exist that include this initial path segment.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 December 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEPTH, the depth of the tree.
//
{
  int *b;
  int direction;
  int DOWNLEFT = 1;
  int i;
  int k;
  int p;
  int UP = 3;
  int UPDOWNRIGHT = 2;

  if ( depth < 1 )
  {
    return;
  }

  b = new int[depth+1];

  for ( i = 0; i <= depth; i++ )
  {
    b[i] = 0;
  }
  p = 0;
  direction = DOWNLEFT;
  k = 0;

  for ( ; ; )
  {
//
//  Try going in direction DOWNLEFT.
//
    if ( direction == DOWNLEFT )
    {
      p = p + 1;
      b[p-1] = 0;
      if ( p < depth )
      {
        cout << "           ";
        for ( i = 0; i < p; i++ )
        {
          cout << b[i];
        }
        cout << "\n";
      }
      else
      {
        cout << "  (  " << setw(4) << k << "  ";
        for ( i = 0; i < p; i++ )
        {
          cout << b[i];
        }
        cout << "\n";
        k = k + 1;
        direction = UPDOWNRIGHT;
      }
    }
//
//  Try going in direction UPDOWNRIGHT.
//
    if ( direction == UPDOWNRIGHT )
    {
      b[p-1] = + 1;
      if ( p < depth )
      {
        cout << "           ";
        for ( i = 0; i < p; i++ )
        {
          cout << b[i];
        }
        cout << "\n";
        direction = DOWNLEFT;
      }
      else
      {
        cout << "  )  " << setw(4) << k << "  ";
        for ( i = 0; i < p; i++ )
        {
          cout << b[i];
        }
        cout << "\n";
        k = k + 1;
        direction = UP;
      }
    }
//
//  Try going in direction UP.
//
    if ( direction == UP )
    {
      p = p - 1;
      if ( 1 <= p )
      {
        cout << "           ";
        for ( i = 0; i < p; i++ )
        {
          cout << b[i];
        }
        cout << "\n";
        if ( b[p-1] == 0 )
        {
          direction = UPDOWNRIGHT;
        }
      }
      else
      {
        break;
      }
    }
  }

  delete [] b;

  return;
}
//****************************************************************************80

int graph_adj_edge_count ( int adj[], int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ADJ_EDGE_COUNT counts the number of edges in a graph.
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
//  Parameters:
//
//    Input, int ADJ[NNODE*NNODE], the adjacency information.
//    ADJ(I,J) is 1 if there is an edge from node I to node J.
//
//    Input, int NNODE, the number of nodes.
//
//    Output, int GRAPH_ADJ_EDGE_COUNT, the number of edges in the graph.
//
{
  int i;
  int j;
  int nedge;

  nedge = 0;

  for ( i = 0; i < nnode; i++ )
  {
    for ( j = 0; j < nnode; j++ )
    {
      if ( i == j )
      {
        nedge = nedge + 2 * adj[i+j*nnode];
      }
      else
      {
        nedge = nedge + adj[i+j*nnode];
      }
    }
  }

  nedge = nedge / 2;

  return nedge;
}
//****************************************************************************80

int graph_adj_is_node_connected ( int adj[], int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ADJ_IS_NODE_CONNECTED determines if a graph is nodewise connected.
//
//  Discussion:
//
//    A graph is nodewise connected if, from every node, there is a path
//    to any other node.
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
//  Parameters:
//
//    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the 
//    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
//
//    Input, int NNODE, the number of nodes.
//
//    Output, int GRAPH_ADJ_IS_NODE_CONNECTED.
//    0, the graph is not nodewise connected.
//    1, the graph is nodewise connected.
//
{
  int *found;
  int i;
  int ihi;
  int ii;
  int ilo;
  int j;
  int jhi;
  int jlo;
  int *list;
  int result;
//
//  FOUND(I) is 1 if node I has been reached.
//  LIST(I) contains a list of the nodes as they are reached.
//
  list = new int[nnode];
  found = new int[nnode];

  for ( i = 0; i < nnode; i++ )
  {
    list[i] = 0;
    found[i] = 0;
  }
//
//  Start at node 1.
//
  found[1-1] = 1;
  list[1-1] = 1;
  ilo = 1;
  ihi = 1;
//
//  From the batch of nodes found last time, LIST(ILO:IHI),
//  look for unfound neighbors, and store their indices in LIST(JLO:JHI).
//
  for ( ; ; )
  {
    jlo = ihi + 1;
    jhi = ihi;

    for ( ii = ilo; ii <= ihi; ii++ )
    {
      i = list[ii-1];

      for ( j = 1; j <= nnode; j++ )
      {
        if ( adj[i-1+(j-1)*nnode] != 0 || adj[j-1+(i-1)*nnode] != 0 )
        {
          if ( found[j-1] == 0 )
          {
            jhi = jhi + 1;
            list[jhi-1] = j;
            found[j-1] = 1;
          }
        }
      }
    }
//
//  If no neighbors were found, exit.
//
    if ( jhi < jlo )
    {
      break;
    }
//
//  If neighbors were found, then go back and find THEIR neighbors.
//
    ilo = jlo;
    ihi = jhi; 
  }
//
//  No more neighbors were found.  Have we reached all nodes?
//
  if ( ihi == nnode )
  {
    result = 1;
  }
  else
  {
    result = 0;
  }

  delete [] found;
  delete [] list;

  return result;
}
//****************************************************************************80

int graph_adj_is_tree ( int adj[], int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ADJ_IS_TREE determines whether a graph is a tree.
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
//  Parameters:
//
//    Input, int ADJ[NNODE*NNODE], the adjacency matrix for the 
//    graph.  ADJ(I,J) is nonzero if there is an edge from node I to node J.
//
//    Input, int NNODE, the number of nodes.
//
//    Output, int GRAPH_ADJ_IS_TREE.
//    0, the graph is not a tree.
//    1, the graph is a tree.
//
{
  int nedge;
  int result;

  if ( nnode <= 1 )
  {
    result = 1;
    return result;
  }
//
//  Every node must be connected to every other node.
//
  result = graph_adj_is_node_connected ( adj, nnode );

  if ( result == 0 )
  {
    return result;
  }
//
//  There must be exactly NNODE-1 edges.
//
  nedge = graph_adj_edge_count ( adj, nnode );

  if ( nedge == nnode - 1 )
  {
    result = 1;
  }
  else
  {
    result = 0;
  }

  return result;
}
//****************************************************************************80

int *graph_arc_degree ( int nnode, int nedge, int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_DEGREE determines the degree of the nodes of a graph.
//
//  Discussion:
//
//    The degree of a node is the number of edges that include the node.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE], the pairs of nodes
//    that form the edges.
//
//    Output, int GRAPH_ARC_DEGREE[NNODE], the degree of each node,
//    that is, the number of edges that include the node.
//
{
  int *degree;
  int i;
  int n;

  degree = new int[nnode];

  for ( i = 0; i < nnode; i++ )
  {
    degree[i] = 0;
  }

  for ( i = 0; i < nedge; i++ )
  {
    n = inode[i];

    if ( 1 <= n && n <= nnode )
    {
      degree[n-1] = degree[n-1] + 1;
    }
    else
    {
      cerr << "\n";
      cerr << "GRAPH_ARC_DEGREE - Fatal error!\n";
      cerr << "  Out-of-range node value = " << n << "\n";
      exit ( 1 );
    }

    n = jnode[i];
    if ( 1 <= n && n <= nnode )
    {
      degree[n-1] = degree[n-1] + 1;
    }
    else
    {
      cerr << "\n";
      cerr << "GRAPH_ARC_DEGREE - Fatal error!\n";
      cerr << "  Out-of-range node value = " << n << "\n";
      exit ( 1 );
    }
  }
  return degree;
}
//****************************************************************************80

int graph_arc_is_tree ( int nedge, int inode[], int jnode[], int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_IS_TREE determines whether a graph is a tree.
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
//  Parameters:
//
//    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
//    JNODE(I) are the start and end nodes of the I-th edge of the graph G.
//
//    Input, int NEDGE, the number of edges in the graph G.
//
//    Input, int NNODE, the number of nodes.
//
//    Output, int GRAPH_ARC_IS_TREE.
//    0, the graph is not a tree.
//    1, the graph is a tree.
//
{
  int *adj;
  int result;

  adj = graph_arc_to_graph_adj ( nedge, inode, jnode );

  result = graph_adj_is_tree ( adj, nnode );

  delete [] adj;

  return result;
}
//****************************************************************************80

int graph_arc_node_count ( int nedge, int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_NODE_COUNT counts the number of nodes in a graph.
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
//  Parameters:
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
//    JNODE(I) are the start and end nodes of the I-th edge.
//
//    Output, int GRAPH_ARC_NODE_COUNT, the number of distinct nodes.
//
{
  int i;
  int *knode;
  int nnode;
//
//  Copy all the node labels into KNODE,
//  sort KNODE,
//  count the unique entries.  
//
//  That's NNODE.
//
  knode = new int[2*nedge];

  for ( i = 0; i < nedge; i++ )
  {
    knode[i] = inode[i];
    knode[i+nedge] = jnode[i];
  }

  i4vec_sort_heap_a ( 2*nedge, knode );

  nnode = i4vec_sorted_unique_count ( 2*nedge, knode );

  delete [] knode;

  return nnode;
}
//****************************************************************************80

int graph_arc_node_max ( int nedge, int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_NODE_MAX determines the maximum node label.
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
//  Parameters:
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE].  INODE(I) and 
//    JNODE(I) are the start and end nodes of the I-th edge.
//
//    Output, int GRAPH_ARC_NODE_MAX, the maximum node index.
//
{
  int i;
  int node_max;

  node_max = 0;
  for ( i = 0; i < nedge; i++ )
  {
    if ( node_max < inode[i] )
    {
      node_max = inode[i];
    }
    if ( node_max < jnode[i] )
    {
      node_max = jnode[i];
    }
  }
  return node_max;
}
//****************************************************************************80

void graph_arc_print ( int nedge, int inode[], int jnode[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_PRINT prints out a graph from an edge list.
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
//  Parameters:
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE], the beginning 
//    and end nodes of the edges.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i < nedge; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(6) << inode[i]
         << "  " << setw(6) << jnode[i] << "\n";
  }

  return;
}
//****************************************************************************80

int *graph_arc_to_graph_adj ( int nedge, int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRAPH_ARC_TO_GRAPH_ADJ converts an arc list graph to an adjacency graph.
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
//  Parameters:
//
//    Input, int NEDGE, the number of edges.
//
//    Input, int INODE[NEDGE], JNODE[NEDGE], the edge array for 
//    an undirected graph.  The I-th edge connects nodes INODE(I) and JNODE(I).
//
//    Output, int GRAPH_ARC_TO_GRAPH_ADJ[NNODE*NNODE], the adjacency information.
//
{
  int *adj;
  int i;
  int j;
  int k;
  int nnode;
//
//  Determine the number of nodes.
//
  nnode = graph_arc_node_count ( nedge, inode, jnode );

  adj = new int[nnode*nnode];

  for ( j = 0; j < nnode; j++ )
  {
    for ( i = 0; i < nnode; i++ )
    {
      adj[i+j*nnode] = 0;
    }
  }

  for ( k = 0; k < nedge; k++ )
  {
    i = inode[k] - 1;
    j = jnode[k] - 1;
    adj[i+j*nnode] = 1;
    adj[j+i*nnode] = 1;
  }

  return adj;
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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2012
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void i4mat_print ( int m, int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT.
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
//    10 September 2009
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
  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
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
//    20 August 2010
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
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << "  " << setw(6) << j - 1;
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ":";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << "  " << setw(6) << a[i-1+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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
  for ( i = ( n / 2 ) - 1; 0 <= i; i-- )
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

int *i4vec_indicator_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR_NEW sets an I4VEC to the indicator vector.
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
//    03 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR_NEW[N], the array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }
  return a;
}
//****************************************************************************80

int i4vec_max ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MAX returns the value of the maximum element in an I4VEC.
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
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, int A[N], the array to be checked.
//
//    Output, int I4VEC_MAX, the value of the maximum element.  This
//    is set to 0 if N <= 0.
//
{
  int i;
  int value;

  if ( n <= 0 )
  {
    return 0;
  }

  value = a[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < a[i] )
    {
      value = a[i];
    }
  }

  return value;
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
         << ": " << setw(8) << a[i]  << "\n";
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

int i4vec_sorted_unique_count ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORTED_UNIQUE_COUNT counts unique elements in a sorted I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    Because the array is sorted, this algorithm is O(N).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, int A[N], the sorted array to examine.
//
//    Output, int I4VEC_SORTED_UNIQUE_COUNT, the number of unique elements of A.
//
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n < 1 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i-1] != a[i] )
    {
      unique_num = unique_num + 1;
    }
  }

  return unique_num;
}
//****************************************************************************80

void pruefer_to_tree_arc ( int nnode, int iarray[], int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRUEFER_TO_TREE_ARC is given a Pruefer code, and computes the tree.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer Verlag, New York, 1986.
//
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int IARRAY[NNODE-2], the Pruefer code of the tree.
//
//    Output, int INODE[NNODE-1], JNODE[NNODE-1], the edge
//    array of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
//
{
  int i;
  int ii;
  int *iwork;
  int j;
//
//  Initialize IWORK(I) to count the number of neighbors of node I.
//  The Pruefer code uses each node one less time than its total
//  number of neighbors.
//
  iwork = new int[nnode];

  for ( i = 0; i < nnode; i++ )
  {
    iwork[i] = 1;
  }
  for ( i = 0; i < nnode - 2; i++ )
  {
    iwork[iarray[i]-1] = iwork[iarray[i]-1] + 1;
  }

  for ( i = 0; i < nnode - 1; i++ )
  {
    inode[i] = -1;
    jnode[i] = -1;
  }
//
//  Now process each entry in the Pruefer code.
//
  for ( i = 0; i < nnode - 2; i++ )
  {
    ii = -1;
    for ( j = 0; j < nnode; j++ )
    {
      if ( iwork[j] == 1 )
      {
        ii = j;
      }
    }
    inode[i] = ii + 1;
    jnode[i] = iarray[i];
    iwork[ii] = 0;
    iwork[iarray[i]-1] = iwork[iarray[i]-1] - 1;
  }

  inode[nnode-2] = iarray[nnode-3];
 
  if ( iarray[nnode-3] != 1 )
  {
    jnode[nnode-2] = 1;
  }
  else
  {
    jnode[nnode-2] = 2;
  }

  delete [] iwork;

  return;
}
//****************************************************************************80

void pruefer_to_tree_2 ( int nnode, int iarray[], int itree[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRUEFER_TO_TREE_2 produces the edge list of a tree from its Pruefer code.
//
//  Discussion:
//
//    One can thus exhibit all trees on N nodes, produce
//    one at random, find the M-th one on the list, etc, by
//    manipulating the Pruefer codes.
//
//    For every labeled tree on N nodes, there is a unique N-2 tuple
//    of integers A1 through AN-2, with each A between 1 and N.  There
//    are N^(N-2) such sequences, and each one is associated with exactly
//    one tree.
//
//    This routine apparently assumes that the Pruefer code is
//    generated by taking the LOWEST labeled terminal node each time.
//    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis. Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int NNODE, number of nodes in desired tree.
//
//    Input, int IARRAY[NNODE].  IARRAY(I), I = 1, NNODE-2 
//    is the Pruefer code for the tree.
//
//    Output, int ITREE[NNODE]; the I-th edge of the tree
//    joins nodes I and ITREE(I).
//
{
  int i;
  int ir;
  int j;
  int k;
  int kp;
  int l;

  for ( i = 0; i < nnode; i++ )
  {
    itree[i] = 0;
  }
 
  for ( i = nnode - 2; 1 <= i; i-- )
  {
    l = iarray[i-1];
 
    if ( itree[l-1] == 0 )
    {
      iarray[i-1] = - l;
      itree[l-1] = - 1;
    }
  }
  iarray[nnode-2] = nnode;
//
//  Find next index K so that ITREE(K) is 0.
//
  k = 1;

 
  while ( itree[k-1] != 0 )
  {
    k = k + 1;
  }
 
  j = 0;
  kp = k;
 
  for ( ; ; )
   {
    j = j + 1;
    ir = abs ( iarray[j-1] );
    itree[kp-1] = ir;
 
    if ( j == nnode - 1 )
    {
      break;
    }
 
    if ( 0 < iarray[j-1] )
    {
      while ( itree[k-1] != 0 )
      {
        k = k + 1;
      }
      kp = k;
      continue;
    }
 
    if ( k < ir )
    {
      itree[ir-1] = 0;
      while ( itree[k-1] != 0 )
      {
        k = k + 1;
      }
      kp = k;
      continue;
    } 
    kp = ir;
  }
//
//  Restore the signs of IARRAY.
//
  for ( i = 0; i < nnode - 2; i++ )
  {
    iarray[i] = abs ( iarray[i] );
  }
 
  return;
}
//****************************************************************************80

int *pruefer_to_tree_2_new ( int nnode, int iarray[] )

//****************************************************************************80
//
//  Purpose:
//
//    PRUEFER_TO_TREE_2_NEW produces the edge list of a tree from its Pruefer code.
//
//  Discussion:
//
//    One can thus exhibit all trees on N nodes, produce
//    one at random, find the M-th one on the list, etc, by
//    manipulating the Pruefer codes.
//
//    For every labeled tree on N nodes, there is a unique N-2 tuple
//    of integers A1 through AN-2, with each A between 1 and N.  There
//    are N^(N-2) such sequences, and each one is associated with exactly
//    one tree.
//
//    This routine apparently assumes that the Pruefer code is
//    generated by taking the LOWEST labeled terminal node each time.
//    This is not consistent with PRUEFER_TO_TREE and TREE_TO_PRUEFER.
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Albert Nijenhuis. Herbert Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int NNODE, number of nodes in desired tree.
//
//    Input, int IARRAY[NNODE].  IARRAY(I), I = 1, NNODE-2 
//    is the Pruefer code for the tree.
//
//    Output, int PRUEFER_TO_TREE_2_NEW[NNODE]; the I-th edge of the tree
//    joins nodes I and ITREE(I).
//
{
  int *itree;

  itree = new int[nnode];

  pruefer_to_tree_2 ( nnode, iarray, itree );
 
  return itree;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

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
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
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
//****************************************************************************80

void tree_arc_center ( int nnode, int inode[], int jnode[], int center[], 
  int &eccent, int &parity )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ARC_CENTER computes the center, eccentricity, and parity of a tree.
//
//  Discussion:
//
//    A tree is an undirected graph of N nodes, which uses N-1 edges,
//    and is connected.  
//
//    A graph with N-1 edges is not guaranteed to be a tree, and so this
//    routine must first check that condition before proceeding.
//
//    The edge distance between two nodes I and J is the minimum number of
//    edges that must be traversed in a path from I and J.
//
//    The eccentricity of a node I is the maximum edge distance between
//    node I and the other nodes J in the graph.
//
//    The radius of a graph is the minimum eccentricity over all nodes
//    in the graph.
//
//    The diameter of a graph is the maximum eccentricity over all nodes
//    in the graph.
//
//    The center of a graph is the set of nodes whose eccentricity is 
//    equal to the radius, that is, the set of nodes of minimum eccentricity.
//
//    For a tree, the center is either a single node, or a pair of
//    neighbor nodes.
//
//    The parity of the tree is 1 if the center is a single node, or 2 if
//    the center is 2 nodes.
//
//    The center of a tree can be found by removing all "leaves", that is,
//    nodes of degree 1.  This step is repeated until only 1 or 2 nodes
//    are left.
//
//    Thanks to Alexander Sax for pointing out that a previous version of the
//    code was failing when the tree had an odd parity, that is, a single
//    center node, 15 April 2013.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges of
//    the tree.  Edge I connects nodes INODE(I) and JNODE(I).
//
//    Output, int CENTER[2].  CENTER(1) is the index of the
//    first node in the center.  CENTER(2) is 0 if there is only one node
//    in the center, or else the index of the second node.
//
//    Output, int &ECCENT, the eccentricity of the nodes in 
//    the center, and the radius of the the tree.
//
//    Output, int &PARITY, the parity of the tree, which is
//    normally 1 or 2.
//
{
  int *degree;
  int i;
  int iedge;
  int ileaf;
  int j;
  int *list;
  int nedge;
  int nleaf;
  int nnode2;
  int result;

  eccent = 0;
  center[0] = 0;
  center[1] = 0;
  parity = 0;

  if ( nnode <= 0 )
  {
    cerr << "\n";
    cerr << "TREE_ARC_CENTER - Fatal error!\n";
    cerr << "  NNODE <= 0.\n";
    exit ( 1 );
  }
  else if ( nnode == 1 )
  {
    eccent = 0;
    center[0] = 1;
    center[1] = 0;
    parity = 1;
    return;
  }
  else if ( nnode == 2 )
  {
    eccent = 1;
    center[0] = 1;
    center[1] = 2;
    parity = 2;
    return;
  }
//
//  Is this graph really a tree?
//
  nedge = nnode - 1;
  result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

  if ( result == 0 )
  {
    cerr << "\n";
    cerr << "TREE_ARC_CENTER - Fatal error!\n";
    cerr << "  This graph is NOT a tree.\n";
    exit ( 1 );
  }
//
//  Compute the degrees.
//
  degree = graph_arc_degree ( nnode, nedge, inode, jnode );
//
//  Defoliate the tree.
//
  nnode2 = nnode;
  list = new int[nnode];

  for ( ; ; )
  {
    eccent = eccent + 1;
//
//  Find and mark the leaves.
//
    nleaf = 0;

    for ( i = 1; i <= nnode; i++ )
    {
      if ( degree[i-1] == 1 )
      {
        nleaf = nleaf + 1;
        list[nleaf-1] = i;
      }
    }
//
//  Delete the leaves.
//
    for ( ileaf = 1; ileaf <= nleaf; ileaf++ )
    {
      i = list[ileaf-1];

      iedge = 0;
      j = 0;

      for ( ; ; )
      {
        iedge = iedge + 1;

        if ( nedge < iedge )
        {
          cerr << "\n";
          cerr << "TREE_ARC_CENTER - Fatal error!\n";
          cerr << "  Data or algorithm failure.\n";
          exit ( 1 );
        }

        if ( inode[iedge-1] == i )
        {
          j = jnode[iedge-1];
          inode[iedge-1] = - inode[iedge-1];
          jnode[iedge-1] = - jnode[iedge-1];
        }
        else if ( jnode[iedge-1] == i )
        {
          j = inode[iedge-1];
          inode[iedge-1] = - inode[iedge-1];
          jnode[iedge-1] = - jnode[iedge-1];
        }

        if ( j != 0 )
        {
          break;
        }
      }

      degree[i-1] = -1;
      nnode2 = nnode2 - 1;
      degree[j-1] = degree[j-1] - 1;
//
//  If the other node has degree 0, we must have just finished
//  stripping all leaves from the tree, leaving a single node.
//  Don't kill it here.  It is our odd center.
//
//     if ( degree(j) == 0 )
//     {
//       nnode2 = nnode2 - 1
//     }
//
    }
//
//  Find the remaining nodes.
//
    nnode2 = 0;

    for ( i = 1; i <= nnode; i++ )
    {
      if ( 0 <= degree[i-1] )
      {
        nnode2 = nnode2 + 1;
        list[nnode2-1] = i;
      }
    }
//
//  If at least 3, more pruning is required.
//
    if ( nnode2 < 3 )
    {
      break;
    }
  }
//
//  If only one or two nodes left, we are done.
//
  parity = nnode2;

  for ( i = 0; i < nnode2; i++ )
  {
    center[i] = list[i];
  }
  for ( i = 0; i < nedge; i++ )
  {
    inode[i] = abs ( inode[i] );
    jnode[i] = abs ( jnode[i] );
  }

  delete [] list;

  return;
}
//****************************************************************************80

void tree_arc_diam ( int nnode, int inode[], int jnode[], int &diam, 
  int label[], int &n1, int &n2 )

//****************************************************************************80
/*
  Purpose:

    TREE_ARC_DIAM computes the "diameter" of a tree.

  Discussion:

    A tree is an undirected graph of N nodes, which uses N-1 edges,
    and is connected.  

    A graph with N-1 edges is not guaranteed to be a tree, and so this
    routine must first check that condition before proceeding.

    The diameter of a graph is the length of the longest possible
    path that never repeats an edge.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    05 August 2013

  Author:

    John Burkardt

  Parameters:

    Input, int NNODE, the number of nodes.

    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edges 
    of the tree.  Edge I connects nodes INODE(I) and JNODE(I).

    Output, int &DIAM, the length of the longest path 
    in the tree.

    Output, int LABEL[NNODE], marks the path between 
    nodes N1 and N2.  Node I is in this path if LABEL(I) is 1.

    Output, int &N1, &N2, the indices of two nodes in the 
    tree which are separated by DIAM edges.
*/
{
  int *degree;
  int i;
  int invals;
  int j;
  int k;
  int kstep;
  int nabe;
  int nedge;
  int result;

  if ( nnode <= 0 )
  {
    diam = 0;
    cerr << "\n";
    cerr << "TREE_ARC_DIAM - Fatal error!\n";
    cerr << "  NNODE <= 0.\n";
    exit ( 1 );
  }

  if ( nnode == 1 )
  {
    diam = 0;
    n1 = 1;
    n2 = 1;
    return;
  }
//
//  Is this graph really a tree?
//
  nedge = nnode - 1;
  result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

  if ( result == 0 )
  {
    cerr << "\n";
    cerr << "TREE_ARC_DIAM - Fatal error!\n";
    cerr << "  This graph is NOT a tree.\n";
    exit ( 1 );
  }

  label = i4vec_indicator_new ( nnode );
//
//  On step K:
//
//    Identify the terminal and interior nodes.
//
//    If there are no interior nodes left, 
//
//      then there are just two nodes left at all.  The diameter is 2*K-1, 
//      and a maximal path extends between the nodes whose labels are 
//      contained in the two remaining terminal nodes.
//
//    Else
//
//      The label of each terminal node is passed to its interior neighbor.
//      If more than one label arrives, take any one.
//
//      The terminal nodes are removed.
//
  kstep = 0;
  degree = new int[nnode];

  for ( ; ; )
  {
    kstep = kstep + 1;
//
//  Compute the degree of each node.
//
    for ( j = 1; j <= nnode; j++ )
    {
      degree[j-1] = 0;   
    }
    for ( j = 1; j <= nedge; j++ )
    {
      k = inode[j-1];
      if ( 0 < k )
      {
        degree[k-1] = degree[k-1] + 1;
      }

      k = jnode[j-1];
      if ( 0 < k )
      {
        degree[k-1] = degree[k-1] + 1;
      }
    }
//
//  Count the number of interior nodes.
//
    invals = 0;
    for ( i = 1; i <= nnode; i++ )
    {
      if ( 1 < degree[i-1] )
      {
        invals = invals + 1;
      }
    }
//
//  If there are 1 or 0 interior nodes, it's time to stop.
//
    if ( invals == 1 )
    {
      diam = 2 * kstep;
      break;
    }
    else if ( invals == 0 )
    {
      diam = 2 * kstep - 1;
      break;
    }
//
//  If there are at least two interior nodes, then chop off the 
//  terminal nodes and pass their labels inward.
//
    for ( k = 1; k <= nnode; k++ )
    {
      if ( degree[k-1] == 1 )
      {
        for ( j = 1; j <= nedge; j++ )
        {
          if ( inode[j-1] == k )
          {
            nabe = jnode[j-1];
            label[nabe-1] = label[k-1];
            inode[j-1] = - inode[j-1];
            jnode[j-1] = - jnode[j-1];
          }
          else if ( jnode[j-1] == k )
          {
            nabe = inode[j-1];
            label[nabe-1] = label[k-1];
            inode[j-1] = - inode[j-1];
            jnode[j-1] = - jnode[j-1];
          }
        }
      }
    }
  }
//
//  Now get the labels from two of the remaining terminal nodes.
//  The nodes represented by these labels will be a diameter apart.
//
  n1 = 0;
  n2 = 0;

  for ( i = 1; i <= nnode; i++ );
  {
    if ( degree[i-1] == 1 )
    {
      if ( n1 == 0 )
      {
        n1 = label[i-1];
      }
      else if ( n2 == 0 )
      {
        n2 = label[i-1];
      }
    }
  }
//
//  Set the labels of the interior node (if any) and nodes marked
//  N1 and N2 to 1, and all others to 0.  This will label the nodes on the path.
//
  if ( invals == 1 )
  {
    for ( i = 1; i <= nnode; i++ )
    {
      if ( 1 < degree[i-1] )
      {
        label[i-1] = 1;
      }
    }
  }

  for ( i = 1; i <= nnode; i++ )
  {
    if ( label[i-1] == n1 || label[i-1] == n2 )
    {
      label[i-1] = 1;
    }
    else
    {
      label[i-1] = 0;
    }
  }
//
//  Clean up the arrays.
//
  for ( j = 1; j <= nedge; j++ )
  {
    inode[j-1] = abs ( inode[j-1] );
    k = inode[j-1];
    jnode[j-1] = abs ( jnode[j-1] );
    k = jnode[j-1];
  }

  delete [] degree;

  return;
}
//****************************************************************************80

void tree_arc_random ( int nnode, int &seed, int code[], int inode[], 
  int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ARC_RANDOM selects a random labeled tree and its Pruefer code.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input/output, int &SEED, a seed for the random number 
//    generator.
//
//    Output, int CODE[NNODE-2], the Pruefer code for the 
//    labeled tree.
//
//    Output, int INODE[NNODE-1], JNODE[NNODE-1], the edge 
//    array for the tree.
//
{
  int i;

  if ( nnode <= 2 )
  {
    return;
  }

  vec_random ( nnode-2, nnode, seed, code );
 
  for ( i = 0; i < nnode - 2; i++ )
  {
    code[i] = code[i] + 1;
  }
  pruefer_to_tree_arc ( nnode, code, inode, jnode );
 
  return;
}
//****************************************************************************80

int *tree_arc_to_pruefer ( int nnode, int inode[], int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ARC_TO_PRUEFER is given a labeled tree, and computes its Pruefer code.
//
//  Discussion:
//
//    A tree is an undirected graph of N nodes, which uses N-1 edges,
//    and is connected.  
//
//    A graph with N-1 edges is not guaranteed to be a tree, and so this
//    routine must first check that condition before proceeding.
//
//    The Pruefer code is a correspondence between all labeled trees of
//    N nodes, and all list of N-2 integers between 1 and N (with repetition
//    allowed).  The number of labeled trees on N nodes is therefore N^(N-2).
//
//    The Pruefer code is constructed from the tree as follows:
//
//    A terminal node on the tree is defined as a node with only one neighbor.
//
//    Consider the set of all terminal nodes on the tree.  Take the one
//    with the highest label, I.  Record the label of its neighbor, J.
//    Delete node I and the edge between node I and J.
//
//    J is the first entry in the Pruefer code for the tree.   Repeat
//    the operation a total of N-2 times to get the complete code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Reference:
//
//    Dennis Stanton, Dennis White,
//    Constructive Combinatorics,
//    Springer Verlage, New York, 1986.
//
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int INODE[NNODE-1], JNODE[NNODE-1], the edge array 
//    of the tree.  The I-th edge joins nodes INODE(I) and JNODE(I).
//
//    Output, int TREE_ARC_TO_PRUEFER[NNODE-2], the Pruefer code of the tree.
//
{
  int *code;
  int *degree;
  int i;
  int i2;
  int iterm;
  int j;
  int jsave;
  int nedge;
  int result;
//
//  Is this graph really a tree?
//
  nedge = nnode - 1;
  result = graph_arc_is_tree ( nedge, inode, jnode, nnode );

  if ( result == 0 )
  {
    cerr << "\n";
    cerr << "TREE_ARC_TO_PRUEFER - Fatal error!\n";
    cerr << "  This graph is NOT a tree.\n";
    exit ( 1 );
  }

  code = new int[nnode-2];
//
//  Compute the degree of each node.
//
  nedge = nnode - 1;
  degree =  graph_arc_degree ( nnode, nedge, inode, jnode );
//
//  Compute the next term of the Pruefer code.
//
  for ( i = 1; i <= nnode - 2; i++ )
  {
//
//  Find the terminal node with the highest label.
//
    iterm = 0;
 
    for ( j = 1; j <= nnode; j++ )
    {
      if ( degree[j-1] == 1 )
      {
        iterm = j;
      }
    }
//
//  Find the edge that includes this node, and note the
//  index of the other node.
//
    for ( j = 1; j < nnode - 1; j++ )
    {
      jsave = j;
 
      if ( inode[j-1] == iterm )
      {
        i2 = 2;
        break;
      }
      else if ( jnode[j-1] == iterm )
      {
        i2 = 1;
        break;
      }
    }
//
//  Delete the edge from the tree.
//
    degree[inode[jsave-1]-1] = degree[inode[jsave-1]-1] - 1;
    degree[jnode[jsave-1]-1] = degree[jnode[jsave-1]-1] - 1;
//
//  Add the neighbor of the node to the Pruefer code.
//
    if ( i2 == 1 )
    {
      code[i-1] = inode[jsave-1];
    }
    else
    {
      code[i-1] = jnode[jsave-1];
    }
//
//  Negate the nodes in the edge list to mark them as used.
//
    inode[jsave-1] = - inode[jsave-1];
    jnode[jsave-1] = - jnode[jsave-1];
  }
//
//  Before returning, restore the original form of the edge list.
//
  for ( i = 1; i <= nnode - 1; i++ )
  {
    inode[i-1] = abs ( inode[i-1] );
    jnode[i-1] = abs ( jnode[i-1] );
  }

  delete [] degree;

  return code;
}
//****************************************************************************80

int tree_enum ( int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ENUM enumerates the labeled trees on NNODE nodes.
//
//  Discussion:
//
//    The formula is due to Cauchy.
//
//  Example:
//
//    NNODE      NTREE
//
//    0              1
//    1              1
//    2              1
//    3              3
//    4             16
//    5            125
//    6           1296
//    7          16807
//    8         262144
//    9        4782969
//   10      100000000
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes in each tree.
//    NNODE must normally be at least 3, but for this routine,
//    any value of NNODE is allowed.  Values of NNODE greater than 10
//    will probably overflow.
//
//    Output, int TREE_ENUM, the number of distinct labeled trees.
//
{
  int ntree;

  if ( nnode < 0 )
  {
    ntree = 0;
  }
  else if ( nnode == 0 )
  {
    ntree = 1;
  }
  else if ( nnode == 1 )
  {
    ntree = 1;
  }
  else if ( nnode == 2 )
  {
    ntree = 1;
  }
  else
  {
    ntree = i4_power ( nnode, nnode - 2 );
  }

  return ntree;
}
//****************************************************************************80

void tree_parent_next ( int nnode, int code[], int itree[], int &more )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_PARENT_NEXT generates, one at a time, all labeled trees.
//
//  Discussion:
//
//    The routine also returns the corresponding Pruefer codes.
//
//  Formula:
//
//    There are N^(N-2) labeled trees on N nodes (Cayley's formula).
//
//    The number of trees in which node I has degree D(I) is the
//    multinomial coefficient: ( N-2; D(1)-1, D(2)-1, ..., D(N)-1 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Parameters:
//
//    Input, int NNODE, the number of nodes to be used in 
//    the trees.
//
//    Input/output, int CODE[NNODE].  The first NNODE-2 entries 
//    of CODE contain the Pruefer code for the given labeled tree.
//
//    Output, int ITREE[NNODE].  The first NNODE-1 entries 
//    of ITREE describe the edges that go between the nodes.  Each pair
//    (I, ITREE(I)) represents an edge.  Thus if ITREE(5) = 3,
//    there is an edge from node 3 to node 5.
//
//    Input/output, int &MORE.  On the first call only, the
//    user is required to set MORE = .FALSE.  Then call the routine, and
//    it will return information about the first tree
//    as well as setting MORE to the value .TRUE.
//    Keep calling to get another tree until MORE is .FALSE.
//    on return, at which point there are no more trees.
//
{
  int i;

  if ( more )
  {
    for ( i = 0; i < nnode - 2; i++ )
    {
      code[i] = code[i] - 1;
     }
  }

  vec_next ( nnode-2, nnode, code, more );
 
  for ( i = 0; i < nnode - 2; i++ )
  {
    code[i] = code[i] + 1;
  }

  pruefer_to_tree_2 ( nnode, code, itree );
 
  return;
}
//****************************************************************************80

void tree_parent_to_arc ( int nnode, int parent[], int &nedge, int inode[], 
  int jnode[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_PARENT_TO_ARC converts a tree from parent to arc representation.
//
//  Discussion:
//
//    Parent representation lists the parent node of each node.  For a
//    tree of N nodes, one node has a parent of 0, representing a null link.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes in the tree.
//
//    Input, int PARENT[NNODE], the parent node representation 
//    of the tree.
//
//    Output, int &NEDGE, the number of edges, normally NNODE-1.
//
//    Output, int INODE[NEDGE], JNODE[NEDGE], pairs of nodes
//    that define the links.
//
{
  int i;

  nedge = 0;

  for ( i = 1; i <= nnode; i++ )
  {
    if ( parent[i-1] != 0 )
    {
      nedge = nedge + 1;
      inode[nedge-1] = i;
      jnode[nedge-1] = parent[i-1];
    }
  }

  return;
}
//****************************************************************************80

int tree_rb_enum ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_RB_ENUM returns the number of rooted binary trees with N nodes.
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
//  Parameters:
//
//    Input, int N, the number of nodes in the rooted 
//    binary tree.  N should be odd.
//
//    Output, int TREE_RB_ENUM, the number of rooted binary trees 
//    with N nodes.
//
{
  int *c;
  int m;
  int num;

  if ( n < 0 )
  {
    num = 0;
  }
  else if ( n == 0 )
  {
    num = 1;
  }
  else if ( ( n % 2 ) == 0 )
  {
    num = 0;
  }
  else
  {
    m = ( n - 1 ) / 2;
    c = catalan ( m );
    num = c[m];
    delete [] c;
  }

  return num;
}
//****************************************************************************80

void tree_rb_lex_next ( int n, int a[], int &more )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_RB_LEX_NEXT generates rooted binary trees in lexicographic order.
//
//  Discussion:
//
//    The information definining the tree of N nodes is stored in a vector 
//    of 0's and 1's, in preorder traversal form.  Essentially, the
//    shape of the tree is traced out with a pencil that starts at the root,
//    and stops at the very last null leaf.  The first time that a (non-null) 
//    node is encountered, a 1 is added to the vector, and the left 
//    descendant of the node is visited next.  When the path returns from
//    the first descendant, the second descendant is visited.  When then path
//    returns again, the path passes back up from the node to its parent.
//    A null leaf is encountered only once, and causes a zero to be added to 
//    the vector, and the path goes back up to the parent node.  
//
//    The lexicographic order is used to order the vectors of 1's and 0's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Reference:
//
//    Frank Ruskey,
//    Combinatorial Generation,
//    To appear.
//
//  Parameters:
//
//    Input, int N, the number of nodes in the rooted binary
//    tree.  N should be odd.
//
//    Input/output, int A[N], the preorder traversal form for
//    the previous/next rooted binary tree.
//
//    Output, logical &MORE, is TRUE if the next rooted binary tree was
//    returned on this call, or FALSE if there are no more rooted binary
//    trees, and the output of the previous call was the last one.
//
{
  int i;
  int k;
  int p;
  int q;

  if ( ! more )
  {
    for ( i = 1; i <= n - 2; i = i + 2 )
    {
      a[i-1] = 1;
    }
    for ( i = 2; i <= n - 1; i = i + 2 )
    {
      a[i-1] = 0;
    }
    a[n-1] = 0;
    more = 1;
    return;
  }
//
//  Find the last 1 in A.
//
  k = n;
  while ( a[k-1] == 0 )
  {
    k = k - 1;
  }
  q = n - k - 1;
//
//  Find the last 0 preceding the last 1 in A.
//  If there is none, then we are done, because 11...1100..00 
//  is the final element.
//
  for ( ; ; )
  {
    if ( k == 1 )
    {
      more = 0;
      return;
    }

    if ( a[k-1] == 0 )
    {
      break;
    }
    k = k - 1;
  }

  p = n - k - q - 1;

  a[k-1] = 1;
  for ( i = k + 1; i <= n - 2 * p + 1; i++ )
  {
    a[i-1] = 0;
  }
  for ( i = n - 2 * p + 2; i <= n - 2; i = i + 2 )
  {
    a[i-1] = 1;
  }
  for ( i = n - 2 * p + 3; i <= n - 1; i = i + 2 )
  {
    a[i-1] = 0;
  }
  a[n-1] = 0;

  return;
}
//****************************************************************************80

int *tree_rb_to_parent ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_RB_TO_PARENT converts rooted binary tree to parent node representation.
//
//  Discussion:
//
//    Parent node representation of a tree assigns to each node a "parent" node,
//    which represents the first link of the path between the node and the 
//    root node.  The root node itself is assigned a parent of 0.
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
//  Parameters:
//
//    Input, int N, the number of nodes in the tree.
//
//    Input, int A[N], the preorder traversal form for the
//    rooted binary tree.
//
//    Output, int TREE_RB_TO_PARENT[N], the parent node representation 
//    of the tree.
//
{
  int dad;
  int k;
  int node;
  int node_num;
  int *parent;
  int *use;

  parent = new int[n];
  use = new int[n];

  node = 0;
  node_num = 0;

  for ( k = 1; k <= n; k++ )
  {
    dad = node;
    node_num = node_num + 1;
    node = node_num;
    parent[node-1] = dad;

    if ( a[k-1] == 1 )
    {
      use[node-1] = 0;
    }
    else
    {
      use[node-1] = 2;

      while ( use[node-1] == 2 )
      {
        node = dad;
        if ( node == 0 )
        {
          break;
        }
        use[node-1] = use[node-1] + 1;
        dad = parent[node-1];
      }
    }
  }
  delete [] use;

  return parent;
}
//****************************************************************************80

void tree_rb_yule ( int &n, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_RB_YULE adds two nodes to a rooted binary tree using the Yule model.
//
//  Discussion:
//
//    The Yule model is a simulation of how an evolutionary family tree
//    develops.  We start with a root node.  The internal nodes of the tree 
//    are inactive and never change.  Each pendant or leaf node of the
//    tree represents a biological family that can spontaneously "fission",
//    developing two new distinct sub families.  In graphical terms, the node
//    becomes internal, with two new leaf nodes depending from it.
//
//    The tree is stored in inorder traversal form.
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
//  Parameters:
//
//    Input/output, int &N, the number of nodes in the input
//    tree.  On output, this number has been increased, usually by 2.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Input/output, int A[*], the preorder traversal form 
//    for the rooted binary tree.  The number of entries in A is N.
//
{
  int i;
  int ileaf;
  int j;
  int jleaf;
  int nleaf;

  if ( n <= 0 )
  {
    n = 1;
    a[0] = 0;
    return;
  }
//
//  Count the expected number of leaves, which are the 0 values.
//
  nleaf = ( n + 1 ) / 2;
//
//  Choose a random number between 1 and NLEAF.
//
  ileaf = i4_uniform_ab ( 1, nleaf, seed );
//
//  Locate leaf number ILEAF.
//
  j = 0;
  jleaf = 0;
  for ( i = 1; i <= n; i++ )
  {
    if ( a[i-1] == 0 )
    {
      jleaf = jleaf + 1;
    }
    if ( jleaf == ileaf )
    {
      j = i;
      break;
    }
  }
//
//  Replace '0' by '100'
//
  for ( i = n; j <= i; i-- )
  {
    a[i+1] = a[i-1];
  }
  a[j-1] = 1;
  a[j] = 0;

  n = n + 2;

  return;
}
//****************************************************************************80

int *tree_rooted_code ( int nnode, int parent[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ROOTED_CODE returns the code of a rooted tree.
//
//  Discussion:
//
//    This code for a rooted tree depends on the node ordering, so it's actually
//    the code for a labeled rooted tree.  To eliminate the effects of node
//    labeling, one could choose as the code for a tree the maximum of all
//    the codes associated with the different possible labelings of the tree.
//    There are more effective ways of arriving at this code than simply
//    generating all possible codes and comparing them.  
//
//    For a tree with NNODES, the code is a list of 2*NNODE 0's and 1's,
//    describing a traversal of the tree starting at an imaginary node 0,
//    moving "down" to the root (a code entry of 1), and then moving
//    "down" (1) or "up" (0) as the tree is traversed in a depth first
//    manner.  The final move must be from the root up to the imaginary
//    node 0.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int PARENT[NNODE], is the parent node of each node.
//    The node with parent 0 is the root.
//
//    Output, int TREE_ROOTED_CODE[2*NNODE], the code for the tree.
//
{
  int *code;
  int father;
  int i;
  int k;
  int son;

  code = new int[2*nnode];
//
//  Find the root.
//
  father = 0;
  for ( i = 1; i <= nnode; i++ )
  {
    if ( parent[i-1] == 0 )
    {
      k = 1;
      code[0] = 1;
      father = i;
      break;
    }
  }

  if ( father == 0 )
  {
    cerr << "\n";
    cerr << "TREE_ROOTED_CODE - Fatal error!\n";
    cerr << "  Could not find the root.\n";
    exit ( 1 );
  }

  while ( father != 0 ) 
  {
    k = k + 1;
    code[k-1] = 0;
    for ( son = 1; son <= nnode; son++ )
    {
      if ( parent[son-1] == father )
      {
        code[k-1] = 1;
        father = son;
        break;
      }
    }

    if ( code[k-1] == 0 )
    {
      parent[father-1] = - parent[father-1];
      father = - parent[father-1];
    }
  }

  for ( i = 0; i < nnode; i++ )
  {
    parent[i] = - parent[i];
  }

  return code;
}
//****************************************************************************80

int tree_rooted_code_compare ( int nnode, int npart, int code1[], int code2[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ROOTED_CODE_COMPARE compares a portion of the code for two rooted trees.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int NPART, the number of nodes for which the code
//    has been determined.  This determines the portion of the codes to be
//    compared.  We expect 0 <= NPART <= NNODE.
//
//    Input, int CODE1[2*NNODE], CODE2[2*NNODE], the two 
//    rooted tree codes to be compared.
//
//    Output, int TREE_ROOTED_CODE_COMPARE, the result of the comparison.
//    -1, CODE1 < CODE2,
//     0, CODE1 = CODE2,
//    +1, CODE1 > CODE2.
//
{
  int i;
  int ihi;
  int result;

  result = 0;

  if ( npart <= 0 )
  {
    return result;
  }

  ihi = 2 * nnode;
  if ( npart < nnode )
  {
    ihi = 2 * npart;
  }

  for ( i = 0; i < ihi; i++ )
  {
    if ( code1[i] < code2[i] )
    {
      result = -1;
      return result;
    }
    else if ( code2[i] < code1[i] )
    {
      result = +1;
      return result;
    }
  }

  return result;
}
//****************************************************************************80

void tree_rooted_depth ( int nnode, int parent[], int &depth, int depth_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ROOTED_DEPTH returns the depth of a rooted tree.
//
//  Discussion:
//
//    The depth of any node of a rooted tree is the number of edges in 
//    the shortest path from the root to the node.
//
//    The depth of the rooted tree is the maximum of the depths
//    of all the nodes.
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
//  Parameters:
//
//    Input, int NNODE, the number of nodes.
//
//    Input, int PARENT[NNODE], is the parent node of each node.
//    The node with parent 0 is the root.
//
//    Output, int &DEPTH, the depth of the tree.
//
//    Output, int DEPTH_NODE[NNODE], the depth of each node.
//
{
  int i;
  int j;
  int root;
//
//  Find the root.
//
  root = -1;
  for ( i = 1; i <= nnode; i++ )
  {
    if ( parent[i-1] == 0 )
    {
      root = i;
      break;
    }
  }

  if ( root == -1 )
  {
    cerr << "\n";
    cerr << "TREE_ROOTED_DEPTH - Fatal error!\n";
    cerr << "  Could not find the root.\n";
    exit ( 1 );
  }
//
//  Determine the depth of each node by moving towards the node.
//  If you reach a node whose depth is already known, stop early.
//
  for ( i = 0; i < nnode; i++ )
  {
    depth_node[i] = 0;
  }

  for ( i = 1; i <= nnode; i++ )
  {
    j = i;

    while ( j != root )
    {
      depth_node[i-1] = depth_node[i-1] + 1;
      j = parent[j-1];

      if ( 0 < depth_node[j-1] )
      {
        depth_node[i-1] = depth_node[i-1] + depth_node[j-1];
        break;
      }
    }
  }
//
//  Determine the maximum depth.
//
  depth = i4vec_max ( nnode, depth_node );

  return;
}


//****************************************************************************80

int *tree_rooted_enum ( int nnode )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ROOTED_ENUM counts the number of unlabeled rooted trees.
//
//  Example:
//
//    Input    Output
//
//      1         1
//      2         1
//      3         2
//      4         4
//      5         9
//      6        20
//      7        48
//      8       115
//      9       286
//     10       719
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int NNODE, the number of nodes.
//
//    Output, int TREE_ROOTED_ENUM[NNODE].  NTREE(I) is the number of 
//    rooted, unlabeled trees on I nodes, for I = 1, 2, ... NNODE.
//
{
  int i;
  int id;
  int isum;
  int itd;
  int j;
  int nlast;
  int *ntree;

  ntree = new int[nnode];

  ntree[0] = 1;
 
  for ( nlast = 2; nlast <= nnode; nlast++ )
  {
    isum = 0;
 
    for ( id = 1; id <= nlast - 1; id++ )
    {
      i = nlast;
      itd = ntree[id-1] * id;
 
      for ( j = 1; j <= nlast - 1; j++ )
      {
        i = i - id;

        if ( i <= 0 )
        {
          break;
        }
        isum = isum + ntree[i-1] * itd;
      }
    }
    ntree[nlast-1] = isum / ( nlast - 1 );
  }
  return ntree;
}
//****************************************************************************80

int *tree_rooted_random ( int nnode, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    TREE_ROOTED_RANDOM selects a random unlabeled rooted tree.
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
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int NNODE, the number of nodes.
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Output, int TREE_ROOTED_RANDOM[NNODE].  (I,ITREE(I)) is the I-th edge
//    of the output tree for I = 2,NNODE.  ITREE(1)=0.
//
{
  int i;
  int id;
  int is1;
  int is2;
  int itd;
  int *itree;
  int iz;
  int j;
  int l;
  int ll;
  int ls;
  int m;
  int *ntree;
  int nval;
  double r;
  int *stack;

  if ( nnode <= 0  )
  {
    cerr << "\n";;
    cerr << "TREE_ROOTED_RANDOM - Fatal error!\n";
    cerr << "  NNODE = " << nnode << "\n";
    cerr << "  but NNODE must be at least 1.\n";
    exit ( 1 );
  }

  itree = new int[nnode];
  stack = new int[2*nnode];
//
//  Compute a table of the number of such trees for a given number of nodes.
//
  ntree = tree_rooted_enum ( nnode );
//
//  Now select one such tree at random.
//
  l = 0;

  nval = nnode;
  is1 = 0;
  is2 = 0;
  
  for ( ; ; )
  {
    while ( 2 < nval )
    {
      r = r8_uniform_01 ( seed );

      iz = ( int ) ( ( nval - 1 ) * ntree[nval-1] * r );

      id = 0;
  
      id = id + 1;
      itd = id * ntree[id-1];
      m = nval;
      j = 0;
 
      for ( ; ; )
      {
        j = j + 1;
        m = m - id;

        if ( m < 1 )
        {
          id = id + 1;
          itd = id * ntree[id-1];
          m = nval;
          j = 0;
          continue;
        }

        iz = iz - ntree[m-1] * itd;
        if ( iz < 0 )
        {
          break;
        }
      }
      is1 = is1 + 1;
      stack[0+(is1-1)*2] = j;
      stack[1+(is1-1)*2] = id;
      nval = m;
    }
 
    itree[is2] = l;
    l = is2 + 1;
    is2 = is2 + nval;

    if ( 1 < nval )
    {
      itree[is2-1] = is2 - 1;
    }
 
    for ( ; ; )
    {
      nval = stack[1+(is1-1)*2];
 
      if ( nval != 0 )
      {
        stack[1+(is1-1)*2] = 0;
        break;
      }
 
      j = stack[0+(is1-1)*2];
      is1 = is1 - 1;
      m = is2 - l + 1;
      ll = itree[l-1];
      ls = l + ( j - 1 ) * m - 1;
 
      if ( j != 1 )
      {
        for ( i = l; i <= ls; i++ )
        {
          itree[i+m-1] = itree[i-1] + m;
          if ( ( ( i - l ) % m ) == 0 )
          {
            itree[i+m-1] = ll;
          }
        }
      }
      is2 = ls + m;
 
      if ( is2 == nnode )
      {
        delete [] ntree;
        delete [] stack;
        return itree;
      }
      l = ll;
    }
  }
 
  delete [] ntree;
  delete [] stack;

  return itree;
}
//****************************************************************************80

void vec_next ( int n, int ibase, int iarray[], int &more )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_NEXT generates all N-vectors of integers modulo a given base.
//
//  Discussion:
//
//    The items are produced one at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Parameters:
//
//    Input, int N, the size of the vectors to be used.
//
//    Input, int IBASE, the base to be used.  IBASE = 2 will
//    give vectors of 0's and 1's, for instance.
//
//    Input/output, int IARRAY[N].  On each return,
//    IARRAY will contain entries in the range 0 to IBASE-1.
//
//    Input/output, int &MORE.  Set this variable 0 before
//    the first call.  Normally, MORE will be returned 1 but
//    once all the vectors have been generated, MORE will be
//    reset 0 and you should stop calling the program.
//
{
  int i;
  static int kount = 0;
  static int last = 0;
  int nn;

  if ( ! more )
  {
    kount = 1;
    last = i4_power ( ibase, n );
    more = 1;
    for ( i = 0; i < n; i++ )
    {
      iarray[i] = 0;
    }
  }
  else
  {
    kount = kount + 1;

    if ( kount == last )
    {
      more = 0;
    }

    iarray[n-1] = iarray[n-1] + 1;
 
    for ( i = 1; i <= n; i++ )
    {
      nn = n - i;

      if ( iarray[nn] < ibase )
      {
        return;
      }

      iarray[nn] = 0;

      if ( nn != 0 )
      {
        iarray[nn-1] = iarray[nn-1] + 1;
      }
    }
  }
  return;
}
//****************************************************************************80

void vec_random ( int n, int base, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_RANDOM selects a random N-vector of integers modulo a given base.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the vector to be generated.
//
//    Input, int BASE, the base to be used.
//
//    Input/output, int &SEED, a random number seed.
//
//    Output, int A[N], a list of N random values between
//    0 and BASE-1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = i4_uniform_ab ( 0, base-1, seed );
  }

  return;
}
//****************************************************************************80

int *vec_random_new ( int n, int base, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    VEC_RANDOM_NEW selects a random N-vector of integers modulo a given base.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the vector to be generated.
//
//    Input, int BASE, the base to be used.
//
//    Input/output, int &SEED, a random number seed.
//
//    Output, int VEC_RANDOM_NEW[N], a list of N random values between
//    0 and BASE-1.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i4_uniform_ab ( 0, base-1, seed );
  }

  return a;
}
