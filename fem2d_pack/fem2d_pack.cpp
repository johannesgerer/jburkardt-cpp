# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>
# include <cstring>

using namespace std;

# include "fem2d_pack.hpp"

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

void bandwidth_var ( int element_order, int element_num, int element_node[],
  int node_num, int var_node[], int var_num, int var[], int *ml, int *mu, 
  int *m )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH_VAR determines the bandwidth for finite element variables.
//
//  Discussion:
//
//    We assume that, attached to each node in the finite element mesh
//    there are a (possibly zero) number of finite element variables.
//    We wish to determine the bandwidth necessary to store the stiffness
//    matrix associated with these variables.
//
//    An entry K(I,J) of the stiffness matrix must be zero unless the
//    variables I and J correspond to nodes N(I) and N(J) which are
//    common to some element.
//
//    In order to determine the bandwidth of the stiffness matrix, we
//    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
//    is maximized.
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
//    We assume the finite element variable adjacency relationship is 
//    symmetric, so we are guaranteed that ML = MU.
//
//    Note that the user is free to number the variables in any way
//    whatsoever, and to associate variables to nodes in any way,
//    so that some nodes have no variables, some have one, and some
//    have several.  
//
//    The storage of the indices of the variables is fairly simple.
//    In VAR, simply list all the variables associated with node 1, 
//    then all those associated with node 2, and so on.  Then set up
//    the pointer array VAR_NODE so that we can jump to the section of
//    VAR where the list begins for any particular node.
//
//    The routine does not check that each variable is only associated
//    with a single node.  This would normally be the case in a finite
//    element setting.
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
  int node_global_i;
  int node_global_j;
  int node_local_i;
  int node_local_j;
  int var_global_i;
  int var_global_j;
  int var_local_i;
  int var_local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( node_local_i = 0; node_local_i < element_order; node_local_i++ )
    {
      node_global_i = element_node[node_local_i+element*element_order];

      for ( var_local_i = var_node[node_global_i-1]; 
            var_local_i <= var_node[node_global_i]-1; var_local_i++ )
      {
        var_global_i = var[var_local_i-1];

        for ( node_local_j = 0; node_local_j < element_order; node_local_j++ )
        {
          node_global_j = element_node[node_local_j+element*element_order];

          for ( var_local_j = var_node[node_global_j-1]; 
                var_local_j <= var_node[node_global_j]-1; var_local_j++ )
          {
            var_global_j = var[var_local_j-1];

            *mu = i4_max ( *mu, var_global_j - var_global_i );
            *ml = i4_max ( *ml, var_global_i - var_global_j );
          }
        }
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
//****************************************************************************80

void basis_11_t3 ( double t[2*3], int i, double p[2], double *qi, 
  double *dqidx, double *dqidy )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T3: one basis at one point for a T3 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the nodes of a triangle.
//
//           3
//          / \
//         /   \
//        /     \
//       1-------2
//
//    It evaluates the linear basis function Q(I)(X,Y) associated with
//    node I, which has the property that it is a linear function
//    which is 1 at node I and zero at the other two nodes.
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
//    Input, double T[2*3], the coordinates of the nodes.
//
//    Input, int I, the index of the desired basis function.
//    I should be between 1 and 3.
//
//    Input, double P[2], the coordinates of the point where 
//    the basis function is to be evaluated.
//
//    Output, double *QI, *DQIDX, *DQIDY, the value of the I-th basis function
//    and its X and Y derivatives.
//
{
  double area;
  int ip1;
  int ip2;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    cout << "\n";
    cout << "BASIS_11_T3 - Fatal error!\n";
    cout << "  Element has zero area.\n";
    cout << "  Area = " << area << "\n";
    exit ( 1 );
  }

  if ( i < 1 || 3 < i )
  {
    cout << "\n";
    cout << "BASIS_11_T3 - Fatal error!\n";
    cout << "  Basis index I is not between 1 and 3.\n";
    cout << "  I = " << i << "\n";
    exit ( 1 );
  }

  ip1 = i4_wrap ( i + 1, 1, 3 );
  ip2 = i4_wrap ( i + 2, 1, 3 );

  *qi = ( ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) 
        * ( p[1]           - t[1+(ip1-1)*2] ) 
        - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) 
        * ( p[0]           - t[0+(ip1-1)*2] ) ) / area;

  *dqidx = - ( t[1+(ip2-1)*2] - t[1+(ip1-1)*2] ) / area;
  *dqidy =   ( t[0+(ip2-1)*2] - t[0+(ip1-1)*2] ) / area;

  return;
}
//****************************************************************************80

void basis_11_t3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T3_TEST verifies BASIS_11_T3.
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
//    None
//
{
# define NODE_NUM 3

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  cout << "\n";
  cout << "BASIS_11_T3_TEST:\n";
  cout << "  Verify basis functions for element T3.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(10) << t[0+j*2] << "  "
         << setw(10) << t[1+j*2] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      cout << "  " << setw(10) << qj;
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t3 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    cout << "  "
         << setw(10) << sum_x << "  "
         << setw(10) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_11_t4 ( double t[2*4], int i, double p[], double *phi, 
  double *dphidx, double *dphidy )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T4: one basis at one point for a T4 element.
//
//  Discussion:
//
//    The T4 element is the cubic bubble triangle.
//
//    The routine is given the coordinates of the vertices of a triangle.
//    It works directly with these coordinates, and does not refer to a 
//    reference element.
//
//    The sides of the triangle DO NOT have to lie along a coordinate
//    axis.
//
//    The routine evaluates the basis functions associated with each vertex,
//    and their derivatives with respect to X and Y.
//
//  Physical Element T4: 
//       
//            3
//           / \
//          /   \
//         /  4  \
//        /       \
//       1---------2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*4], the coordinates of the vertices
//    of the triangle, and the coordinates of the centroid.  
//    It is common to list the first three points in counter clockwise
//    order.
//
//    Input, int I, the index of the basis function.
//
//    Input, double P[2], the point where the basis function
//    is to be evaluated.
//
//    Output, double *PHI, the value of the basis function
//    at the evaluation point.
//
//    Output, double *DPHIDX, *DPHIDY, the value of the 
//    derivatives at the evaluation point.
//
//  Local parameters:
//
//    Local, double AREA, is (twice) the area of the triangle.
//
{
  double area;
  double dpsidx[4];
  double dpsidy[4];
  int j;
  int n = 4;
  double psi[4];

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  psi[0] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1] - t[1+1*2] )     
                 - ( t[1+2*2] - t[1+1*2] ) * ( p[0] - t[0+1*2] ) );
  dpsidx[0] =    - ( t[1+2*2] - t[1+1*2] );
  dpsidy[0] =      ( t[0+2*2] - t[0+1*2] );

  psi[1] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1] - t[1+2*2] )     
                 - ( t[1+0*2] - t[1+2*2] ) * ( p[0] - t[0+2*2] ) );
  dpsidx[1] =    - ( t[1+0*2] - t[1+2*2] );
  dpsidy[1] =      ( t[0+0*2] - t[0+2*2] );

  psi[2] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1] - t[1+0*2] )     
                 - ( t[1+1*2] - t[1+0*2] ) * ( p[0] - t[0+0*2] ) );
  dpsidx[2] =    - ( t[1+1*2] - t[1+0*2] );
  dpsidy[2] =      ( t[0+1*2] - t[0+0*2] );
//
//  Normalize the first three functions.
//
  for ( j = 0; j < 3; j++ )
  {
    psi[j]    = psi[j]    / area;
    dpsidx[j] = dpsidx[j] / area;
    dpsidy[j] = dpsidy[j] / area;
  }
//
//  Compute the cubic bubble function.
//
  psi[3] = 27.0 * psi[0] * psi[1] * psi[2];

  dpsidx[3] = 27.0 * (
                  dpsidx[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidx[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidx[2] );

  dpsidy[3] = 27.0 * (
                  dpsidy[0] *    psi[1] *    psi[2]
                  +  psi[0] * dpsidy[1] *    psi[2]
                  +  psi[0] *    psi[1] * dpsidy[2] );
//
//  Subtract 1/3 of the cubic bubble function from each of the three linears.
//
  for ( j = 0; j < 3; j++ )
  {
       psi[j] =    psi[j] -    psi[3] / 3.0;
    dpsidx[j] = dpsidx[j] - dpsidx[3] / 3.0;
    dpsidy[j] = dpsidy[j] - dpsidy[3] / 3.0;
  }

  *phi = psi[i-1];
  *dphidx = dpsidx[i-1];
  *dphidy = dpsidy[i-1];

  return;
}
//****************************************************************************80

void basis_11_t4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T4_TEST verifies BASIS_11_T4.
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
//  Parameters:
//
//    None
//
{
# define NODE_NUM 4

  double dqjdx;
  double dqjdy;
  int i;
  int j;
  double qj;
  double p[2];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    0.0, 0.0 };
//
//  The node associated with the fourth basis function is the centroid.
//
  t[0+3*2] = ( t[0+0*2] + t[0+1*2] + t[0+2*2] ) / 3.0;
  t[1+3*2] = ( t[1+0*2] + t[1+1*2] + t[1+2*2] ) / 3.0;

  cout << "\n";
  cout << "BASIS_11_T4_TEST:\n";
  cout << "  Verify basis functions for element T4.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(10) << t[0+j*2] << "  "
         << setw(10) << t[1+j*2] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      cout << "  " << setw(10) << qj;
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";

  for ( i = 0; i < NODE_NUM; i++ )
  {
    p[0] = t[0+i*2];
    p[1] = t[1+i*2];

    sum_x = 0.0;
    sum_y = 0.0;
    for ( j = 0; j < NODE_NUM; j++ )
    {
      basis_11_t4 ( t, j+1, p, &qj, &dqjdx, &dqjdy );
      sum_x = sum_x + dqjdx;
      sum_y = sum_y + dqjdy;
    }
    cout << "  "
         << setw(10) << sum_x << "  "
         << setw(10) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_11_t6 ( double t[2*6], int i, double p[], double *bi, 
  double *dbidx, double *dbidy )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T6: one basis at one point for the T6 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the nodes of a triangle. 
//       
//           3
//          / \
//         6   5
//        /     \
//       1---4---2
//
//    It evaluates the quadratic basis function B(I)(X,Y) associated with
//    node I, which has the property that it is a quadratic function
//    which is 1 at node I and zero at the other five nodes.
//
//    This routine assumes that the sides of the triangle are straight,
//    so that the midside nodes fall on the line between two vertices.
//
//    This routine relies on the fact that each basis function can be
//    written as the product of two linear factors, which are easily
//    computed and normalized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the coordinates of the nodes.
//
//    Input, int I, the index of the desired basis function.
//    I should be between 1 and 6.
//
//    Input, double P[2], the coordinates of a point at which the basis
//    function is to be evaluated.
//
//    Output, double *BI, *DBIDX, *DBIDY, the values of the basis function
//    and its X and Y derivatives.
//
{
  double gf;
  double gn;
  double hf;
  double hn;
  int j1;
  int j2;
  int k1;
  int k2;

  if ( i < 1 || 6 < i )
  {
    cout << "\n";
    cout << "BASIS_11_T6 - Fatal error!\n";
    cout << "  Basis index I is not between 1 and 6.\n";
    cout << "  I = " << i << "\n";
    exit ( 1 );
  }
//
//  Determine the pairs of nodes.
//
  if ( i <= 3 )
  {
    j1 = i4_wrap ( i + 1, 1, 3 );
    j2 = i4_wrap ( i + 2, 1, 3 );
    k1 = i + 3;
    k2 = i4_wrap ( i + 5, 4, 6 );
  }
  else
  {
    j1 = i - 3;
    j2 = i4_wrap ( i - 3 + 2, 1, 3 );
    k1 = i4_wrap ( i - 3 + 1, 1, 3 );
    k2 = i4_wrap ( i - 3 + 2, 1, 3 );
  }
//
//  For C++ indexing, it is helpful to knock the indices down by one.
//
  i  = i  - 1;
  j1 = j1 - 1;
  j2 = j2 - 1;
  k1 = k1 - 1;
  k2 = k2 - 1;
//
//  Evaluate the two linear factors GF and HF, 
//  and their normalizers GN and HN.
//
  gf = ( p[0]      - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( p[1]      - t[1+j1*2] );

  gn = ( t[0+i*2]  - t[0+j1*2] ) * ( t[1+j2*2] - t[1+j1*2] ) 
     - ( t[0+j2*2] - t[0+j1*2] ) * ( t[1+i*2]  - t[1+j1*2] );

  hf = ( p[0]      - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( p[1]      - t[1+k1*2] );

  hn = ( t[0+i*2]  - t[0+k1*2] ) * ( t[1+k2*2] - t[1+k1*2] ) 
     - ( t[0+k2*2] - t[0+k1*2] ) * ( t[1+i*2]  - t[1+k1*2] );
//
//  Construct the basis function and its derivatives.
//
  *bi =        ( gf                      / gn ) 
             * ( hf                      / hn );

  *dbidx =   ( ( t[1+j2*2] - t[1+j1*2] ) / gn ) 
           * (   hf                      / hn )
           + (   gf                      / gn ) 
           * ( ( t[1+k2*2] - t[1+k1*2] ) / hn );

  *dbidy = - ( ( t[0+j2*2] - t[0+j1*2] ) / gn ) 
           * (   hf                      / hn )
           - (   gf                      / gn ) 
           * ( ( t[0+k2*2] - t[0+k1*2] ) / hn );

  return;
}
//****************************************************************************80

void basis_11_t6_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_11_T6_TEST verifies BASIS_11_T6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2006
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
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double p[2];
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };
  double v1;
  double v2;
  double v3;

  cout << "\n";
  cout << "BASIS_11_T6_TEST:\n";
  cout << "  Verify basis functions for element T6.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(6) << j << "  "
         << setw(7) << t[0+j*2] << "  "
         << setw(7) << t[1+j*2] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";
  for ( i = 1; i <= NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      p[0] = t[0+j*2];
      p[1] = t[1+j*2];

      basis_11_t6 ( t, i, p, &v1, &v2, &v3 );

      phi[i-1+j*NODE_NUM] = v1;
      dphidx[i-1+j*NODE_NUM] = v2;
      dphidy[i-1+j*NODE_NUM] = v3;
    }
  }
  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      cout << "  " << setw(7) << phi[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    cout << "  "
         << setw(14) << sum_x << "  "
         << setw(14) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_mn_q4 ( double q[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_Q4: all bases at N points for a Q4 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a quadrilateral.
//    It works directly with these coordinates, and does not refer to a 
//    reference element.
//
//    The sides of the element are presumed to lie along coordinate axes.
//
//    The routine evaluates the basis functions associated with each corner,
//    and their derivatives with respect to X and Y.
//
//  Physical Element Q4:
//
//    |
//    |  4-----3
//    |  |     |
//    |  |     |
//    Y  |     |
//    |  |     |
//    |  |     |
//    |  1-----2
//    |
//    +-----X------>
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
//  Parameters:
//
//    Input, double Q[2*4], the coordinates of the vertices.
//    It is common to list these points in counter clockwise order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the evaluation points.
//
//    Output, double PHI[4*N], the bases at the evaluation points.
//
//    Output, double DPHIDX[4*N], DPHIDY[4*N], the derivatives of the
//    bases at the evaluation points.
//
{
  double area;
  int i;
  int j;

  area =        ( q[0+2*2]         - q[0+0*2] ) 
              * ( q[1+2*2]         - q[1+0*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =      ( q[0+2*2] - p[0+j*2]             ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[1+j*4] =      (            p[0+j*2]  - q[0+0*2] ) 
                    * ( q[1+2*2] - p[1+j*2]             );
    phi[2+j*4] =      (            p[0+j*2]  - q[0+0*2] )
                    * (            p[1+j*2]  - q[1+0*2] );
    phi[3+j*4] =      ( q[0+2*2] - p[0+j*2]             )
                    * (            p[1+j*2]  - q[1+0*2] );

    dphidx[0+j*4] = - ( q[1+2*2] - p[1+j*2]             );
    dphidx[1+j*4] =   ( q[1+2*2] - p[1+j*2]             );
    dphidx[2+j*4] =   (            p[1+j*2]  - q[1+0*2] );
    dphidx[3+j*4] = - (            p[1+j*2]  - q[1+0*2] );

    dphidy[0+j*4] = - ( q[0+2*2] - p[0+j*2]             );
    dphidy[1+j*4] = - (            p[0+j*2]  - q[0+0*2] );
    dphidy[2+j*4] =   (            p[0+j*2]  - q[0+0*2] );
    dphidy[3+j*4] =   ( q[0+2*2] - p[0+j*2]             );
  }
//
//  Normalize.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 4; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
  return;
}
//****************************************************************************80

void basis_mn_q4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_Q4_TEST verifies BASIS_MN_Q4.
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
//  Parameters:
//
//    None
//
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double q[2*NODE_NUM] = {
    3.0, 1.0,
    5.0, 1.0,
    5.0, 4.0, 
    3.0, 4.0 };
  double sum_x;
  double sum_y;

  cout << "\n";
  cout << "BASIS_MN_Q4_TEST:\n";
  cout << "  Verify basis functions for element Q4.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( i = 0; i < NODE_NUM; i++ )
  {
    cout << "  "
         << setw(10) << q[0+i*2] << "  "
         << setw(10) << q[1+i*2] << "\n";
  }
 
  basis_mn_q4 ( q, NODE_NUM, q, phi, dphidx, dphidy );

  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
    cout << "  " << setw(10) << phi[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    cout << "  "
         << setw(10) << sum_x << "  "
         << setw(10) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_mn_t3 ( double t[2*3], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T3: all bases at N points for a T3 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a triangle.
//    It works directly with these coordinates, and does not refer to a 
//    reference element.
//
//    The sides of the triangle DO NOT have to lie along a coordinate
//    axis.
//
//    The routine evaluates the basis functions associated with each vertex,
//    and their derivatives with respect to X and Y.
//
//  Physical Element T3: 
//       
//            3
//           / \
//          /   \
//         /     \
//        /       \
//       1---------2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices
//    of the triangle.  It is common to list these points in counter clockwise
//    order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the points where the basis functions 
//    are to be evaluated.
//
//    Output, double PHI[3*N], the value of the basis functions 
//    at the evaluation points.
//
//    Output, double DPHIDX[3*N], DPHIDY[3*N], the value of the 
//    derivatives at the evaluation points.
//
//  Local parameters:
//
//    Local, double AREA, is (twice) the area of the triangle.
//
{
  double area;
  int i;
  int j;
  double temp;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  if ( area == 0.0 )
  {
    cout << "\n";
    cout << "BASIS_MN_T3 - Fatal error!\n";
    cout << "  Element has zero area.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*3] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*3] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*3] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*3] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*3] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*3] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*3] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*3] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*3] =      ( t[0+1*2] - t[0+0*2] );
  }
//
//  Normalize.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*3]    = phi[i+j*3]    / area;
      dphidx[i+j*3] = dphidx[i+j*3] / area;
      dphidy[i+j*3] = dphidy[i+j*3] / area;
    }
  }

  return;
}
//****************************************************************************80

void basis_mn_t3_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T3_TEST verifies BASIS_MN_T3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None.
//
{
# define NODE_NUM 3

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0 };

  cout << "\n";
  cout << "BASIS_MN_T3_TEST:\n";
  cout << "  Verify basis functions for element T3.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(10) << t[0+j*2] << "  "
         << setw(10) << t[1+j*2] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  basis_mn_t3 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      cout << "  " << setw(10) << phi[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    cout << "  "
         << setw(10) << sum_x << "  "
         << setw(10) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_mn_t4 ( double t[2*4], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T4: all bases at N points for a T4 element.
//
//  Discussion:
//
//    The T4 element is the cubic bubble triangle.
//
//    The routine is given the coordinates of the vertices of a triangle.
//    It works directly with these coordinates, and does not refer to a 
//    reference element.
//
//    The sides of the triangle DO NOT have to lie along a coordinate
//    axis.
//
//    The routine evaluates the basis functions associated with each vertex,
//    and their derivatives with respect to X and Y.
//
//  Physical Element T4: 
//       
//            3
//           / \
//          /   \
//         /  4  \
//        /       \
//       1---------2
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
//  Parameters:
//
//    Input, double T[2*4], the coordinates of the vertices
//    of the triangle, and the coordinates of the centroid.  
//    It is common to list the first three points in counter clockwise
//    order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the points where the basis functions 
//    are to be evaluated.
//
//    Output, double PHI[4*N], the value of the basis functions 
//    at the evaluation points.
//
//    Output, double DPHIDX[4*N], DPHIDY[4*N], the value of the 
//    derivatives at the evaluation points.
//
//  Local parameters:
//
//    Local, double AREA, is (twice) the area of the triangle.
//
{
  double area;
  int i;
  int j;

  area = t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) 
       + t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) 
       + t[0+2*2] * ( t[1+0*2] - t[1+1*2] );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] =     (   ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] )     
                       - ( t[1+2*2] - t[1+1*2] ) * ( p[0+j*2] - t[0+1*2] ) );
    dphidx[0+j*4] =    - ( t[1+2*2] - t[1+1*2] );
    dphidy[0+j*4] =      ( t[0+2*2] - t[0+1*2] );

    phi[1+j*4] =     (   ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] )     
                       - ( t[1+0*2] - t[1+2*2] ) * ( p[0+j*2] - t[0+2*2] ) );
    dphidx[1+j*4] =    - ( t[1+0*2] - t[1+2*2] );
    dphidy[1+j*4] =      ( t[0+0*2] - t[0+2*2] );

    phi[2+j*4] =     (   ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] )     
                       - ( t[1+1*2] - t[1+0*2] ) * ( p[0+j*2] - t[0+0*2] ) );
    dphidx[2+j*4] =    - ( t[1+1*2] - t[1+0*2] );
    dphidy[2+j*4] =      ( t[0+1*2] - t[0+0*2] );
  }
//
//  Normalize the first three functions.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phi[i+j*4]    = phi[i+j*4]    / area;
      dphidx[i+j*4] = dphidx[i+j*4] / area;
      dphidy[i+j*4] = dphidy[i+j*4] / area;
    }
  }
//
//  Compute the cubic bubble function.
//
  for ( j = 0; j < n; j++ )
  {
       phi[3+j*4] = 27.0 * phi[0+j*4] * phi[1+j*4] * phi[2+j*4];

    dphidx[3+j*4] = 27.0 * (
                    dphidx[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidx[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidx[2+j*4] );

    dphidy[3+j*4] = 27.0 * (
                    dphidy[0+j*4] *    phi[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] * dphidy[1+j*4] *    phi[2+j*4]
                    +  phi[0+j*4] *    phi[1+j*4] * dphidy[2+j*4] );
  }
//
//  Subtract 1/3 of the cubic bubble function from each of the three linears.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
         phi[i+j*4] =    phi[i+j*4] -    phi[3+j*4] / 3.0;
      dphidx[i+j*4] = dphidx[i+j*4] - dphidx[3+j*4] / 3.0;
      dphidy[i+j*4] = dphidy[i+j*4] - dphidy[3+j*4] / 3.0;
    }
  }

  return;
}
//****************************************************************************80

void basis_mn_t4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T4_TEST verifies BASIS_MN_T4.
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
//  Parameters:
//
{
# define NODE_NUM 4

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = { 
    2.0, 0.0,
    4.0, 2.0,
    0.0, 4.0,
    2.0, 2.0 };

  cout << "\n";
  cout << "BASIS_MN_T4_TEST:\n";
  cout << "  Verify basis functions for element T4.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(10) << t[0+j*2] << "  "
         << setw(10) << t[1+j*2] << "\n";
  }

  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  basis_mn_t4 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      cout << "  " << setw(10) << phi[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";

  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    cout << "  "
         << setw(10) << sum_x << "  "
         << setw(10) << sum_y << "\n";
  }
  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_mn_t6 ( double t[2*6], int n, double p[], double phi[], 
  double dphidx[], double dphidy[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T6: all bases at N points for a T6 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices and midside
//    nodes of a triangle.  It works directly with these coordinates, and does 
//    not refer to a reference element.
//
//    This routine requires that the midside nodes be "in line"
//    with the vertices, that is, that the sides of the triangle be
//    straight.  However, the midside nodes do not actually have to
//    be halfway along the side of the triangle.  
//
//  The physical element T6:
//
//    This picture indicates the assumed ordering of the six nodes
//    of the triangle.
//
//    |
//    |   
//    |        3
//    |       / \
//    |      /   \
//    Y     6     5
//    |    /       \
//    |   /         \
//    |  1-----4-----2
//    |
//    +--------X-------->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the nodal oordinates of the element.
//    It is common to list these points in counter clockwise order.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[2*N], the coordinates of the point where
//    the basis functions are to be evaluated.
//
//    Output, double PHI[6*N], the value of the basis functions at P.
//
//    Output, double DPHIDX[6*N], DPHIDY[6*N], the value of the X 
//    and Y derivatives of the basis functions at P.
//
{
  double gn;
  double gx;
  double hn;
  double hx;
  int j;

  for ( j = 0; j < n; j++ )
  {
//
//  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( p[1+j*2] - t[1+3*2] );

    hn = ( t[0+0*2] - t[0+3*2] ) * ( t[1+5*2] - t[1+3*2] ) 
       - ( t[0+5*2] - t[0+3*2] ) * ( t[1+0*2] - t[1+3*2] );

    phi[0+j*6] =     ( gx * hx ) / ( gn * hn );
    dphidx[0+j*6] =  (      ( t[1+2*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+5*2] - t[1+3*2] ) ) / ( gn * hn );
    dphidy[0+j*6] = -(      ( t[0+2*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+5*2] - t[0+3*2] ) ) / ( gn * hn );
//
//  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
//
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+1*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( p[1+j*2] - t[1+4*2] );

    hn = ( t[0+1*2] - t[0+4*2] ) * ( t[1+3*2] - t[1+4*2] ) 
       - ( t[0+3*2] - t[0+4*2] ) * ( t[1+1*2] - t[1+4*2] );

    phi[1+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[1+j*6] =  (      ( t[1+2*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+3*2] - t[1+4*2] ) ) / ( gn * hn );
    dphidy[1+j*6] = -(      ( t[0+2*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+3*2] - t[0+4*2] ) ) / ( gn * hn );
//
//  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+2*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( p[1+j*2] - t[1+5*2] );

    hn = ( t[0+2*2] - t[0+5*2] ) * ( t[1+4*2] - t[1+5*2] ) 
       - ( t[0+4*2] - t[0+5*2] ) * ( t[1+2*2] - t[1+5*2] );

    phi[2+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[2+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+4*2] - t[1+5*2] ) ) / ( gn * hn );
    dphidy[2+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+4*2] - t[0+5*2] ) ) / ( gn * hn );
//
//  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
//
    gx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    gn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] ) 
       - ( t[0+0*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    hx = ( p[0+j*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( p[1+j*2] - t[1+2*2] );

    hn = ( t[0+3*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
       - ( t[0+1*2] - t[0+2*2] ) * ( t[1+3*2] - t[1+2*2] );

    phi[3+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[3+j*6] =  (      ( t[1+0*2] - t[1+2*2] ) * hx 
                     + gx * ( t[1+1*2] - t[1+2*2] ) ) / ( gn * hn );
    dphidy[3+j*6] = -(      ( t[0+0*2] - t[0+2*2] ) * hx 
                     + gx * ( t[0+1*2] - t[0+2*2] ) ) / ( gn * hn );
//
//  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
//
    gx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    gn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) 
       - ( t[0+1*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    hx = ( p[0+j*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( p[1+j*2] - t[1+0*2] );

    hn = ( t[0+4*2] - t[0+0*2] ) * ( t[1+2*2] - t[1+0*2] ) 
       - ( t[0+2*2] - t[0+0*2] ) * ( t[1+4*2] - t[1+0*2] );

    phi[4+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[4+j*6] =  (      ( t[1+1*2] - t[1+0*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+0*2] ) ) / ( gn * hn );
    dphidy[4+j*6] = -(      ( t[0+1*2] - t[0+0*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+0*2] ) ) / ( gn * hn );
//
//  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
//
    gx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    gn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+0*2] - t[1+1*2] ) 
       - ( t[0+0*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    hx = ( p[0+j*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( p[1+j*2] - t[1+1*2] );

    hn = ( t[0+5*2] - t[0+1*2] ) * ( t[1+2*2] - t[1+1*2] ) 
       - ( t[0+2*2] - t[0+1*2] ) * ( t[1+5*2] - t[1+1*2] );

    phi[5+j*6] = ( gx * hx ) / ( gn * hn );
    dphidx[5+j*6] =  (      ( t[1+0*2] - t[1+1*2] ) * hx 
                     + gx * ( t[1+2*2] - t[1+1*2] ) ) / ( gn * hn );
    dphidy[5+j*6] = -(      ( t[0+0*2] - t[0+1*2] ) * hx 
                     + gx * ( t[0+2*2] - t[0+1*2] ) ) / ( gn * hn );
  }

  return;
}
//****************************************************************************80

void basis_mn_t6_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T6_TEST verifies BASIS_MN_T6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2006
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
# define NODE_NUM 6

  double dphidx[NODE_NUM*NODE_NUM];
  double dphidy[NODE_NUM*NODE_NUM];
  int i;
  int j;
  double phi[NODE_NUM*NODE_NUM];
  double sum_x;
  double sum_y;
  double t[2*NODE_NUM] = {
    2.0, 0.0,
    4.0, 3.0,
    0.0, 4.0,
    3.0, 1.5,
    2.0, 3.5,
    1.0, 2.0 };

  cout << "\n";
  cout << "BASIS_MN_T6_TEST:\n";
  cout << "  Verify basis functions for element T6.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  cout << "\n";
  cout << "  Physical Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(6) << j << "  "
         << setw(7) << t[0+j*2] << "  "
         << setw(7) << t[1+j*2] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  basis_mn_t6 ( t, NODE_NUM, t, phi, dphidx, dphidy );

  for ( i = 0; i < NODE_NUM; i++ )
  {
    for ( j = 0; j < NODE_NUM; j++ )
    {
      cout << "  " << setw(7) << phi[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The X and Y derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dPhidX sum, dPhidY sum:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    sum_x = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_x = sum_x + dphidx[i+j*NODE_NUM];
    }
    sum_y = 0.0;
    for ( i = 0; i < NODE_NUM; i++ )
    {
      sum_y = sum_y + dphidy[i+j*NODE_NUM];
    }
    cout << "  "
         << setw(14) << sum_x << "  "
         << setw(14) << sum_y << "\n";
  }

  return;
# undef NODE_NUM
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

double degrees_to_radians ( double angle )

//****************************************************************************80
//
//  Purpose: 
//
//    DEGREES_TO_RADIANS converts an angle from degrees to radians.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double ANGLE, an angle in degrees.
//
//    Output, double DEGREES_TO_RADIANS, the equivalent angle
//    in radians.
//
{
# define PI 3.141592653589793

  return ( angle * PI / 180.0 );

# undef PI
}
//****************************************************************************80

void derivative_average_t3 ( int node_num, double node_xy[], int element_num,
  int element_node[], double c[], double dcdx[], double dcdy[] )

//****************************************************************************80
//
//  Purpose:
//
//    DERIVATIVE_AVERAGE_T3 averages derivatives at the nodes of a T3 mesh.
//
//  Discussion:
//
//    This routine can be used to compute an averaged nodal value of any
//    quantity associated with the finite element function.  At a node
//    that is shared by several elements, the fundamental function
//    U will be continuous, but its spatial derivatives, for instance,
//    will generally be discontinuous.  This routine computes the
//    value of the spatial derivatives in each element, and averages
//    them, to make a reasonable assignment of a nodal value.
//
//    Note that the ELEMENT_NODE array is assumed to be 1-based, rather
//    than 0-based.  Thus, entries from this array must be decreased by
//    1 before being used!
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
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[3*ELEMENT_NUM], the element->node data.
//
//    Input, double C[NODE_NUM], the finite element coefficient vector.
//
//    Output, double DCDX[NODE_NUM], DCDY[NODE_NUM], the averaged
//    values of dCdX and dCdY at the nodes.
//
{
# define OFFSET 1

  int dim;
  double dphidx[3*3];
  double dphidy[3*3];
  int element;
  int j;
  int node;
  int node_count[node_num];
  int node_global1;
  int node_global2;
  int node_local1;
  int node_local2;
  double phi[3*3];
  double t[2*3];

  for ( node = 0; node < node_num; node++ )
  {
    node_count[node] = 0;
    dcdx[node] = 0.0;
    dcdy[node] = 0.0;
  }
//
//  Consider every element.
//
  for ( element = 0; element < element_num; element++ )
  {
//
//  Get the coordinates of the nodes of the element.
//
    for ( j = 0; j < 3; j++ )
    {
      for ( dim = 0; dim < 2; dim++ )
      {
        t[dim+2*j] = node_xy[dim+(element_node[j+element*3]-OFFSET)];
      }
    }
//
//  Evaluate the X and Y derivatives of the 3 basis functions at the
//  3 nodes.
//
    basis_mn_t3 ( t, 3, t, phi, dphidx, dphidy );
//
//  Evaluate dCdX and dCdY at each node in the element, and add
//  them to the running totals.
//
    for ( node_local1 = 0; node_local1 < 3; node_local1++ )
    {
      node_global1 = element_node[node_local1+element*3]-OFFSET;

      for ( node_local2 = 0; node_local2 < 3; node_local2++ )
      {
        node_global2 = element_node[node_local2+element*3]-OFFSET;

        dcdx[node_global1] = dcdx[node_global1]
          + c[node_global2] * dphidx[node_local2+node_local1*3];

        dcdy[node_global1] = dcdy[node_global1] 
          + c[node_global2] * dphidy[node_local2+node_local1*3];
      }
      node_count[node_global1] = node_count[node_global1] + 1;
    }
  }
//
//  Average the running totals.
//
  for ( node = 0; node < node_num; node++ )
  {
    dcdx[node] = dcdx[node] / ( double ) node_count[node];
    dcdy[node] = dcdy[node] / ( double ) node_count[node];
  }
  return;
# undef OFFSET
}
//****************************************************************************80

void div_q4 ( int m, int n, double u[], double v[], double xlo, double xhi, 
  double ylo, double yhi, double div[], double vort[] )

//****************************************************************************80
//
//  Purpose: 
//
//    DIV_Q4 estimates the divergence and vorticity of a discrete field.
//
//  Discussion:
//
//    The routine is given the values of a vector field ( U(X,Y), V(X,Y) ) at
//    an array of points ( X(1:M), Y(1:N) ).
//
//    The routine models the vector field over the interior of this region using
//    a bilinear interpolant.  It then uses the interpolant to estimate the
//    value of the divergence:
//
//      DIV(X,Y) = dU/dX + dV/dY
//
//    and the vorticity:
//
//      VORT(X,Y) = dV/dX - dU/dY
//
//    at the center point of each of the bilinear elements.
//
//        |       |       |
//      (3,1)---(3,2)---(3,3)---
//        |       |       |
//        | [2,1] | [2,2] |
//        |       |       |
//      (2,1)---(2,2)---(2,3)---
//        |       |       |
//        | [1,1] | [1,2] |
//        |       |       |
//      (1,1)---(1,2)---(1,3)---
//
//    Here, the nodes labeled with parentheses represent the points at
//    which the original (U,V) data is given, while the nodes labeled
//    with square brackets represent the centers of the bilinear
//    elements, where the approximations to the divergence and vorticity
//    are made.
//
//    The reason for evaluating the divergence and vorticity in this way
//    is that the bilinear interpolant to the (U,V) data is not
//    differentiable at the boundaries of the elements, nor especially at
//    the nodes, but is an (infinitely differentiable) bilinear function
//    in the interior of each element.  If a value at the original nodes
//    is strongly desired, then the average at the four surrounding
//    central nodes may be taken.
//
//  Element Q4:
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
//    02 February 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of data rows.  M must be at least 2.
//
//    Input, int N, the number of data columns.  N must be at least 2.
//
//    Input, double U[M*N], V[M*N], the value of the components 
//    of a vector quantity whose divergence and vorticity are desired. 
//    A common example would be that U and V are the horizontal and 
//    vertical velocity component of a flow field.
//
//    Input, double XLO, XHI, the minimum and maximum X coordinates.
//
//    Input, double YLO, YHI, the minimum and maximum Y coordinates.
//
//    Output, double DIV[(M-1)*(N-1)], an estimate for
//    the divergence in the bilinear element that lies between
//    data rows I and I+1, and data columns J and J+1.
//
//    Output, double VORT[(M-1)*(N-1)], an estimate for
//    the vorticity in the bilinear element that lies between
//    data rows I and I+1, and data columns J and J+1.
//
{
  double dphidx[4];
  double dphidy[4];
  int i;
  int j;
  double p[2];
  double phi[4];
  double q[2*4];
  double xl;
  double xr;
  double yb;
  double yt;

  if ( m <= 1 )
  {
    cout << "\n";
    cout << "DIV_Q4 - Fatal error!\n";
    cout << "  M must be at least 2,\n";
    cout << "  but the input value of M is " << m << "\n";
    exit ( 1 );
  }

  if ( n <= 1 )
  {
    cout << "\n";
    cout << "DIV_Q4 - Fatal error!\n";
    cout << "  N must be at least 2,\n";
    cout << "  but the input value of N is " << n << "\n";
    exit ( 1 );
  }

  if ( xhi == xlo )
  {
    cout << "\n";
    cout << "DIV_Q4 - Fatal error!\n";
    cout << "  XHI and XLO must be distinct,\n";
    cout << "  but the input value of XLO is " << xlo << "\n";
    cout << "  and the input value of XHI is " << xhi << "\n";
    exit ( 1 );
  }

  if ( yhi == ylo )
  {
    cout << "\n";
    cout << "DIV_Q4 - Fatal error!\n";
    cout << "  YHI and YLO must be distinct,\n";
    cout << "  but the input value of YLO is " << ylo << "\n";
    cout << "  and the input value of YHI is " << yhi << "\n";
    exit ( 1 );
  }

  for ( i = 1; i <= m-1; i++ )
  {
    yb = ( ( double ) ( 2 * m - 2 * i     ) * ylo   
         + ( double ) (         2 * i - 2 ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );
    p[1] = ( ( double ) ( 2 * m - 2 * i - 1 ) * ylo   
           + ( double ) (         2 * i - 1 ) * yhi ) 
           / ( double ) ( 2 * m         - 2 );
    yt = ( ( double ) ( 2 * m - 2 * i - 2 ) * ylo   
         + ( double ) (         2 * i     ) * yhi ) 
         / ( double ) ( 2 * m         - 2 );

    q[1+0*2] = yb;
    q[1+1*2] = yb;
    q[1+2*2] = yt;
    q[1+3*2] = yt;

    for ( j = 1; j <= n-1; j++ )
    {
      xl = ( ( double ) ( 2 * n - 2 * j     ) * xlo   
           + ( double ) (         2 * j - 2 ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );
      p[0] = ( ( double ) ( 2 * n - 2 * j - 1 ) * xlo   
             + ( double ) (         2 * j - 1 ) * xhi ) 
             / ( double ) ( 2 * n         - 2 );
      xr = ( ( double ) ( 2 * n - 2 * j - 2 ) * xlo   
           + ( double ) (         2 * j     ) * xhi ) 
           / ( double ) ( 2 * n         - 2 );

      q[0+0*2] = xl;
      q[0+1*2] = xr;
      q[0+2*2] = xr;
      q[0+3*2] = xl;
//
//  Evaluate the basis function and derivatives at the center of the element.
//
      basis_mn_q4 ( q, 1, p, phi, dphidx, dphidy );
//
//  Note the following formula for the value of U and V at the same
//  point that the divergence and vorticity are being evaluated.
//
//         umid =  u(i  ,j  ) * phi[0] &
//               + u(i  ,j+1) * phi[1] &
//               + u(i+1,j+1) * phi[2] &
//               + u(i+1,j  ) * phi[3] 
//
//         vmid =  v(i  ,j  ) * phi[0] &
//               + v(i  ,j+1) * phi[1] &
//               + v(i+1,j+1) * phi[2] &
//               + v(i+1,j  ) * phi[3] 
//
      div[i-1+(j-1)*(m-1)] =  
                  u[i-1+(j-1)*m] * dphidx[0] + v[i-1+(j-1)*m] * dphidy[0] 
                + u[i-1+(j  )*m] * dphidx[1] + v[i-1+(j  )*m] * dphidy[1] 
                + u[i  +(j  )*m] * dphidx[2] + v[i  +(j  )*m] * dphidy[2] 
                + u[i  +(j-1)*m] * dphidx[3] + v[i  +(j-1)*m] * dphidy[3];

      vort[i-1+(j-1)*(m-1)] =  
                   v[i-1+(j-1)*m] * dphidx[0] - u[i-1+(j-1)*m] * dphidy[0] 
                 + v[i-1+(j  )*m] * dphidx[1] - u[i-1+(j  )*m] * dphidy[1] 
                 + v[i  +(j  )*m] * dphidx[2] - u[i  +(j  )*m] * dphidy[2] 
                 + v[i  +(j-1)*m] * dphidx[3] - u[i  +(j-1)*m] * dphidy[3]; 
    }
  }

  return;
}
//****************************************************************************80

string element_code ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    ELEMENT_CODE returns the code for each element.
//
//  List:
//
//    I  ELEMENT_CODE   Definition
//    -  ------------   ----------
//    1  Q4             4 node linear Lagrange/serendipity quadrilateral;
//    2  Q8             8 node quadratic serendipity quadrilateral;
//    3  Q9             9 node quadratic Lagrange quadrilateral;
//    4  Q12            12 node cubic serendipity quadrilateral;
//    5  Q16            16 node cubic Lagrange quadrilateral;
//    6  QL             6 node linear/quadratic quadrilateral;
//    7  T3             3 node linear triangle;
//    8  T4             4 node cubic bubble triangle
//    9  T6             6 node quadratic triangle;
//   10  T10            10 node cubic triangle.
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
//  Parameters:
//
//    Input, int I, the number of the element.  
//
//    Output, string ELEMENT_CODE, the code for the element.
//
{
  string value;

  if ( i == 1 )
  {
    value = "Q4";
  }
  else if ( i == 2 )
  {
    value = "Q8";
  }
  else if ( i == 3 )
  {
    value = "Q9";
  }
  else if ( i == 4 )
  {
    value = "Q12";
  }
  else if ( i == 5 )
  {
    value = "Q16";
  }
  else if ( i == 6 )
  {
    value = "QL";
  }
  else if ( i == 7 )
  {
    value = "T3";
  }
  else if ( i == 8 )
  {
    value = "T4";
  }
  else if ( i == 9 )
  {
    value = "T6";
  }
  else if ( i == 10 )
  {
    value = "T10";
  }
  else
  {
    value = "????";
  }

  return value;
}
//****************************************************************************80

void elements_eps ( string file_name, int node_num, double node_xy[], string code, 
  int element_num, bool element_mask[], int element_node[], int node_show, 
  int element_show )

//****************************************************************************80
//
//  Purpose: 
//
//    ELEMENTS_EPS creates an EPS file image of the elements of a grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, string CODE, the code for the element.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, bool ELEMENT_MASK[ELEMENT_NUM], a mask for the elements.
//    Only elements with a TRUE mask will be shown.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes making up
//    each element.
//
//    Input, int NODE_SHOW:
//    0, do not show nodes;
//    1, show nodes;
//    2, show nodes and label them.
//
//    Input, int TRIANGLE_SHOW:
//    0, do not show triangles;
//    1, show triangles;
//    2, show triangles and label them.
//
{
  double ave_x;
  double ave_y;
  int circle_size = 3;
  int delta;
  int e;
  int element;
  int element_order;
  ofstream file_unit;
  int i;
  int j;
  int local;
  int node;
  bool *node_mask;
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

  element_order = order_code ( code );
//
//  Determine which nodes are visible, controlled by which elements are visible.
//
  node_mask = new bool[node_num];
  for ( node = 0; node < node_num; node++ )
  {
    node_mask[node] = false;
  }
  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element*element_order]-1;
        node_mask[node] = true;
      }
    }
  }
//
//  We need to do some figuring here, so that we can determine
//  the range of the data, and hence the height and width
//  of the piece of paper.
//
  x_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
     if ( node_mask[node] )
     {
       if ( x_max < node_xy[0+node*2] )
       {
         x_max = node_xy[0+node*2];
       }
    }
  }
  x_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
       if ( node_xy[0+node*2] < x_min )
       {
         x_min = node_xy[0+node*2];
       }
    }
  }
  x_scale = x_max - x_min;

  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = -r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( y_max < node_xy[1+node*2] )
      {
        y_max = node_xy[1+node*2];
      }
    }
  }
  y_min = r8_huge ( );
  for ( node = 0; node < node_num; node++ )
  {
    if ( node_mask[node] )
    {
      if ( node_xy[1+node*2] < y_min )
      {
        y_min = node_xy[1+node*2];
      }
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

  file_unit.open ( file_name.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "ELEMENTS_EPS - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: elements_eps.C\n";
  file_unit << "%%Title: " << file_name << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  if ( 1 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Draw filled dots at each node:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.150  0.750  setrgbcolor\n";
    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    {
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        file_unit << "newpath  " 
          << x_ps << "  " 
          << y_ps << "  "
          << circle_size << " 0 360 arc closepath fill\n";
      }
    }
  }
//
//  Label the nodes.
//
  if ( 2 <= node_show )
  {
    file_unit << "%\n";
    file_unit << "%  Label the nodes:\n";
    file_unit << "%\n";
    file_unit << "%  Set the color to darker blue:\n";
    file_unit << "%\n";
    file_unit << "0.000  0.250  0.850  setrgbcolor\n";
    file_unit << "/Times-Roman findfont\n";
    file_unit << "0.20 inch scalefont\n";
    file_unit << "setfont\n";

    file_unit << "%\n";

    for ( node = 0; node < node_num; node++ )
    { 
      if ( node_mask[node] )
      {
        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        file_unit << "newpath  " 
          << x_ps     << "  " 
          << y_ps + 5 << "  moveto ("
          << node+1   << ") show\n";
      }
    }
  }
//
//  Draw the elements.
//
  file_unit << "%\n";
  file_unit << "%  Draw the element sides:\n";
  file_unit << "%\n";
  file_unit << " 9.0000 0.0000 0.0000 setrgbcolor\n";

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      local = 1;
      node = element_node[local-1+element*element_order] - 1;

      x_ps = ( int ) (
        ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
        + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
        / ( x_max                     - x_min ) );

      y_ps = ( int ) (
        ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
        + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
        / ( y_max                     - y_min ) );

      file_unit << "newpath  " 
        << x_ps << "  " 
        << y_ps << "  moveto\n";

      for ( ; ; )
      {
        local = next_boundary_node ( local, code );
        node = element_node[local-1+element*element_order] - 1;

        x_ps = ( int ) (
          ( ( x_max - node_xy[0+node*2]         ) * ( double ) ( x_ps_min )  
          + (       + node_xy[0+node*2] - x_min ) * ( double ) ( x_ps_max ) ) 
          / ( x_max                     - x_min ) );

        y_ps = ( int ) (
          ( ( y_max - node_xy[1+node*2]         ) * ( double ) ( y_ps_min )  
          + (         node_xy[1+node*2] - y_min ) * ( double ) ( y_ps_max ) ) 
          / ( y_max                     - y_min ) );

        file_unit << "  " 
          << x_ps << "  " 
          << y_ps << "  lineto\n";

        if ( local == 1 )
        {
          break;
        }
      }
      file_unit << "stroke\n";
    }
  }
//
//  Label the elements.
//
  file_unit << "%\n";
  file_unit << "%  Label the elements:\n";
  file_unit << "%\n";
  file_unit << " 1.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.30 inch scalefont setfont\n";

  for ( element = 0; element < element_num; element++ )
  {
    if ( element_mask[element] )
    {
      ave_x = 0.0;
      ave_y = 0.0;

      for ( local = 0; local < element_order; local++ )
      {
        node = element_node[local+element_order*element] - 1;

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

      file_unit << "newpath  " 
        << x_ps     << "  " 
        << y_ps + 5 << "  moveto ("
        << element+1   << ") show\n";
    }
  }
//
//  Finish up.
//
  file_unit << "%\n";
  file_unit << "restore  showpage\n";
  file_unit << "%\n";
  file_unit << "%  End of page.\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  delete [] node_mask;

  return;
}
//****************************************************************************80

int *grid_element ( string code, int element_order, int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_ELEMENT returns the element grid associated with any available element.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
//    "T4", "T6" and "T10".
//
//    Input, int ELEMENT_ORDER, the number of nodes per element.
//
//    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY for quadrilaterals, or 2 * NELEMX * NELEMY for
//    triangles.
//
//    Output, int GRID_ELEMENT[ELEMENT_ORDER*ELEMENT_NUM], the nodes 
//    that form each element.  
//
{
  int *element_node;

  if ( code == "Q4" )
  {
    element_node = grid_q4_element ( nelemx, nelemy );
  }
  else if ( code == "Q8" )
  {
    element_node = grid_q8_element ( nelemx, nelemy );
  }
  else if ( code == "Q9" )
  {
    element_node = grid_q9_element ( nelemx, nelemy );
  }
  else if ( code == "Q12" )
  {
    element_node = grid_q12_element ( nelemx, nelemy );
  }
  else if ( code == "Q16" )
  {
    element_node = grid_q16_element ( nelemx, nelemy );
  }
  else if ( code == "QL" )
  {
    element_node = grid_ql_element ( nelemx, nelemy );
  }
  else if ( code == "T3" )
  {
    element_node = grid_t3_element ( nelemx, nelemy );
  }
  else if ( code == "T4" )
  {
    element_node = grid_t4_element ( nelemx, nelemy );
  }
  else if ( code == "T6" )
  {
    element_node = grid_t6_element ( nelemx, nelemy );
  }
  else if ( code == "T10" )
  {
    element_node = grid_t10_element ( nelemx, nelemy );
  }
  else
  {
    element_node = NULL;
    cout << "\n";
    cout << "GRID_ELEMENT - Fatal error!\n";
    cout << "  Illegal value of CODE = \"" << code << "\".\n";
    exit ( 1 );
  }

  return element_node;
}
//****************************************************************************80

int grid_element_num ( string code, int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_ELEMENT_NUM returns the number of elements in a grid.
//
//  Discussion:
//
//    The number of elements generated will be NELEMX * NELEMY for
//    quadrilaterals, or 2 * NELEMX * NELEMY for triangles.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
//    'T4', 'T6' and 'T10'.
//
//    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
//    X and Y directions.  
//
//    Output, int GRID_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  if ( code == "Q4" )
  {
    element_num = grid_q4_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q8" )
  {
    element_num = grid_q8_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q9" )
  {
    element_num = grid_q9_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q12" )
  {
    element_num = grid_q12_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q16" )
  {
    element_num = grid_q16_element_num ( nelemx, nelemy );
  }
  else if ( code == "QL" )
  {
    element_num = grid_ql_element_num ( nelemx, nelemy );
  }
  else if ( code == "T3" )
  {
    element_num = grid_t3_element_num ( nelemx, nelemy );
  }
  else if ( code == "T4" )
  {
    element_num = grid_t4_element_num ( nelemx, nelemy );
  }
  else if ( code == "T6" )
  {
    element_num = grid_t6_element_num ( nelemx, nelemy );
  }
  else if ( code == "T10" )
  {
    element_num = grid_t10_element_num ( nelemx, nelemy );
  }
  else
  {
    cout << "\n";
    cout << "GRID_ELEMENT_NUM - Fatal error!\n";
    cout << "  Illegal value of CODE = \"" << code << "\".\n";
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
//****************************************************************************80

int grid_node_num ( string code, int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_NODE_NUM returns the number of nodes in a grid.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
//    'T4', 'T6' and 'T10'.
//
//    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
//    X and Y directions.  
//
//    Output, int GRID_NODE_NUM, the number of elements in the grid.
//
{
  int node_num;

  if ( code == "Q4" )
  {
    node_num = grid_q4_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q8" )
  {
    node_num = grid_q8_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q9" )
  {
    node_num = grid_q9_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q12" )
  {
    node_num = grid_q12_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q16" )
  {
    node_num = grid_q16_node_num ( nelemx, nelemy );
  }
  else if ( code == "QL" )
  {
    node_num = grid_ql_node_num ( nelemx, nelemy );
  }
  else if ( code == "T3" )
  {
    node_num = grid_t3_node_num ( nelemx, nelemy );
  }
  else if ( code == "T4" )
  {
    node_num = grid_t4_node_num ( nelemx, nelemy );
  }
  else if ( code == "T6" )
  {
    node_num = grid_t6_node_num ( nelemx, nelemy );
  }
  else if ( code == "T10" )
  {
    node_num = grid_t10_node_num ( nelemx, nelemy );
  }
  else
  {
    cout << "\n";
    cout << "GRID_NODE_NUM - Fatal error!\n";
    cout << "  Illegal value of CODE = \"" << code << "\".\n";
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
//****************************************************************************80

double *grid_nodes_01 ( int x_num, int y_num )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_NODES_01 returns an equally spaced rectangular grid of nodes in the unit square.
//
//  Example:
//
//    X_NUM = 5
//    Y_NUM = 3
//
//    NODE_XY =
//    ( 0, 0.25, 0.5, 0.75, 1, 0,   0.25, 0.5, 0.75, 1,   0, 0.25, 0.5, 0.75, 1;
//      0, 0,    0,   0,    0, 0.5, 0.5,  0.5, 0.5,  0.5, 1, 1.0,  1.0, 1.0,  1 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, Y_NUM, the number of nodes in the X and Y directions.
//
//    Output, double GRID_NODES_01[2*X_NUM*Y_NUM], the coordinates of
//    the nodes.
//
{
  int i;
  int j;
  int k;
  int node_num;
  double *node_xy;
  double value;

  node_num = x_num * y_num;

  node_xy = new double[2*node_num];

  if ( x_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[0+k*2] = 0.5;
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      value = ( double ) ( i ) / ( double ) ( x_num - 1 );
      for ( j = i; j < node_num; j = j + x_num )
      {
        node_xy[0+j*2] = value;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( k = 0; k < node_num; k++ )
    {
      node_xy[1+k*2] = 0.5;
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      value = ( double ) ( j ) / ( double ) ( y_num - 1 );
      for ( i = j*x_num; i < ( j + 1 ) * x_num; i++ )
      {
        node_xy[1+i*2] = value;
      }
    }
  }

  return node_xy;
}
//****************************************************************************80

void grid_print ( int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_PRINT prints the elements that form a grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the number of nodes per element.
//
//    Input, int ELEMENT_NUM, the number of elements.
// 
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
//    each element.
//
{
  int element;
  int i;

  cout << "\n";
  cout << "  GRID_PRINT: Element -> Node table.\n";
  cout << "\n";
  cout << "  Element order = " << element_order << "\n";
  cout << "  Number of elements = " << element_num << "\n";
  cout << "\n";
  cout << "    #   ";
  for ( i = 0; i < element_order; i++ )
  {
    cout << setw(3) << i;
  }
  cout << "\n";
  cout << "\n";

  for ( element = 0; element < element_num; element++ )
  {
    cout << "  " << setw(3) << element << "   ";
    for ( i = 0; i < element_order; i++ )
    {
      cout <<  setw(3) << element_node[i+element*element_order];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

int *grid_q4_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_Q4_ELEMENT produces a grid of 4 node quadrilaterals.
//
//  Discussion:
//
//    For each element, the nodes are listed in counter-clockwise order.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1, 2,  6,  5;
//         2, 3,  7,  6;
//         3, 4,  8,  7;
//         5, 6, 10,  9;
//         6, 7, 11, 10;
//         7, 8, 12, 11.
//
//  Grid:
//
//    9---10---11---12
//    |    |    |    |
//    |    |    |    |
//    |  4 |  5 |  6 |
//    |    |    |    |
//    5----6----7----8
//    |    |    |    |
//    |    |    |    |
//    |  1 |  2 |  3 |
//    |    |    |    |
//    1----2----3----4
//
//  Element Q4:
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
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q4[4*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = new int[element_order*nelemx*nelemy];
//
//  Node labeling:
//
//    NW---NE
//     |    |
//    SW---SE
//
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );
  
      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_q4_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q4_ELEMENT_NUM counts the elements in a grid of 4 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_q4_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q4_NODE_NUM counts the nodes in a grid of 4 node quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_Q4_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
//****************************************************************************80

int *grid_q8_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_Q8_ELEMENT produces a grid of 8 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  3, 14, 12,  2,  9, 13,  8;
//         3,  5, 16, 14,  4, 10, 15,  9;
//         5,  7, 18, 16,  6, 11, 17, 10;
//        12, 14, 25, 23, 13, 20, 24, 19;
//        14, 16, 27, 25, 15, 21, 26, 20;
//        16, 18, 29, 27, 17, 22, 28, 21.
//
//  Diagram:
//
//   23---24---25---26---27---28---29
//    |         |         |         |
//    |         |         |         |
//   19        20        21        22
//    |         |         |         |
//    | 4       | 5       | 6       |
//   12---13---14---15---16---17---18
//    |         |         |         |
//    |         |         |         |
//    8         9        10        11
//    |         |         |         |
//    | 1       | 2       | 3       |
//    1----2----3----4----5----6----7
//
//  Element Q8:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8     6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q8[8*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int e;
  int element;
  int *element_node;
  int element_order = 8;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = new int[element_order*nelemx*nelemy];
//
//  Node labeling:
//
//    NW----N----NE
//     |          |
//     W   (C)    E
//     |          |
//    SW----S----SE
//

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = ( j - 1 )  * ( 3 * nelemx + 2 ) + 2 * i - 1;
      w  = sw + 2 * nelemx + 2 - i;
      nw = sw + 3 * nelemx + 2;

      s =  sw + 1;
      n =  sw + ( 3 * nelemx + 2 ) + 1;

      se = sw + 2;
      e  = sw + 2 * nelemx + 2 - i + 1;
      ne = sw + ( 3 * nelemx + 2 ) + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_q8_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q8_ELEMENT_NUM counts the elements in a grid of 8 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q8_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_q8_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q8_NODE_NUM counts the nodes in a grid of 8 node quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_Q8_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 3 * nelemx * nelemy + 2 * nelemx + 2 * nelemy + 1;

  return node_num;
}
//****************************************************************************80

int *grid_q9_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_Q9_ELEMENT produces a grid of 9 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  3, 17, 15,  2, 10, 16,  8,  9;
//         3,  5, 19, 17,  4, 12, 18, 10, 11;
//         5,  7, 21, 19,  6, 14, 20, 12, 13;
//        15, 17, 31, 29, 16, 24, 30, 22, 23;
//        17, 19, 33, 31, 18, 26, 32, 24, 25;
//        19, 21, 35, 33, 20, 28, 34, 26, 27.
//
//  Grid:
//
//   29---30---31---32---33---34---35
//    |    .    |    .    |    .    |
//    |    .    |    .    |    .    |
//   22 . 23 . 24 . 25 . 26 . 27 . 28
//    |    .    |    .    |    .    |
//    | 4  .    | 5  .    | 6  .    |
//   15---16---17---18---19---20---21
//    |    .    |    .    |    .    |
//    |    .    |    .    |    .    |
//    8 .  9 . 10 . 11 . 12 . 13 . 14
//    |    .    |    .    |    .    |
//    | 1  .    | 2  .    | 3  .    |
//    1----2----3----4----5----6----7
//
//  Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q9[9*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 9;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = new int[element_order*nelemx*nelemy];
//
//  Node labeling:
//
//    NW----N----NE
//     |          |
//     W    C     E
//     |          |
//    SW----S----SE
//
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = ne;
      element_node[3+element*element_order] = nw;
      element_node[4+element*element_order] = s;
      element_node[5+element*element_order] = e;
      element_node[6+element*element_order] = n;
      element_node[7+element*element_order] = w;
      element_node[8+element*element_order] = c;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_q9_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q9_ELEMENT_NUM counts the elements in a grid of 9 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_q9_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q9_NODE_NUM counts the nodes in a grid of 9 node quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_Q9_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
//****************************************************************************80

int *grid_q12_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_Q12_ELEMENT produces a grid of 12 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  2,  3,  4, 11, 12, 15, 16, 19, 20, 21, 22;
//         4,  5,  6,  7, 12, 13, 16, 17, 22, 23, 24, 25;
//         7,  8,  9, 10, 13, 14, 17, 18, 25, 26, 27, 28;
//        19, 20, 21, 22, 29, 30, 33, 34, 37, 38, 39, 40;
//        22, 23, 24, 25, 30, 31, 34, 35, 40, 41, 42, 43;
//        25, 26, 27, 28, 31, 32, 35, 36, 43, 44, 45, 46.
//
//  Grid:
//
//   37-38-39-40-41-42-43-44-45-46
//    |        |        |        |
//   33       34       35       36
//    |        |        |        |
//   29       30       31       32
//    | 4      | 5      | 6      |
//   19-20-21-22-23-24-25-26-27-28
//    |        |        |        |
//   15       16       17       18
//    |        |        |        |
//   11       12       13       14
//    | 1      | 2      | 3      |
//    1--2--3--4--5--6--7--8--9-10
//
//  Element Q12:
//
//    |
//    1  9-10-11-12
//    |  |        |
//    |  7        8
//    S  |        |
//    |  5        6
//    |  |        |
//    0  1--2--3--4
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q12[12*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int base;
  int element;
  int *element_node;
  int element_order = 12;
  int i;
  int j;

  element_node = new int[element_order*nelemx*nelemy];

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 5 * nelemx + 3 ) + 1;

      element_node[ 0+element*element_order] =  base + ( i - 1 ) * 3;
      element_node[ 1+element*element_order] =  base + ( i - 1 ) * 3 + 1;
      element_node[ 2+element*element_order] =  base + ( i - 1 ) * 3 + 2;
      element_node[ 3+element*element_order] =  base + ( i - 1 ) * 3 + 3;

      element_node[ 4+element*element_order] =  base + 3 * nelemx + i;
      element_node[ 5+element*element_order] =  base + 3 * nelemx + i + 1;

      element_node[ 6+element*element_order] =  base + 4 * nelemx + i + 1;
      element_node[ 7+element*element_order] =  base + 4 * nelemx + i + 2;

      element_node[ 8+element*element_order] =  base + 5 * nelemx + 3 * i;
      element_node[ 9+element*element_order] = base + 5 * nelemx + 3 * i + 1;
      element_node[10+element*element_order] = base + 5 * nelemx + 3 * i + 2;
      element_node[11+element*element_order] = base + 5 * nelemx + 3 * i + 3;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_q12_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q12_ELEMENT_NUM counts the elements in a grid of 12 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q12_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_q12_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q12_NODE_NUM counts the nodes in a grid of 12 node quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_Q12_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 5 * nelemx * nelemy + 3 * nelemx + 3 * nelemy + 1;

  return node_num;
}
//****************************************************************************80

int *grid_q16_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_Q16_ELEMENT produces a grid of 16 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 2, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  2,  3,  4,  8,  9, 10, 11, 15, 16, 17, 18, 22, 23, 24, 25;
//         4,  5,  6,  7, 11, 12, 13, 14, 18, 19, 20, 21, 25, 26, 27, 28;
//        22, 23, 24, 25, 29, 30, 31, 32, 36, 37, 38, 39, 43, 44, 45, 46;
//        25, 26, 27, 28, 32, 33, 34, 35, 39, 40, 41, 42, 46, 47, 48, 49. 
//        
//  Grid:
//
//   43-44-45-46-47-48-49
//    |        |        |
//    |        |        |
//   36 37 38 39 40 41 42
//    |        |        |
//    |        |        |
//   29 30 31 32 33 34 35
//    |        |        |
//    | 3      | 4      |
//   22-23-24-25-26-27-28
//    |        |        |
//    |        |        |
//   15 16 17 18 19 20 21
//    |        |        |
//    |        |        |
//    8  9 10 11 12 13 14
//    |        |        |
//    | 1      | 2      |
//    1--2--3--4--5--6--7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q16[16*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int base;
  int element;
  int *element_node;
  int element_order = 16;
  int i;
  int j;

  element_node = new int[element_order*nelemx*nelemy];

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[ 0+element*element_order] = base;
      element_node[ 1+element*element_order] = base                          + 1;
      element_node[ 2+element*element_order] = base                          + 2;
      element_node[ 3+element*element_order] = base                          + 3;
      element_node[ 4+element*element_order] = base +     ( 3 * nelemx + 1 );
      element_node[ 5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 1;
      element_node[ 6+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[ 7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[ 8+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[ 9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[10+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element_node[11+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[12+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[13+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[14+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[15+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_q16_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q16_ELEMENT_NUM counts the elements in a grid of 16 node quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_q16_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_Q16_NODE_NUM counts the nodes in a grid of 16 node quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_Q16_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
//****************************************************************************80

int *grid_ql_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_QL_ELEMENT produces a grid of 6 node quadratics/linears.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  2,  3,  8,  9, 10;
//         3,  4,  5, 10, 11, 12;
//         5,  6,  7, 12, 13, 14;
//         8,  9, 10, 15, 16, 17;
//        10, 11, 12, 17, 18, 19;
//        12, 13, 14, 19, 20, 21.
//
//  Grid:
//
//   15---16---17---18---19---20---21
//    |         |         |         |
//    |         |         |         |
//    |    4    |    5    |    6    |
//    |         |         |         |
//    |         |         |         |
//    8----9---10---11---12---13---14
//    |         |         |         |
//    |         |         |         |
//    |    1    |    2    |    3    |
//    |         |         |         |
//    |         |         |         |
//    1----2----3----4----5----6----7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.  X will the the "quadratic direction", and
//    Y will be the "linear direction".
//
//    Output, int GRID_QL[6*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int base;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;

  element_node = new int[element_order*nelemx*nelemy];

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * i - 1;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base + 1;
      element_node[2+element*element_order] = base + 2;
      element_node[3+element*element_order] = base + ( 2 * nelemx + 1 );
      element_node[4+element*element_order] = base + ( 2 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base + ( 2 * nelemx + 1 ) + 2;

      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_ql_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_QL_ELEMENT_NUM counts the elements in a grid of quadratic/linear quadrilaterals.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int GRID_QL_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_ql_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_QL_NODE_NUM counts the nodes in a grid of quadratic/linear quadrilaterals.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_QL_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 2 * nelemx * nelemy + 2 * nelemx + nelemy + 1;

  return node_num;
}
//****************************************************************************80

void grid_shape_2d ( int n, double a[], int *n1, int *n2 )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_SHAPE_2D guesses the shape N1 by N2 of a vector of data.
//
//  Discussion:
//
//    The data vector A is assumed to contain N1 * N2 values, with
//    where each of N2 values is repeated N1 times.
//
//  Example:
//
//    Input:
//
//      A = ( 2, 2, 2, 7, 7, 7 )
//
//    Output:
//
//      N1 = 3, N2 = 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data values.
//
//    Input, double A[N], the data, which should have the properties
//    described above.
//
//    Output, int *N1, *N2, the "shape" of the data in the array.
//
{
  int i;
//
//  Make a guess for N1.
//
  i = 1;
  *n1 = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[0] )
    {
      break;
    }
    *n1 = *n1 + 1;
  }
//
//  Guess that N2 = N / N1.
//
  *n2 = n / (*n1);

  return;
}
//****************************************************************************80

int *grid_t3_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_T3_ELEMENT produces a grid of pairs of 3 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  2,  5;
//         6,  5,  2;
//         2,  3,  6;
//         7,  6,  3;
//         3,  4,  7;
//         8,  7,  4;
//         5,  6,  9;
//        10,  9,  6;
//         6,  7, 10;
//        11, 10,  7;
//         7,  8, 11;
//        12, 11,  8.
//
//  Grid:
//
//    9---10---11---12
//    |\ 8 |\10 |\12 |
//    | \  | \  | \  |
//    |  \ |  \ |  \ |
//    |  7\|  9\| 11\|
//    5----6----7----8
//    |\ 2 |\ 4 |\ 6 |
//    | \  | \  | \  |
//    |  \ |  \ |  \ |
//    |  1\|  3\|  5\|
//    1----2----3----4
//
//  Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T3[3*2*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int element;
  int *element_node;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_node = new int[element_order*2*nelemx*nelemy];
//
//  Node labeling:
//
//    NW--NE
//     |\ |
//     | \|
//    SW--SE
//
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( nelemx + 1 );
      se = i + 1 + ( j - 1 ) * ( nelemx + 1 );
      nw = i     +   j       * ( nelemx + 1 );
      ne = i + 1 +   j       * ( nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_t3_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T3_ELEMENT_NUM counts the elements in a grid of 3 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T3_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_t3_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T3_NODE_NUM counts the nodes in a grid of 3 node triangles.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_T3_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 );

  return node_num;
}
//****************************************************************************80

int *grid_t4_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_T4_ELEMENT produces a grid of pairs of 4 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  2,  11,  5;
//        12, 11,   2,  8;
//         2,  3,  12,  6;
//        13, 12,   3,  9;
//         3   4   13,  7;
//        14, 13,   4,  10;
//        11, 12,  21,  15;
//        22, 21,  12,  18;
//        12, 13,  22,  16;
//        23, 22,  13,  19;
//        13  14   23,  17;
//        24, 23,  14,  20;
//
//  Grid:
//
//   21---22---23---24
//    |\18 |\19 |\20 |
//    | \  | \  | \  |
//    |  \ |  \ |  \ |
//    | 15\| 16\| 17\|
//   11---12---13---14
//    |\ 8 |\ 9 |\10 |
//    | \  | \  | \  |
//    |  \ |  \ |  \ |
//    | 5 \|  6\|  7\|
//    1----2----3----4
//
//  Element T4:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  | 4 \
//    |  |    \
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
//    23 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T4[4*2*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int element;
  int *element_node;
  int element_order = 4;
  int i;
  int j;
  int nc;
  int ne;
  int nw;
  int sc;
  int se;
  int sw;

  element_node = new int[element_order*2*nelemx*nelemy];
//
//  Node labeling:
//
//    NW----NE
//     |\   |
//     | \NC|
//     |SC\ |
//     |   \|
//    SW---SE
//
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = i     + ( j - 1 ) * ( 3 * nelemx + 1 );
      se = sw + 1;
      sc = sw +     nelemx + 1;
      nc = sw + 2 * nelemx + 1;
      nw = sw + 3 * nelemx + 1;
      ne = sw + 3 * nelemx + 2;

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = sc;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = nc;
      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_t4_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T4_ELEMENT_NUM counts the elements in a grid of 4 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T4_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_t4_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T4_NODE_NUM counts the nodes in a grid of 4 node triangles.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_T4_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( nelemx + 1 ) * ( nelemy + 1 ) + 2 * nelemx * nelemy;

  return node_num;
}
//****************************************************************************80

int *grid_t6_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_T6_ELEMENT produces a grid of pairs of 6 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  3, 15,  2,  9,  8;
//        17, 15,  3, 16,  9, 10;
//         3,  5, 17,  4, 11, 10;
//        19, 17,  5, 18, 11, 12;
//         5,  7, 19,  6, 13, 12;
//        21, 19,  7, 20, 13, 14;
//        15, 17, 29, 16, 23, 22;
//        31, 29, 17, 30, 23, 24;
//        17, 19, 31, 18, 25, 24;
//        33, 31, 19, 32, 25, 26;
//        19, 21, 33, 20, 27, 26;
//        35, 33, 21, 34, 27, 28.
//
//  Grid:
//
//   29-30-31-32-33-34-35
//    |\ 8  |\10  |\12  |
//    | \   | \   | \   |
//   22 23 24 25 26 27 28
//    |   \ |   \ |   \ |
//    |  7 \|  9 \| 11 \|
//   15-16-17-18-19-20-21
//    |\ 2  |\ 4  |\ 6  |
//    | \   | \   | \   |
//    8  9 10 11 12 13 14
//    |   \ |   \ |   \ |
//    |  1 \|  3 \|  5 \|
//    1--2--3--4--5--6--7
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T6[6*2*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int c;
  int e;
  int element;
  int *element_node;
  int element_order = 6;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_node = new int[element_order*2*nelemx*nelemy];
//
//  Node labeling:
//
//    NW---N--NE
//     | \     |
//     W   C   E
//     |    \  |
//    SW---S--SE
//
  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      sw = 2 * ( j - 1 )  * ( 2 * nelemx + 1 ) + 2 * ( i - 1 ) + 1;
      w  = sw +               2 * nelemx + 1;
      nw = sw +         2 * ( 2 * nelemx + 1 );

      s  = sw + 1;
      c  = sw + 1 +               2 * nelemx + 1;
      n  = sw + 1 +         2 * ( 2 * nelemx + 1 );

      se = sw + 2;
      e  = sw + 2 +               2 * nelemx + 1;
      ne = sw + 2 +         2 * ( 2 * nelemx + 1 );

      element_node[0+element*element_order] = sw;
      element_node[1+element*element_order] = se;
      element_node[2+element*element_order] = nw;
      element_node[3+element*element_order] = s;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = w;
      element = element + 1;

      element_node[0+element*element_order] = ne;
      element_node[1+element*element_order] = nw;
      element_node[2+element*element_order] = se;
      element_node[3+element*element_order] = n;
      element_node[4+element*element_order] = c;
      element_node[5+element*element_order] = e;
      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_t6_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T6_ELEMENT_NUM counts the elements in a grid of 6 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T6_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_t6_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T6_NODE_NUM counts the nodes in a grid of 6 node triangles.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_T6_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( 2 * nelemx + 1 ) * ( 2 * nelemy + 1 );

  return node_num;
}
//****************************************************************************80

int *grid_t10_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_T10_ELEMENT produces a grid of pairs of 10 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 2, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE = 
//         1,  2,  3,  4, 10, 16, 22, 15,  8,  9;
//        25, 24, 23, 22, 16, 10,  4, 11, 18, 17;
//         4,  5,  6,  7, 13, 19, 25, 18, 11, 12;
//        28, 27, 26, 25, 19, 13,  7, 14, 21, 20;
//        22, 23, 24, 25, 31, 37, 43, 36, 29, 30;
//        46, 45, 44, 43, 37, 31, 25, 32, 39, 38;
//        25, 26, 27, 28, 34, 40, 46, 39, 31, 33;
//        49, 48, 47, 46, 40, 34, 28, 35, 42, 41.
//        
//  Grid:
//
//   43-44-45-46-47-48-49
//    |\     6 |\     8 |
//    | \      | \      |
//   36 37 38 39 40 41 42
//    |   \    |   \    |
//    |    \   |    \   |
//   29 30 31 32 33 34 35
//    |      \ |      \ |
//    | 5     \| 7     \|
//   22-23-24-25-26-27-28
//    |\     2 |\     4 |
//    | \      | \      |
//   15 16 17 18 19 20 21
//    |   \    |   \    |
//    |    \   |    \   |
//    8  9 10 11 12 13 14
//    |      \ |      \ |
//    | 1     \| 3     \|
//    1--2--3--4--5--6--7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T10[10*2*NELEMX*NELEMY], the nodes that form
//    each element.
//
{
  int base;
  int element;
  int *element_node;
  int element_order = 10;
  int i;
  int j;

  element_node = new int[element_order*2*nelemx*nelemy];

  element = 0;
 
  for ( j = 1; j <= nelemy; j++ )
  {
    for ( i = 1; i <= nelemx; i++ )
    {
      base = ( j - 1 ) * 3 * ( 3 * nelemx + 1 ) + 3 * i - 2;

      element_node[0+element*element_order] = base;
      element_node[1+element*element_order] = base                          + 1;
      element_node[2+element*element_order] = base                          + 2;
      element_node[3+element*element_order] = base                          + 3;
      element_node[4+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[5+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[6+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[7+element*element_order] = base + 2 * ( 3 * nelemx + 1 );
      element_node[8+element*element_order] = base +     ( 2 * nelemx + 1 ) + 2;
      element_node[9+element*element_order] = base +     ( 2 * nelemx + 1 ) + 3;
      element = element + 1;

      element_node[0+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 3;
      element_node[1+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 2;
      element_node[2+element*element_order] = base + 3 * ( 3 * nelemx + 1 ) + 1;
      element_node[3+element*element_order] = base + 3 * ( 3 * nelemx + 1 );
      element_node[4+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 1;
      element_node[5+element*element_order] = base +     ( 3 * nelemx + 1 ) + 2;
      element_node[6+element*element_order] = base                          + 3;
      element_node[7+element*element_order] = base +     ( 3 * nelemx + 1 ) + 3;
      element_node[8+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 3;
      element_node[9+element*element_order] = base + 2 * ( 3 * nelemx + 1 ) + 2;
      element = element + 1;
    }
  }

  return element_node;
}
//****************************************************************************80

int grid_t10_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T10_ELEMENT_NUM counts the elements in a grid of 10 node triangles.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * NELEMY = 12
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    2 * NELEMX * NELEMY.
//
//    Output, int GRID_T10_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int grid_t10_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_T10_NODE_NUM counts the nodes in a grid of 10 node triangles.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, int GRID_T10_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = ( 3 * nelemx + 1 ) * ( 3 * nelemy + 1 );

  return node_num;
}
//****************************************************************************80

void grid_test ( string code )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_TEST tests the grid routines.
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
//  Parameters:
//
//    Input, string CODE, the code for the element.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T4', 'T6' and 'T10'.
//
{
  int *element_node;
  int element_num;
  int element_order;
  int nelemx;
  int nelemy;
  int width;
//
//  NODE is defined as a vector rather than a two dimensional array,
//  so that we can handle the various cases using a single array.
//
  cout << "\n";
  cout << "  GRID_TEST: Test the grid routine for element " << code << "\n";

  nelemx = 3;
  nelemy = 2;

  if ( code == "Q4" || 
       code == "Q8" || 
       code == "Q9" ||
       code == "Q12" ||
       code == "Q16" ||
       code == "QL" )
  {
    element_num = nelemx * nelemy;
  }
  else if ( code == "T3" || code == "T4" || code == "T6" || code == "T10" )
  {
    element_num = 2 * nelemx * nelemy;
  }

  element_order = order_code ( code );

  element_node = grid_element ( code, element_order, nelemx, nelemy );

  grid_print ( element_order, element_num, element_node );

  width = grid_width ( element_order, element_num, element_node );

  cout << "\n";
  cout << "  Grid width is " << width << "\n";

  delete [] element_node;

  return;
}
//****************************************************************************80

int grid_width ( int element_order, int element_num, int element_node[] )

//****************************************************************************80
//
//  Purpose: 
//
//    GRID_WIDTH computes the width of a given grid.
//
//  Definition:
//
//    The grid width is defined to the maximum absolute
//    difference of global indices of nodes in the same element.
//
//  Example:
//
//    For the following grid, the grid width is 13.
//
//   23---24---25---26---27---28---29
//    |         |         |         |
//    |         |         |         |
//   19        20        21        22
//    |         |         |         |
//    | 4       | 5       | 6       |
//   12---13---14---15---16---17---18
//    |         |         |         |
//    |         |         |         |
//    8         9        10        11
//    |         |         |         |
//    | 1       | 2       | 3       |
//    1----2----3----4----5----6----7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
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
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], the nodes that form
//    each element.
//
//    Output, int GRID_WIDTH, the grid width.
//
{
  int element;
  int ip1;
  int ip2;
  int node1;
  int node2;
  int width;

  width = 0;
 
  for ( element = 0; element < element_num; element++ )
  {
    for ( node1 = 0; node1 < element_order; node1++ )
    {
      ip1 = element_node[node1+element*element_order];
      for ( node2 = 0; node2 < element_order; node2++ )
      {
        ip2 = element_node[node2+element*element_order];
        width = i4_max ( width, abs ( ip1 - ip2 ) );
      }
    }
  }
 
  return width;
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
//    Input, int I1, I2, are two I4's to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
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
//    Input, int I1, I2, two I4's to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

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
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
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
//    I  I4_WRAP
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
//    Input, int IVAL, an I4 value.
//
//    Input, int ILO, IHI, the desired bounds for the value.
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
//    Input, string TITLE, a title to be printed.
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

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
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
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ ) 
  {
    cout << setw(6) << i + 1 << "  " 
         << setw(8) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

void interp ( string code, int element_order, double r, double s, 
  double ubase[], double *u, double *dudr, double *duds )

//****************************************************************************80
//
//  Purpose: 
//
//    INTERP interpolates a quantity in an element from basis node values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T6' and 'T10'.
//
//    Input, int ELEMENT_ORDER, order of the element.
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Input, double UBASE[ELEMENT_ORDER], the value of the quantity 
//    at the basis nodes.
//
//    Output, double *U, *DUDR, *DUDS, the interpolated value of the
//    quantity and its derivatives at the point (R,S).
//
{
  double *dtdr;
  double *dtds;
  int i;
  double *t;

  dtdr = new double[element_order];
  dtds = new double[element_order];
  t = new double[element_order];

  shape ( code, r, s, t, dtdr, dtds );
 
  *u = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *u = *u + ubase[i] * t[i];
  }

  *dudr = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *dudr = *dudr + ubase[i] * dtdr[i];
  }

  *duds = 0.0;
  for ( i = 0; i < element_order; i++ )
  {
    *duds = *duds + ubase[i] * dtds[i];
  }

  delete [] dtdr;
  delete [] dtds;
  delete [] t;
 
  return;
}
//****************************************************************************80

void interp_test ( string code )

//****************************************************************************80
//
//  Purpose:
//
//    INTERP_TEST tests the interpolation property of an element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element.
//    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", 
//    "T3", "T4", "T6" and "T10".
//
{
  double area;
  double dudr;
  double dudr_exact;
  double duds;
  double duds_exact;
  int element_order;
  int i;
  int node;
  double *node_r;
  double *node_s;
  double *node_u;
  double r;
  double r_factor;
  int *rexp;
  double s;
  double s_factor;
  int *sexp;
  int seed;
  int test;
  int test_num = 5;
  double u;
  double u_exact;

  if ( code == "T4" )
  {
    cout << "\n";
    cout << "INTERP_TEST - Warning!\n";
    cout << "  Skipping test for element \"T4\".\n";
    return;
  }

  cout << "\n";
  cout << "INTERP_TEST for element \"" << code << "\".\n";

  element_order = order_code ( code );

  cout << "  Number of nodes = " << element_order << "\n";

  node_r = new double [element_order];
  node_s = new double [element_order];
  node_u = new double [element_order];
  rexp = new int [element_order];
  sexp = new int[element_order];
//
//  Get the coordinates of the reference nodes.
//
  node_reference ( code, node_r, node_s, &area );
//
//  Get the monomial exponents for which the element is exact.
//
  poly ( code, rexp, sexp );

  seed = 123456789;

  for ( i = 0; i < element_order; i++ )
  {
    cout << "\n";
    cout << "  Interpolate R^" << rexp[i] << " * S^" << sexp[i] << "\n";
    cout << "\n";
//
//  Evaluate R**REXP(I) * S**SEXP(I) at the nodes.  This is our data.
//
    for ( node = 0; node < element_order; node++ )
    {
      r = node_r[node];
      s = node_s[node];
      if ( rexp[i] == 0 )
      {
        r_factor = 1.0;
      }
      else
      {
        r_factor = pow ( r, rexp[i] );
      }
      if ( sexp[i] == 0 )
      {
        s_factor = 1.0;
      }
      else
      {
        s_factor = pow ( s, sexp[i] );
      }
      node_u[node] = r_factor * s_factor;
      cout << "  (R,S,U):     " 
           << "  " << setw(12) << r
           << "  " << setw(12) << s
           << "  " << setw(12) << node_u[node] << "\n";
    }
//
//  Now pick random points in the element, and compute the interpolated
//  value of R**REXP(*) * S**SEXP(I) there.  Mathematically, these
//  values should be exact.
//
    for ( test = 1; test <= test_num; test++ )
    {
      reference_sample ( code, &seed, &r, &s );

      cout << "\n";
      cout << "  (R,S):"
           << "  " << setw(12) << r
           << "  " << setw(12) << s << "\n";

      u_exact = r8_power ( r, rexp[i] ) * r8_power ( s, sexp[i] );

      dudr_exact = ( double ) ( rexp[i] ) 
        * r8_power ( r, rexp[i] - 1 ) * r8_power ( s, sexp[i] );

      duds_exact = r8_power ( r, rexp[i] ) * ( double ) ( sexp[i] ) 
        * r8_power ( s, sexp[i] - 1 );

      interp ( code, element_order, r, s, node_u, &u, &dudr, &duds );

      cout << "  (U ,U* ,Error): "
           << "  " << setw(12) << u_exact
           << "  " << setw(12) << u
           << "  " << setw(12) << fabs ( u_exact - u ) << "\n";

      cout << "  (Ur,Ur*,Error): "
           << "  " << setw(12) << dudr_exact
           << "  " << setw(12) << dudr
           << "  " << setw(12) << fabs ( dudr_exact - dudr ) << "\n";

      cout << "  (Us,Us*,Error): "
           << "  " << setw(12) << duds_exact
           << "  " << setw(12) << duds
           << "  " << setw(12) << fabs ( duds_exact - duds ) << "\n";
    }
  }

  delete [] node_r;
  delete [] node_s;
  delete [] node_u;
  delete [] rexp;
  delete [] sexp;

  return;
}
//****************************************************************************80

void legendre_com ( int norder, double xtab[], double weight[] )

//****************************************************************************80
//
//  Purpose: 
//
//    LEGENDRE_COM computes abscissas and weights for Gauss-Legendre quadrature.
//
//  Integration interval:
//
//    [ -1, 1 ]
//
//  Weight function:
//
//    1.
//
//  Integral to approximate:
//
//    Integral ( -1 <= X <= 1 ) F(X) dX.
//
//  Approximate integral:
//
//    sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NORDER, the order of the rule.
//    NORDER must be greater than 0.
//
//    Output, double XTAB[NORDER], the abscissas of the rule.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//    The weights are positive, symmetric, and should sum to 2.
//
{
# define PI 3.141592653589793

  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( norder < 1 )
  {
    cout << "\n";
    cout << "LEGENDRE_COM - Fatal error!\n";
    cout << "  Illegal value of NORDER = " << norder << "\n";
    exit ( 1 );
  }
 
  e1 = ( double ) ( norder * ( norder + 1 ) );
 
  m = ( norder + 1 ) / 2;
 
  for ( i = 1; i <= ( norder + 1 ) / 2; i++ )
  {
    mp1mi = m + 1 - i;
    t = PI * ( double ) ( 4 * i - 1 ) / ( double ) ( 4 * norder + 2 );
    x0 = cos(t) * ( 1.0 - ( 1.0 - 1.0 / 
      ( double ) ( norder ) ) / ( double ) ( 8 * norder * norder ) );
 
    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= norder; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }
 
    d1 = ( double ) ( norder ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
//
//  Initial approximation H:
//
    h = - u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn 
      / ( 3.0 * dpn ) ) ) );
//
//  Refine H using one step of Newton's method:
//
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0 
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;
 
    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;
 
    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0 
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx ); 
  }
 
  if ( ( norder % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
//
//  Shift the data up.
//
  nmove = ( norder + 1 ) / 2;
  ncopy = norder - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = norder + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
//
//  Reflect values for the negative abscissas.
//
  for ( i = 0; i < norder - nmove; i++ )
  {
    xtab[i] = - xtab[norder-1-i];
    weight[i] = weight[norder-1-i];
  }
 
  return;

# undef PI
}
//****************************************************************************80

void legendre_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
//
//  Discussion:
//
//    The integration interval is [ -1, 1 ].
//
//    The weight function w(x-1] = 1.0;
//
//    The integral to approximate:
//
//      Integral ( -1 <= X <= 1 ) F(X) dX
//
//    Quadrature rule:
//
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//
//    The quadrature rule will integrate exactly all polynomials up to
//    X**(2*N-1).
//
//    The abscissas of the rule are the zeroes of the Legendre polynomial
//    P(N)(X).
//
//    The integral produced by a Gauss-Legendre rule is equal to the
//    integral of the unique polynomial of degree N-1 which
//    agrees with the function at the N abscissas of the rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 October 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798.
//
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//    N must be between 1 and 33, 63, 64, 65, 127 or 255.
//
//    Output, double X[N], the abscissas of the rule.
//
//    Output, double W[N], the weights of the rule.
//    The weights are positive, symmetric and should sum to 2.
//
{
  if ( n == 1 )
  {
    x[0] =   0.0;

    w[0] = 2.0;
  }
  else if ( n == 2 )
  {
    x[0] = - 0.577350269189625764509148780502;
    x[1] =   0.577350269189625764509148780502;

    w[0] = 1.0;
    w[1] = 1.0;
  }
  else if ( n == 3 )
  {
    x[0] = - 0.774596669241483377035853079956;
    x[1] =   0.0;
    x[2] =   0.774596669241483377035853079956;

    w[0] = 5.0 / 9.0;
    w[1] = 8.0 / 9.0;
    w[2] = 5.0 / 9.0;
  }
  else if ( n == 4 )
  {
    x[0] = - 0.861136311594052575223946488893;
    x[1] = - 0.339981043584856264802665759103;
    x[2] =   0.339981043584856264802665759103;
    x[3] =   0.861136311594052575223946488893;

    w[0] = 0.347854845137453857373063949222;
    w[1] = 0.652145154862546142626936050778;
    w[2] = 0.652145154862546142626936050778;
    w[3] = 0.347854845137453857373063949222;
  }
  else if ( n == 5 )
  {
    x[0] = - 0.906179845938663992797626878299;
    x[1] = - 0.538469310105683091036314420700;
    x[2] =   0.0;
    x[3] =   0.538469310105683091036314420700;
    x[4] =   0.906179845938663992797626878299;

    w[0] = 0.236926885056189087514264040720;
    w[1] = 0.478628670499366468041291514836;
    w[2] = 0.568888888888888888888888888889;
    w[3] = 0.478628670499366468041291514836;
    w[4] = 0.236926885056189087514264040720;
  }
  else if ( n == 6 )
  {
    x[0] = - 0.932469514203152027812301554494;
    x[1] = - 0.661209386466264513661399595020;
    x[2] = - 0.238619186083196908630501721681;
    x[3] =   0.238619186083196908630501721681;
    x[4] =   0.661209386466264513661399595020;
    x[5] =   0.932469514203152027812301554494;

    w[0] = 0.171324492379170345040296142173;
    w[1] = 0.360761573048138607569833513838;
    w[2] = 0.467913934572691047389870343990;
    w[3] = 0.467913934572691047389870343990;
    w[4] = 0.360761573048138607569833513838;
    w[5] = 0.171324492379170345040296142173;
  }
  else if ( n == 7 )
  {
    x[0] = - 0.949107912342758524526189684048;
    x[1] = - 0.741531185599394439863864773281;
    x[2] = - 0.405845151377397166906606412077;
    x[3] =   0.0;
    x[4] =   0.405845151377397166906606412077;
    x[5] =   0.741531185599394439863864773281;
    x[6] =   0.949107912342758524526189684048;

    w[0] = 0.129484966168869693270611432679;
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489;
    w[3] = 0.417959183673469387755102040816;
    w[4] = 0.381830050505118944950369775489;
    w[5] = 0.279705391489276667901467771424;
    w[6] = 0.129484966168869693270611432679;
  }
  else if ( n == 8 )
  {
    x[0] = - 0.960289856497536231683560868569;
    x[1] = - 0.796666477413626739591553936476;
    x[2] = - 0.525532409916328985817739049189;
    x[3] = - 0.183434642495649804939476142360;
    x[4] =   0.183434642495649804939476142360;
    x[5] =   0.525532409916328985817739049189;
    x[6] =   0.796666477413626739591553936476;
    x[7] =   0.960289856497536231683560868569;

    w[0] = 0.101228536290376259152531354310;
    w[1] = 0.222381034453374470544355994426;
    w[2] = 0.313706645877887287337962201987;
    w[3] = 0.362683783378361982965150449277;
    w[4] = 0.362683783378361982965150449277;
    w[5] = 0.313706645877887287337962201987;
    w[6] = 0.222381034453374470544355994426;
    w[7] = 0.101228536290376259152531354310;
  }
  else if ( n == 9 )
  {
    x[0] = - 0.968160239507626089835576202904;
    x[1] = - 0.836031107326635794299429788070;
    x[2] = - 0.613371432700590397308702039341;
    x[3] = - 0.324253423403808929038538014643;
    x[4] =   0.0;
    x[5] =   0.324253423403808929038538014643;
    x[6] =   0.613371432700590397308702039341;
    x[7] =   0.836031107326635794299429788070;
    x[8] =   0.968160239507626089835576202904;

    w[0] = 0.812743883615744119718921581105E-01;
    w[1] = 0.180648160694857404058472031243;
    w[2] = 0.260610696402935462318742869419;
    w[3] = 0.312347077040002840068630406584;
    w[4] = 0.330239355001259763164525069287;
    w[5] = 0.312347077040002840068630406584;
    w[6] = 0.260610696402935462318742869419;
    w[7] = 0.180648160694857404058472031243;
    w[8] = 0.812743883615744119718921581105E-01;
  }
  else if ( n == 10 )
  {
    x[0] =  - 0.973906528517171720077964012084;
    x[1] =  - 0.865063366688984510732096688423;
    x[2] =  - 0.679409568299024406234327365115;
    x[3] =  - 0.433395394129247190799265943166;
    x[4] =  - 0.148874338981631210884826001130;
    x[5] =    0.148874338981631210884826001130;
    x[6] =    0.433395394129247190799265943166;
    x[7] =    0.679409568299024406234327365115;
    x[8] =    0.865063366688984510732096688423;
    x[9] =   0.973906528517171720077964012084;

    w[0] =  0.666713443086881375935688098933E-01;
    w[1] =  0.149451349150580593145776339658;
    w[2] =  0.219086362515982043995534934228;
    w[3] =  0.269266719309996355091226921569;
    w[4] =  0.295524224714752870173892994651;
    w[5] =  0.295524224714752870173892994651;
    w[6] =  0.269266719309996355091226921569;
    w[7] =  0.219086362515982043995534934228;
    w[8] =  0.149451349150580593145776339658;
    w[9] = 0.666713443086881375935688098933E-01;
  }
  else if ( n == 11 )
  {
    x[0] =  - 0.978228658146056992803938001123;
    x[1] =  - 0.887062599768095299075157769304;
    x[2] =  - 0.730152005574049324093416252031;
    x[3] =  - 0.519096129206811815925725669459;
    x[4] =  - 0.269543155952344972331531985401;
    x[5] =    0.0;
    x[6] =    0.269543155952344972331531985401;
    x[7] =    0.519096129206811815925725669459;
    x[8] =    0.730152005574049324093416252031;
    x[9] =   0.887062599768095299075157769304;
    x[10] =   0.978228658146056992803938001123;

    w[0] =  0.556685671161736664827537204425E-01;
    w[1] =  0.125580369464904624634694299224;
    w[2] =  0.186290210927734251426097641432;
    w[3] =  0.233193764591990479918523704843;
    w[4] =  0.262804544510246662180688869891;
    w[5] =  0.272925086777900630714483528336;
    w[6] =  0.262804544510246662180688869891;
    w[7] =  0.233193764591990479918523704843;
    w[8] =  0.186290210927734251426097641432;
    w[9] = 0.125580369464904624634694299224;
    w[10] = 0.556685671161736664827537204425E-01;
  }
  else if ( n == 12 )
  {
    x[0] =  - 0.981560634246719250690549090149;
    x[1] =  - 0.904117256370474856678465866119;
    x[2] =  - 0.769902674194304687036893833213;
    x[3] =  - 0.587317954286617447296702418941;
    x[4] =  - 0.367831498998180193752691536644;
    x[5] =  - 0.125233408511468915472441369464;
    x[6] =    0.125233408511468915472441369464;
    x[7] =    0.367831498998180193752691536644;
    x[8] =    0.587317954286617447296702418941;
    x[9] =   0.769902674194304687036893833213;
    x[10] =   0.904117256370474856678465866119;
    x[11] =   0.981560634246719250690549090149;

    w[0] =  0.471753363865118271946159614850E-01;
    w[1] =  0.106939325995318430960254718194;
    w[2] =  0.160078328543346226334652529543;
    w[3] =  0.203167426723065921749064455810;
    w[4] =  0.233492536538354808760849898925;
    w[5] =  0.249147045813402785000562436043;
    w[6] =  0.249147045813402785000562436043;
    w[7] =  0.233492536538354808760849898925;
    w[8] =  0.203167426723065921749064455810;
    w[9] = 0.160078328543346226334652529543;
    w[10] = 0.106939325995318430960254718194;
    w[11] = 0.471753363865118271946159614850E-01;
  }
  else if ( n == 13 )
  {
    x[0] =  - 0.984183054718588149472829448807;
    x[1] =  - 0.917598399222977965206547836501;
    x[2] =  - 0.801578090733309912794206489583;
    x[3] =  - 0.642349339440340220643984606996;
    x[4] =  - 0.448492751036446852877912852128;
    x[5] =  - 0.230458315955134794065528121098;
    x[6] =    0.0;
    x[7] =    0.230458315955134794065528121098;
    x[8] =    0.448492751036446852877912852128;
    x[9] =   0.642349339440340220643984606996;
    x[10] =   0.801578090733309912794206489583;
    x[11] =   0.917598399222977965206547836501;
    x[12] =   0.984183054718588149472829448807;

    w[0] =  0.404840047653158795200215922010E-01;
    w[1] =  0.921214998377284479144217759538E-01;
    w[2] =  0.138873510219787238463601776869;
    w[3] =  0.178145980761945738280046691996;
    w[4] =  0.207816047536888502312523219306;
    w[5] =  0.226283180262897238412090186040;
    w[6] =  0.232551553230873910194589515269;
    w[7] =  0.226283180262897238412090186040;
    w[8] =  0.207816047536888502312523219306;
    w[9] = 0.178145980761945738280046691996;
    w[10] = 0.138873510219787238463601776869;
    w[11] = 0.921214998377284479144217759538E-01;
    w[12] = 0.404840047653158795200215922010E-01;
  }
  else if ( n == 14 )
  {
    x[0] =  - 0.986283808696812338841597266704;
    x[1] =  - 0.928434883663573517336391139378;
    x[2] =  - 0.827201315069764993189794742650;
    x[3] =  - 0.687292904811685470148019803019;
    x[4] =  - 0.515248636358154091965290718551;
    x[5] =  - 0.319112368927889760435671824168;
    x[6] =  - 0.108054948707343662066244650220;
    x[7] =    0.108054948707343662066244650220;
    x[8] =    0.319112368927889760435671824168;
    x[9] =   0.515248636358154091965290718551;
    x[10] =   0.687292904811685470148019803019;
    x[11] =   0.827201315069764993189794742650;
    x[12] =   0.928434883663573517336391139378;
    x[13] =   0.986283808696812338841597266704;

    w[0] =  0.351194603317518630318328761382E-01;
    w[1] =  0.801580871597602098056332770629E-01;
    w[2] =  0.121518570687903184689414809072;
    w[3] =  0.157203167158193534569601938624;
    w[4] =  0.185538397477937813741716590125;
    w[5] =  0.205198463721295603965924065661;
    w[6] =  0.215263853463157790195876443316;
    w[7] =  0.215263853463157790195876443316;
    w[8] =  0.205198463721295603965924065661;
    w[9] = 0.185538397477937813741716590125;
    w[10] = 0.157203167158193534569601938624;
    w[11] = 0.121518570687903184689414809072;
    w[12] = 0.801580871597602098056332770629E-01;
    w[13] = 0.351194603317518630318328761382E-01;
  }
  else if ( n == 15 )
  {
    x[0] =  - 0.987992518020485428489565718587;
    x[1] =  - 0.937273392400705904307758947710;
    x[2] =  - 0.848206583410427216200648320774;
    x[3] =  - 0.724417731360170047416186054614;
    x[4] =  - 0.570972172608538847537226737254;
    x[5] =  - 0.394151347077563369897207370981;
    x[6] =  - 0.201194093997434522300628303395;
    x[7] =    0.0;
    x[8] =    0.201194093997434522300628303395;
    x[9] =   0.394151347077563369897207370981;
    x[10] =   0.570972172608538847537226737254;
    x[11] =   0.724417731360170047416186054614;
    x[12] =   0.848206583410427216200648320774;
    x[13] =   0.937273392400705904307758947710;
    x[14] =   0.987992518020485428489565718587;

    w[0] =  0.307532419961172683546283935772E-01;
    w[1] =  0.703660474881081247092674164507E-01;
    w[2] =  0.107159220467171935011869546686;
    w[3] =  0.139570677926154314447804794511;
    w[4] =  0.166269205816993933553200860481;
    w[5] =  0.186161000015562211026800561866;
    w[6] =  0.198431485327111576456118326444;
    w[7] =  0.202578241925561272880620199968;
    w[8] =  0.198431485327111576456118326444;
    w[9] = 0.186161000015562211026800561866;
    w[10] = 0.166269205816993933553200860481;
    w[11] = 0.139570677926154314447804794511;
    w[12] = 0.107159220467171935011869546686;
    w[13] = 0.703660474881081247092674164507E-01;
    w[14] = 0.307532419961172683546283935772E-01;
  }
  else if ( n == 16 )
  {
    x[0] =  - 0.989400934991649932596154173450;
    x[1] =  - 0.944575023073232576077988415535;
    x[2] =  - 0.865631202387831743880467897712;
    x[3] =  - 0.755404408355003033895101194847;
    x[4] =  - 0.617876244402643748446671764049;
    x[5] =  - 0.458016777657227386342419442984;
    x[6] =  - 0.281603550779258913230460501460;
    x[7] =  - 0.950125098376374401853193354250E-01;
    x[8] =    0.950125098376374401853193354250E-01;
    x[9] =   0.281603550779258913230460501460;
    x[10] =   0.458016777657227386342419442984;
    x[11] =   0.617876244402643748446671764049;
    x[12] =   0.755404408355003033895101194847;
    x[13] =   0.865631202387831743880467897712;
    x[14] =   0.944575023073232576077988415535;
    x[15] =   0.989400934991649932596154173450;

    w[0] =  0.271524594117540948517805724560E-01;
    w[1] =  0.622535239386478928628438369944E-01;
    w[2] =  0.951585116824927848099251076022E-01;
    w[3] =  0.124628971255533872052476282192;
    w[4] =  0.149595988816576732081501730547;
    w[5] =  0.169156519395002538189312079030;
    w[6] =  0.182603415044923588866763667969;
    w[7] =  0.189450610455068496285396723208;
    w[8] =  0.189450610455068496285396723208;
    w[9] = 0.182603415044923588866763667969;
    w[10] = 0.169156519395002538189312079030;
    w[11] = 0.149595988816576732081501730547;
    w[12] = 0.124628971255533872052476282192;
    w[13] = 0.951585116824927848099251076022E-01;
    w[14] = 0.622535239386478928628438369944E-01;
    w[15] = 0.271524594117540948517805724560E-01;
  }
  else if ( n == 17 )
  {
    x[0] =  - 0.990575475314417335675434019941;
    x[1] =  - 0.950675521768767761222716957896;
    x[2] =  - 0.880239153726985902122955694488;
    x[3] =  - 0.781514003896801406925230055520;
    x[4] =  - 0.657671159216690765850302216643;
    x[5] =  - 0.512690537086476967886246568630;
    x[6] =  - 0.351231763453876315297185517095;
    x[7] =  - 0.178484181495847855850677493654;
    x[8] =    0.0;
    x[9] =   0.178484181495847855850677493654;
    x[10] =   0.351231763453876315297185517095;
    x[11] =   0.512690537086476967886246568630;
    x[12] =   0.657671159216690765850302216643;
    x[13] =   0.781514003896801406925230055520;
    x[14] =   0.880239153726985902122955694488;
    x[15] =   0.950675521768767761222716957896;
    x[16] =   0.990575475314417335675434019941;

    w[0] =  0.241483028685479319601100262876E-01;
    w[1] =  0.554595293739872011294401653582E-01;
    w[2] =  0.850361483171791808835353701911E-01;
    w[3] =  0.111883847193403971094788385626;
    w[4] =  0.135136368468525473286319981702;
    w[5] =  0.154045761076810288081431594802;
    w[6] =  0.168004102156450044509970663788;
    w[7] =  0.176562705366992646325270990113;
    w[8] =  0.179446470356206525458265644262;
    w[9] = 0.176562705366992646325270990113;
    w[10] = 0.168004102156450044509970663788;
    w[11] = 0.154045761076810288081431594802;
    w[12] = 0.135136368468525473286319981702;
    w[13] = 0.111883847193403971094788385626;
    w[14] = 0.850361483171791808835353701911E-01;
    w[15] = 0.554595293739872011294401653582E-01;
    w[16] = 0.241483028685479319601100262876E-01;
  }
  else if ( n == 18 )
  {
    x[0] =  - 0.991565168420930946730016004706;
    x[1] =  - 0.955823949571397755181195892930;
    x[2] =  - 0.892602466497555739206060591127;
    x[3] =  - 0.803704958972523115682417455015;
    x[4] =  - 0.691687043060353207874891081289;
    x[5] =  - 0.559770831073947534607871548525;
    x[6] =  - 0.411751161462842646035931793833;
    x[7] =  - 0.251886225691505509588972854878;
    x[8] =  - 0.847750130417353012422618529358E-01;
    x[9] =   0.847750130417353012422618529358E-01;
    x[10] =   0.251886225691505509588972854878;
    x[11] =   0.411751161462842646035931793833;
    x[12] =   0.559770831073947534607871548525;
    x[13] =   0.691687043060353207874891081289;
    x[14] =   0.803704958972523115682417455015;
    x[15] =   0.892602466497555739206060591127;
    x[16] =   0.955823949571397755181195892930;
    x[17] =   0.991565168420930946730016004706;

    w[0] =  0.216160135264833103133427102665E-01;
    w[1] =  0.497145488949697964533349462026E-01;
    w[2] =  0.764257302548890565291296776166E-01;
    w[3] =  0.100942044106287165562813984925;
    w[4] =  0.122555206711478460184519126800;
    w[5] =  0.140642914670650651204731303752;
    w[6] =  0.154684675126265244925418003836;
    w[7] =  0.164276483745832722986053776466;
    w[8] =  0.169142382963143591840656470135;
    w[9] = 0.169142382963143591840656470135;
    w[10] = 0.164276483745832722986053776466;
    w[11] = 0.154684675126265244925418003836;
    w[12] = 0.140642914670650651204731303752;
    w[13] = 0.122555206711478460184519126800;
    w[14] = 0.100942044106287165562813984925;
    w[15] = 0.764257302548890565291296776166E-01;
    w[16] = 0.497145488949697964533349462026E-01;
    w[17] = 0.216160135264833103133427102665E-01;
  }
  else if ( n == 19 )
  {
    x[0] =  - 0.992406843843584403189017670253;
    x[1] =  - 0.960208152134830030852778840688;
    x[2] =  - 0.903155903614817901642660928532;
    x[3] =  - 0.822714656537142824978922486713;
    x[4] =  - 0.720966177335229378617095860824;
    x[5] =  - 0.600545304661681023469638164946;
    x[6] =  - 0.464570741375960945717267148104;
    x[7] =  - 0.316564099963629831990117328850;
    x[8] =  - 0.160358645640225375868096115741;
    x[9] =   0.0;
    x[10] =   0.160358645640225375868096115741;
    x[11] =   0.316564099963629831990117328850;
    x[12] =   0.464570741375960945717267148104;
    x[13] =   0.600545304661681023469638164946;
    x[14] =   0.720966177335229378617095860824;
    x[15] =   0.822714656537142824978922486713;
    x[16] =   0.903155903614817901642660928532;
    x[17] =   0.960208152134830030852778840688;
    x[18] =   0.992406843843584403189017670253;

    w[0] =  0.194617882297264770363120414644E-01;
    w[1] =  0.448142267656996003328381574020E-01;
    w[2] =  0.690445427376412265807082580060E-01;
    w[3] =  0.914900216224499994644620941238E-01;
    w[4] =  0.111566645547333994716023901682;
    w[5] =  0.128753962539336227675515784857;
    w[6] =  0.142606702173606611775746109442;
    w[7] =  0.152766042065859666778855400898;
    w[8] =  0.158968843393954347649956439465;
    w[9] = 0.161054449848783695979163625321;
    w[10] = 0.158968843393954347649956439465;
    w[11] = 0.152766042065859666778855400898;
    w[12] = 0.142606702173606611775746109442;
    w[13] = 0.128753962539336227675515784857;
    w[14] = 0.111566645547333994716023901682;
    w[15] = 0.914900216224499994644620941238E-01;
    w[16] = 0.690445427376412265807082580060E-01;
    w[17] = 0.448142267656996003328381574020E-01;
    w[18] = 0.194617882297264770363120414644E-01;
  }
  else if ( n == 20 )
  {
    x[0] =  - 0.993128599185094924786122388471;
    x[1] =  - 0.963971927277913791267666131197;
    x[2] =  - 0.912234428251325905867752441203;
    x[3] =  - 0.839116971822218823394529061702;
    x[4] =  - 0.746331906460150792614305070356;
    x[5] =  - 0.636053680726515025452836696226;
    x[6] =  - 0.510867001950827098004364050955;
    x[7] =  - 0.373706088715419560672548177025;
    x[8] =  - 0.227785851141645078080496195369;
    x[9] = - 0.765265211334973337546404093988E-01;
    x[10] =   0.765265211334973337546404093988E-01;
    x[11] =   0.227785851141645078080496195369;
    x[12] =   0.373706088715419560672548177025;
    x[13] =   0.510867001950827098004364050955;
    x[14] =   0.636053680726515025452836696226;
    x[15] =   0.746331906460150792614305070356;
    x[16] =   0.839116971822218823394529061702;
    x[17] =   0.912234428251325905867752441203;
    x[18] =   0.963971927277913791267666131197;
    x[19] =   0.993128599185094924786122388471;

    w[0] =  0.176140071391521183118619623519E-01;
    w[1] =  0.406014298003869413310399522749E-01;
    w[2] =  0.626720483341090635695065351870E-01;
    w[3] =  0.832767415767047487247581432220E-01;
    w[4] =  0.101930119817240435036750135480;
    w[5] =  0.118194531961518417312377377711;
    w[6] =  0.131688638449176626898494499748;
    w[7] =  0.142096109318382051329298325067;
    w[8] =  0.149172986472603746787828737002;
    w[9] = 0.152753387130725850698084331955;
    w[10] = 0.152753387130725850698084331955;
    w[11] = 0.149172986472603746787828737002;
    w[12] = 0.142096109318382051329298325067;
    w[13] = 0.131688638449176626898494499748;
    w[14] = 0.118194531961518417312377377711;
    w[15] = 0.101930119817240435036750135480;
    w[16] = 0.832767415767047487247581432220E-01;
    w[17] = 0.626720483341090635695065351870E-01;
    w[18] = 0.406014298003869413310399522749E-01;
    w[19] = 0.176140071391521183118619623519E-01;
  }
  else if ( n == 21 )
  {
    x[ 0] =  -0.9937521706203896E+00;
    x[ 1] =  -0.9672268385663063E+00;
    x[ 2] =  -0.9200993341504008E+00;
    x[ 3] =  -0.8533633645833173E+00;
    x[ 4] =  -0.7684399634756779E+00;
    x[ 5] =  -0.6671388041974123E+00;
    x[ 6] =  -0.5516188358872198E+00;
    x[ 7] =  -0.4243421202074388E+00;
    x[ 8] =  -0.2880213168024011E+00;
    x[ 9] =  -0.1455618541608951E+00;
    x[10] =   0.0000000000000000E+00;
    x[11] =   0.1455618541608951E+00;
    x[12] =   0.2880213168024011E+00;
    x[13] =   0.4243421202074388E+00;
    x[14] =   0.5516188358872198E+00;
    x[15] =   0.6671388041974123E+00;
    x[16] =   0.7684399634756779E+00;
    x[17] =   0.8533633645833173E+00;
    x[18] =   0.9200993341504008E+00;
    x[19] =   0.9672268385663063E+00;
    x[20] =   0.9937521706203896E+00;

    w[ 0] =   0.1601722825777420E-01;
    w[ 1] =   0.3695378977085242E-01;
    w[ 2] =   0.5713442542685715E-01;
    w[ 3] =   0.7610011362837928E-01;
    w[ 4] =   0.9344442345603393E-01;
    w[ 5] =   0.1087972991671484E+00;
    w[ 6] =   0.1218314160537285E+00;
    w[ 7] =   0.1322689386333373E+00;
    w[ 8] =   0.1398873947910731E+00;
    w[ 9] =   0.1445244039899700E+00;
    w[10] =   0.1460811336496904E+00;
    w[11] =   0.1445244039899700E+00;
    w[12] =   0.1398873947910731E+00;
    w[13] =   0.1322689386333373E+00;
    w[14] =   0.1218314160537285E+00;
    w[15] =   0.1087972991671484E+00;
    w[16] =   0.9344442345603393E-01;
    w[17] =   0.7610011362837928E-01;
    w[18] =   0.5713442542685715E-01;
    w[19] =   0.3695378977085242E-01;
    w[20] =   0.1601722825777420E-01;
  }
  else if ( n == 22 )
  {
    x[ 0] =  -0.9942945854823994E+00;
    x[ 1] =  -0.9700604978354287E+00;
    x[ 2] =  -0.9269567721871740E+00;
    x[ 3] =  -0.8658125777203002E+00;
    x[ 4] =  -0.7878168059792081E+00;
    x[ 5] =  -0.6944872631866827E+00;
    x[ 6] =  -0.5876404035069116E+00;
    x[ 7] =  -0.4693558379867570E+00;
    x[ 8] =  -0.3419358208920842E+00;
    x [9] =  -0.2078604266882213E+00;
    x[10] =  -0.6973927331972223E-01;
    x[11] =   0.6973927331972223E-01;
    x[12] =   0.2078604266882213E+00;
    x[13] =   0.3419358208920842E+00;
    x[14] =   0.4693558379867570E+00;
    x[15] =   0.5876404035069116E+00;
    x[16] =   0.6944872631866827E+00;
    x[17] =   0.7878168059792081E+00;
    x[18] =   0.8658125777203002E+00;
    x[19] =   0.9269567721871740E+00;
    x[20] =   0.9700604978354287E+00;
    x[21] =   0.9942945854823994E+00;

    w[ 0] =   0.1462799529827203E-01;
    w[ 1] =   0.3377490158481413E-01;
    w[ 2] =   0.5229333515268327E-01;
    w[ 3] =   0.6979646842452038E-01;
    w[ 4] =   0.8594160621706777E-01;
    w[ 5] =   0.1004141444428809E+00;
    w[ 6] =   0.1129322960805392E+00;
    w[ 7] =   0.1232523768105124E+00;
    w[ 8] =   0.1311735047870623E+00;
    w[ 9] =   0.1365414983460152E+00;
    w[10] =   0.1392518728556321E+00;
    w[11] =   0.1392518728556321E+00;
    w[12] =   0.1365414983460152E+00;
    w[13] =   0.1311735047870623E+00;
    w[14] =   0.1232523768105124E+00;
    w[15] =   0.1129322960805392E+00;
    w[16] =   0.1004141444428809E+00;
    w[17] =   0.8594160621706777E-01;
    w[18] =   0.6979646842452038E-01;
    w[19] =   0.5229333515268327E-01;
    w[20] =   0.3377490158481413E-01;
    w[21] =   0.1462799529827203E-01;
  }
  else if ( n == 23 )
  {
    x[ 0] =  -0.9947693349975522E+00;
    x[ 1] =  -0.9725424712181152E+00;
    x[ 2] =  -0.9329710868260161E+00;
    x[ 3] =  -0.8767523582704416E+00;
    x[ 4] =  -0.8048884016188399E+00;
    x[ 5] =  -0.7186613631319502E+00;
    x[ 6] =  -0.6196098757636461E+00;
    x[ 7] =  -0.5095014778460075E+00;
    x[ 8] =  -0.3903010380302908E+00;
    x[ 9] =  -0.2641356809703449E+00;
    x[10] =  -0.1332568242984661E+00;
    x[11] =   0.0000000000000000E+00;
    x[12] =   0.1332568242984661E+00;
    x[13] =   0.2641356809703449E+00;
    x[14] =   0.3903010380302908E+00;
    x[15] =   0.5095014778460075E+00;
    x[16] =   0.6196098757636461E+00;
    x[17] =   0.7186613631319502E+00;
    x[18] =   0.8048884016188399E+00;
    x[19] =   0.8767523582704416E+00;
    x[20] =   0.9329710868260161E+00;
    x[21] =   0.9725424712181152E+00;
    x[22] =   0.9947693349975522E+00;

    w[ 0] =   0.1341185948714167E-01;
    w[ 1] =   0.3098800585697944E-01;
    w[ 2] =   0.4803767173108464E-01;
    w[ 3] =   0.6423242140852586E-01;
    w[ 4] =   0.7928141177671895E-01;
    w[ 5] =   0.9291576606003514E-01;
    w[ 6] =   0.1048920914645414E+00;
    w[ 7] =   0.1149966402224114E+00;
    w[ 8] =   0.1230490843067295E+00;
    w[ 9] =   0.1289057221880822E+00;
    w[10] =   0.1324620394046967E+00;
    w[11] =   0.1336545721861062E+00;
    w[12] =   0.1324620394046967E+00;
    w[13] =   0.1289057221880822E+00;
    w[14] =   0.1230490843067295E+00;
    w[15] =   0.1149966402224114E+00;
    w[16] =   0.1048920914645414E+00;
    w[17] =   0.9291576606003514E-01;
    w[18] =   0.7928141177671895E-01;
    w[19] =   0.6423242140852586E-01;
    w[20] =   0.4803767173108464E-01;
    w[21] =   0.3098800585697944E-01;
    w[22] =   0.1341185948714167E-01;
  }
  else if ( n == 24 )
  {
    x[ 0] =  -0.9951872199970213E+00;
    x[ 1] =  -0.9747285559713095E+00;
    x[ 2] =  -0.9382745520027327E+00;
    x[ 3] =  -0.8864155270044011E+00;
    x[ 4] =  -0.8200019859739029E+00;
    x[ 5] =  -0.7401241915785544E+00;
    x[ 6] =  -0.6480936519369755E+00;
    x[ 7] =  -0.5454214713888396E+00;
    x[ 8] =  -0.4337935076260451E+00;
    x[ 9] =  -0.3150426796961634E+00;
    x[10] =  -0.1911188674736163E+00;
    x[11] =  -0.6405689286260562E-01;
    x[12] =   0.6405689286260562E-01;
    x[13] =   0.1911188674736163E+00;
    x[14] =   0.3150426796961634E+00;
    x[15] =   0.4337935076260451E+00;
    x[16] =   0.5454214713888396E+00;
    x[17] =   0.6480936519369755E+00;
    x[18] =   0.7401241915785544E+00;
    x[19] =   0.8200019859739029E+00;
    x[20] =   0.8864155270044011E+00;
    x[21] =   0.9382745520027327E+00;
    x[22] =   0.9747285559713095E+00;
    x[23] =   0.9951872199970213E+00;

    w[ 0] =   0.1234122979998730E-01;
    w[ 1] =   0.2853138862893375E-01;
    w[ 2] =   0.4427743881741982E-01;
    w[ 3] =   0.5929858491543672E-01;
    w[ 4] =   0.7334648141108031E-01;
    w[ 5] =   0.8619016153195320E-01;
    w[ 6] =   0.9761865210411380E-01;
    w[ 7] =   0.1074442701159656E+00;
    w[ 8] =   0.1155056680537256E+00;
    w[ 9] =   0.1216704729278035E+00;
    w[10] =   0.1258374563468283E+00;
    w[11] =   0.1279381953467521E+00;
    w[12] =   0.1279381953467521E+00;
    w[13] =   0.1258374563468283E+00;
    w[14] =   0.1216704729278035E+00;
    w[15] =   0.1155056680537256E+00;
    w[16] =   0.1074442701159656E+00;
    w[17] =   0.9761865210411380E-01;
    w[18] =   0.8619016153195320E-01;
    w[19] =   0.7334648141108031E-01;
    w[20] =   0.5929858491543672E-01;
    w[21] =   0.4427743881741982E-01;
    w[22] =   0.2853138862893375E-01;
    w[23] =   0.1234122979998730E-01;
  }
  else if ( n == 25 )
  {
    x[ 0] =  -0.9955569697904981E+00;
    x[ 1] =  -0.9766639214595175E+00;
    x[ 2] =  -0.9429745712289743E+00;
    x[ 3] =  -0.8949919978782754E+00;
    x[ 4] =  -0.8334426287608340E+00;
    x[ 5] =  -0.7592592630373577E+00;
    x[ 6] =  -0.6735663684734684E+00;
    x[ 7] =  -0.5776629302412229E+00;
    x[ 8] =  -0.4730027314457150E+00;
    x[ 9] =  -0.3611723058093879E+00;
    x[10] =  -0.2438668837209884E+00;
    x[11] =  -0.1228646926107104E+00;
    x[12] =   0.0000000000000000E+00;
    x[13] =   0.1228646926107104E+00;
    x[14] =   0.2438668837209884E+00;
    x[15] =   0.3611723058093879E+00;
    x[16] =   0.4730027314457150E+00;
    x[17] =   0.5776629302412229E+00;
    x[18] =   0.6735663684734684E+00;
    x[19] =   0.7592592630373577E+00;
    x[20] =   0.8334426287608340E+00;
    x[21] =   0.8949919978782754E+00;
    x[22] =   0.9429745712289743E+00;
    x[23] =   0.9766639214595175E+00;
    x[24] =   0.9955569697904981E+00;

    w[ 0] =   0.1139379850102617E-01;
    w[ 1] =   0.2635498661503214E-01;
    w[ 2] =   0.4093915670130639E-01;
    w[ 3] =   0.5490469597583517E-01;
    w[ 4] =   0.6803833381235694E-01;
    w[ 5] =   0.8014070033500101E-01;
    w[ 6] =   0.9102826198296370E-01;
    w[ 7] =   0.1005359490670506E+00;
    w[ 8] =   0.1085196244742637E+00;
    w[ 9] =   0.1148582591457116E+00;
    w[10] =   0.1194557635357847E+00;
    w[11] =   0.1222424429903101E+00;
    w[12] =   0.1231760537267154E+00;
    w[13] =   0.1222424429903101E+00;
    w[14] =   0.1194557635357847E+00;
    w[15] =   0.1148582591457116E+00;
    w[16] =   0.1085196244742637E+00;
    w[17] =   0.1005359490670506E+00;
    w[18] =   0.9102826198296370E-01;
    w[19] =   0.8014070033500101E-01;
    w[20] =   0.6803833381235694E-01;
    w[21] =   0.5490469597583517E-01;
    w[22] =   0.4093915670130639E-01;
    w[23] =   0.2635498661503214E-01;
    w[24] =   0.1139379850102617E-01;
  }
  else if ( n == 26 )
  {
    x[ 0] =  -0.9958857011456169E+00;
    x[ 1] =  -0.9783854459564710E+00;
    x[ 2] =  -0.9471590666617142E+00;
    x[ 3] =  -0.9026378619843071E+00;
    x[ 4] =  -0.8454459427884981E+00;
    x[ 5] =  -0.7763859488206789E+00;
    x[ 6] =  -0.6964272604199573E+00;
    x[ 7] =  -0.6066922930176181E+00;
    x[ 8] =  -0.5084407148245057E+00;
    x[ 9] =  -0.4030517551234863E+00;
    x[10] =  -0.2920048394859569E+00;
    x[11] =  -0.1768588203568902E+00;
    x[12] =  -0.5923009342931320E-01;
    x[13] =   0.5923009342931320E-01;
    x[14] =   0.1768588203568902E+00;
    x[15] =   0.2920048394859569E+00;
    x[16] =   0.4030517551234863E+00;
    x[17] =   0.5084407148245057E+00;
    x[18] =   0.6066922930176181E+00;
    x[19] =   0.6964272604199573E+00;
    x[20] =   0.7763859488206789E+00;
    x[21] =   0.8454459427884981E+00;
    x[22] =   0.9026378619843071E+00;
    x[23] =   0.9471590666617142E+00;
    x[24] =   0.9783854459564710E+00;
    x[25] =   0.9958857011456169E+00;

    w[ 0] =   0.1055137261734304E-01;
    w[ 1] =   0.2441785109263173E-01;
    w[ 2] =   0.3796238329436282E-01;
    w[ 3] =   0.5097582529714782E-01;
    w[ 4] =   0.6327404632957484E-01;
    w[ 5] =   0.7468414976565967E-01;
    w[ 6] =   0.8504589431348521E-01;
    w[ 7] =   0.9421380035591416E-01;
    w[ 8] =   0.1020591610944255E+00;
    w[ 9] =   0.1084718405285765E+00;
    w[10] =   0.1133618165463197E+00;
    w[11] =   0.1166604434852967E+00;
    w[12] =   0.1183214152792622E+00;
    w[13] =   0.1183214152792622E+00;
    w[14] =   0.1166604434852967E+00;
    w[15] =   0.1133618165463197E+00;
    w[16] =   0.1084718405285765E+00;
    w[17] =   0.1020591610944255E+00;
    w[18] =   0.9421380035591416E-01;
    w[19] =   0.8504589431348521E-01;
    w[20] =   0.7468414976565967E-01;
    w[21] =   0.6327404632957484E-01;
    w[22] =   0.5097582529714782E-01;
    w[23] =   0.3796238329436282E-01;
    w[24] =   0.2441785109263173E-01;
    w[25] =   0.1055137261734304E-01;
  }
  else if ( n == 27 )
  {
    x[ 0] =  -0.9961792628889886E+00;
    x[ 1] =  -0.9799234759615012E+00;
    x[ 2] =  -0.9509005578147051E+00;
    x[ 3] =  -0.9094823206774911E+00;
    x[ 4] =  -0.8562079080182945E+00;
    x[ 5] =  -0.7917716390705082E+00;
    x[ 6] =  -0.7170134737394237E+00;
    x[ 7] =  -0.6329079719464952E+00;
    x[ 8] =  -0.5405515645794569E+00;
    x[ 9] =  -0.4411482517500269E+00;
    x[10] =  -0.3359939036385089E+00;
    x[11] =  -0.2264593654395369E+00;
    x[12] =  -0.1139725856095300E+00;
    x[13] =   0.0000000000000000E+00;
    x[14] =   0.1139725856095300E+00;
    x[15] =   0.2264593654395369E+00;
    x[16] =   0.3359939036385089E+00;
    x[17] =   0.4411482517500269E+00;
    x[18] =   0.5405515645794569E+00;
    x[19] =   0.6329079719464952E+00;
    x[20] =   0.7170134737394237E+00;
    x[21] =   0.7917716390705082E+00;
    x[22] =   0.8562079080182945E+00;
    x[23] =   0.9094823206774911E+00;
    x[24] =   0.9509005578147051E+00;
    x[25] =   0.9799234759615012E+00;
    x[26] =   0.9961792628889886E+00;

    w[ 0] =   0.9798996051294232E-02;
    w[ 1] =   0.2268623159618062E-01;
    w[ 2] =   0.3529705375741969E-01;
    w[ 3] =   0.4744941252061504E-01;
    w[ 4] =   0.5898353685983366E-01;
    w[ 5] =   0.6974882376624561E-01;
    w[ 6] =   0.7960486777305781E-01;
    w[ 7] =   0.8842315854375689E-01;
    w[ 8] =   0.9608872737002842E-01;
    w[ 9] =   0.1025016378177459E+00;
    w[10] =   0.1075782857885332E+00;
    w[11] =   0.1112524883568452E+00;
    w[12] =   0.1134763461089651E+00;
    w[13] =   0.1142208673789570E+00;
    w[14] =   0.1134763461089651E+00;
    w[15] =   0.1112524883568452E+00;
    w[16] =   0.1075782857885332E+00;
    w[17] =   0.1025016378177459E+00;
    w[18] =   0.9608872737002842E-01;
    w[19] =   0.8842315854375689E-01;
    w[20] =   0.7960486777305781E-01;
    w[21] =   0.6974882376624561E-01;
    w[22] =   0.5898353685983366E-01;
    w[23] =   0.4744941252061504E-01;
    w[24] =   0.3529705375741969E-01;
    w[25] =   0.2268623159618062E-01;
    w[26] =   0.9798996051294232E-02;
  }
  else if ( n == 28 )
  {
    x[ 0] =  -0.9964424975739544E+00;
    x[ 1] =  -0.9813031653708728E+00;
    x[ 2] =  -0.9542592806289382E+00;
    x[ 3] =  -0.9156330263921321E+00;
    x[ 4] =  -0.8658925225743951E+00;
    x[ 5] =  -0.8056413709171791E+00;
    x[ 6] =  -0.7356108780136318E+00;
    x[ 7] =  -0.6566510940388650E+00;
    x[ 8] =  -0.5697204718114017E+00;
    x[ 9] =  -0.4758742249551183E+00;
    x[10] =  -0.3762515160890787E+00;
    x[11] =  -0.2720616276351780E+00;
    x[12] =  -0.1645692821333808E+00;
    x[13] =  -0.5507928988403427E-01;
    x[14] =   0.5507928988403427E-01;
    x[15] =   0.1645692821333808E+00;
    x[16] =   0.2720616276351780E+00;
    x[17] =   0.3762515160890787E+00;
    x[18] =   0.4758742249551183E+00;
    x[19] =   0.5697204718114017E+00;
    x[20] =   0.6566510940388650E+00;
    x[21] =   0.7356108780136318E+00;
    x[22] =   0.8056413709171791E+00;
    x[23] =   0.8658925225743951E+00;
    x[24] =   0.9156330263921321E+00;
    x[25] =   0.9542592806289382E+00;
    x[26] =   0.9813031653708728E+00;
    x[27] =   0.9964424975739544E+00;

    w[ 0] =   0.9124282593094672E-02;
    w[ 1] =   0.2113211259277118E-01;
    w[ 2] =   0.3290142778230441E-01;
    w[ 3] =   0.4427293475900429E-01;
    w[ 4] =   0.5510734567571667E-01;
    w[ 5] =   0.6527292396699959E-01;
    w[ 6] =   0.7464621423456877E-01;
    w[ 7] =   0.8311341722890127E-01;
    w[ 8] =   0.9057174439303289E-01;
    w[ 9] =   0.9693065799792999E-01;
    w[10] =   0.1021129675780608E+00;
    w[11] =   0.1060557659228464E+00;
    w[12] =   0.1087111922582942E+00;
    w[13] =   0.1100470130164752E+00;
    w[14] =   0.1100470130164752E+00;
    w[15] =   0.1087111922582942E+00;
    w[16] =   0.1060557659228464E+00;
    w[17] =   0.1021129675780608E+00;
    w[18] =   0.9693065799792999E-01;
    w[19] =   0.9057174439303289E-01;
    w[20] =   0.8311341722890127E-01;
    w[21] =   0.7464621423456877E-01;
    w[22] =   0.6527292396699959E-01;
    w[23] =   0.5510734567571667E-01;
    w[24] =   0.4427293475900429E-01;
    w[25] =   0.3290142778230441E-01;
    w[26] =   0.2113211259277118E-01;
    w[27] =   0.9124282593094672E-02;
  }
  else if ( n == 29 )
  {
    x[ 0] =  -0.9966794422605966E+00;
    x[ 1] =  -0.9825455052614132E+00;
    x[ 2] =  -0.9572855957780877E+00;
    x[ 3] =  -0.9211802329530588E+00;
    x[ 4] =  -0.8746378049201028E+00;
    x[ 5] =  -0.8181854876152524E+00;
    x[ 6] =  -0.7524628517344771E+00;
    x[ 7] =  -0.6782145376026865E+00;
    x[ 8] =  -0.5962817971382278E+00;
    x[ 9] =  -0.5075929551242276E+00;
    x[10] =  -0.4131528881740087E+00;
    x[11] =  -0.3140316378676399E+00;
    x[12] =  -0.2113522861660011E+00;
    x[13] =  -0.1062782301326792E+00;
    x[14] =   0.0000000000000000E+00;
    x[15] =   0.1062782301326792E+00;
    x[16] =   0.2113522861660011E+00;
    x[17] =   0.3140316378676399E+00;
    x[18] =   0.4131528881740087E+00;
    x[19] =   0.5075929551242276E+00;
    x[20] =   0.5962817971382278E+00;
    x[21] =   0.6782145376026865E+00;
    x[22] =   0.7524628517344771E+00;
    x[23] =   0.8181854876152524E+00;
    x[24] =   0.8746378049201028E+00;
    x[25] =   0.9211802329530588E+00;
    x[26] =   0.9572855957780877E+00;
    x[27] =   0.9825455052614132E+00;
    x[28] =   0.9966794422605966E+00;

    w[ 0] =   0.8516903878746365E-02;
    w[ 1] =   0.1973208505612276E-01;
    w[ 2] =   0.3074049220209360E-01;
    w[ 3] =   0.4140206251868281E-01;
    w[ 4] =   0.5159482690249799E-01;
    w[ 5] =   0.6120309065707916E-01;
    w[ 6] =   0.7011793325505125E-01;
    w[ 7] =   0.7823832713576385E-01;
    w[ 8] =   0.8547225736617248E-01;
    w[ 9] =   0.9173775713925882E-01;
    w[10] =   0.9696383409440862E-01;
    w[11] =   0.1010912737599150E+00;
    w[12] =   0.1040733100777293E+00;
    w[13] =   0.1058761550973210E+00;
    w[14] =   0.1064793817183143E+00;
    w[15] =   0.1058761550973210E+00;
    w[16] =   0.1040733100777293E+00;
    w[17] =   0.1010912737599150E+00;
    w[18] =   0.9696383409440862E-01;
    w[19] =   0.9173775713925882E-01;
    w[20] =   0.8547225736617248E-01;
    w[21] =   0.7823832713576385E-01;
    w[22] =   0.7011793325505125E-01;
    w[23] =   0.6120309065707916E-01;
    w[24] =   0.5159482690249799E-01;
    w[25] =   0.4140206251868281E-01;
    w[26] =   0.3074049220209360E-01;
    w[27] =   0.1973208505612276E-01;
    w[28] =   0.8516903878746365E-02;
  }
  else if ( n == 30 )
  {
    x[ 0] =  -0.9968934840746495E+00;
    x[ 1] =  -0.9836681232797472E+00;
    x[ 2] =  -0.9600218649683075E+00;
    x[ 3] =  -0.9262000474292743E+00;
    x[ 4] =  -0.8825605357920526E+00;
    x[ 5] =  -0.8295657623827684E+00;
    x[ 6] =  -0.7677774321048262E+00;
    x[ 7] =  -0.6978504947933158E+00;
    x[ 8] =  -0.6205261829892429E+00;
    x[ 9] =  -0.5366241481420199E+00;
    x[10] =  -0.4470337695380892E+00;
    x[11] =  -0.3527047255308781E+00;
    x[12] =  -0.2546369261678899E+00;
    x[13] =  -0.1538699136085835E+00;
    x[14] =  -0.5147184255531770E-01;
    x[15] =   0.5147184255531770E-01;
    x[16] =   0.1538699136085835E+00;
    x[17] =   0.2546369261678899E+00;
    x[18] =   0.3527047255308781E+00;
    x[19] =   0.4470337695380892E+00;
    x[20] =   0.5366241481420199E+00;
    x[21] =   0.6205261829892429E+00;
    x[22] =   0.6978504947933158E+00;
    x[23] =   0.7677774321048262E+00;
    x[24] =   0.8295657623827684E+00;
    x[25] =   0.8825605357920526E+00;
    x[26] =   0.9262000474292743E+00;
    x[27] =   0.9600218649683075E+00;
    x[28] =   0.9836681232797472E+00;
    x[29] =   0.9968934840746495E+00;

    w[ 0] =   0.7968192496166648E-02;
    w[ 1] =   0.1846646831109099E-01;
    w[ 2] =   0.2878470788332330E-01;
    w[ 3] =   0.3879919256962704E-01;
    w[ 4] =   0.4840267283059405E-01;
    w[ 5] =   0.5749315621761905E-01;
    w[ 6] =   0.6597422988218052E-01;
    w[ 7] =   0.7375597473770516E-01;
    w[ 8] =   0.8075589522942023E-01;
    w[ 9] =   0.8689978720108314E-01;
    w[10] =   0.9212252223778619E-01;
    w[11] =   0.9636873717464424E-01;
    w[12] =   0.9959342058679524E-01;
    w[13] =   0.1017623897484056E+00;
    w[14] =   0.1028526528935587E+00;
    w[15] =   0.1028526528935587E+00;
    w[16] =   0.1017623897484056E+00;
    w[17] =   0.9959342058679524E-01;
    w[18] =   0.9636873717464424E-01;
    w[19] =   0.9212252223778619E-01;
    w[20] =   0.8689978720108314E-01;
    w[21] =   0.8075589522942023E-01;
    w[22] =   0.7375597473770516E-01;
    w[23] =   0.6597422988218052E-01;
    w[24] =   0.5749315621761905E-01;
    w[25] =   0.4840267283059405E-01;
    w[26] =   0.3879919256962704E-01;
    w[27] =   0.2878470788332330E-01;
    w[28] =   0.1846646831109099E-01;
    w[29] =   0.7968192496166648E-02;
  }
  else if ( n == 31 )
  {
    x[ 0] =  -0.99708748181947707454263838179654;    
    x[ 1] =  -0.98468590966515248400211329970113;    
    x[ 2] =  -0.96250392509294966178905249675943;    
    x[ 3] =  -0.93075699789664816495694576311725;    
    x[ 4] =  -0.88976002994827104337419200908023;    
    x[ 5] =  -0.83992032014626734008690453594388;    
    x[ 6] =  -0.78173314841662494040636002019484;    
    x[ 7] =  -0.71577678458685328390597086536649;    
    x[ 8] =  -0.64270672292426034618441820323250;    
    x[ 9] =  -0.56324916140714926272094492359516;    
    x[10] =  -0.47819378204490248044059403935649;    
    x[11] =  -0.38838590160823294306135146128752;    
    x[12] =  -0.29471806998170161661790389767170;    
    x[13] =  -0.19812119933557062877241299603283;    
    x[14] =  -0.99555312152341520325174790118941E-01;
    x[15] =   0.00000000000000000000000000000000;   
    x[16] =   0.99555312152341520325174790118941E-01;
    x[17] =   0.19812119933557062877241299603283;    
    x[18] =   0.29471806998170161661790389767170;    
    x[19] =   0.38838590160823294306135146128752;    
    x[20] =   0.47819378204490248044059403935649;    
    x[21] =   0.56324916140714926272094492359516;    
    x[22] =   0.64270672292426034618441820323250;    
    x[23] =   0.71577678458685328390597086536649;    
    x[24] =   0.78173314841662494040636002019484;    
    x[25] =   0.83992032014626734008690453594388;    
    x[26] =   0.88976002994827104337419200908023;    
    x[27] =   0.93075699789664816495694576311725;    
    x[28] =   0.96250392509294966178905249675943;    
    x[29] =   0.98468590966515248400211329970113;    
    x[30] =   0.99708748181947707454263838179654;
 
    w[ 0] =   0.74708315792487746093913218970494E-02;
    w[ 1] =   0.17318620790310582463552990782414E-01;
    w[ 2] =   0.27009019184979421800608642617676E-01;
    w[ 3] =   0.36432273912385464024392008749009E-01;
    w[ 4] =   0.45493707527201102902315857856518E-01;
    w[ 5] =   0.54103082424916853711666259085477E-01;
    w[ 6] =   0.62174786561028426910343543686657E-01;
    w[ 7] =   0.69628583235410366167756126255124E-01;
    w[ 8] =   0.76390386598776616426357674901331E-01;
    w[ 9] =   0.82392991761589263903823367431962E-01;
    w[10] =   0.87576740608477876126198069695333E-01;
    w[11] =   0.91890113893641478215362871607150E-01;
    w[12] =   0.95290242912319512807204197487597E-01;
    w[13] =   0.97743335386328725093474010978997E-01;
    w[14] =   0.99225011226672307874875514428615E-01;
    w[15] =   0.99720544793426451427533833734349E-01;
    w[16] =   0.99225011226672307874875514428615E-01;
    w[17] =   0.97743335386328725093474010978997E-01;
    w[18] =   0.95290242912319512807204197487597E-01;
    w[19] =   0.91890113893641478215362871607150E-01;
    w[20] =   0.87576740608477876126198069695333E-01;
    w[21] =   0.82392991761589263903823367431962E-01;
    w[22] =   0.76390386598776616426357674901331E-01;
    w[23] =   0.69628583235410366167756126255124E-01;
    w[24] =   0.62174786561028426910343543686657E-01;
    w[25] =   0.54103082424916853711666259085477E-01;
    w[26] =   0.45493707527201102902315857856518E-01;
    w[27] =   0.36432273912385464024392008749009E-01;
    w[28] =   0.27009019184979421800608642617676E-01;
    w[29] =   0.17318620790310582463552990782414E-01;
    w[30] =   0.74708315792487746093913218970494E-02;
  }
  else if ( n == 32 )
  {
    x[0] =  - 0.997263861849481563544981128665;
    x[1] =  - 0.985611511545268335400175044631;
    x[2] =  - 0.964762255587506430773811928118;
    x[3] =  - 0.934906075937739689170919134835;
    x[4] =  - 0.896321155766052123965307243719;
    x[5] =  - 0.849367613732569970133693004968;
    x[6] =  - 0.794483795967942406963097298970;
    x[7] =  - 0.732182118740289680387426665091;
    x[8] =  - 0.663044266930215200975115168663;
    x[9] =  - 0.587715757240762329040745476402;
    x[10] = - 0.506899908932229390023747474378;
    x[11] = - 0.421351276130635345364119436172;
    x[12] = - 0.331868602282127649779916805730;
    x[13] = - 0.239287362252137074544603209166;
    x[14] = - 0.144471961582796493485186373599;
    x[15] = - 0.483076656877383162348125704405E-01;
    x[16] =   0.483076656877383162348125704405E-01;
    x[17] =   0.144471961582796493485186373599;
    x[18] =   0.239287362252137074544603209166;
    x[19] =   0.331868602282127649779916805730;
    x[20] =   0.421351276130635345364119436172;
    x[21] =   0.506899908932229390023747474378;
    x[22] =   0.587715757240762329040745476402;
    x[23] =   0.663044266930215200975115168663;
    x[24] =   0.732182118740289680387426665091;
    x[25] =   0.794483795967942406963097298970;
    x[26] =   0.849367613732569970133693004968;
    x[27] =   0.896321155766052123965307243719;
    x[28] =   0.934906075937739689170919134835;
    x[29] =   0.964762255587506430773811928118;
    x[30] =   0.985611511545268335400175044631;
    x[31] =   0.997263861849481563544981128665;

    w[0] =  0.701861000947009660040706373885E-02;
    w[1] =  0.162743947309056706051705622064E-01;
    w[2] =  0.253920653092620594557525897892E-01;
    w[3] =  0.342738629130214331026877322524E-01;
    w[4] =  0.428358980222266806568786466061E-01;
    w[5] =  0.509980592623761761961632446895E-01;
    w[6] =  0.586840934785355471452836373002E-01;
    w[7] =  0.658222227763618468376500637069E-01;
    w[8] =  0.723457941088485062253993564785E-01;
    w[9] =  0.781938957870703064717409188283E-01;
    w[10] = 0.833119242269467552221990746043E-01;
    w[11] = 0.876520930044038111427714627518E-01;
    w[12] = 0.911738786957638847128685771116E-01;
    w[13] = 0.938443990808045656391802376681E-01;
    w[14] = 0.956387200792748594190820022041E-01;
    w[15] = 0.965400885147278005667648300636E-01;
    w[16] = 0.965400885147278005667648300636E-01;
    w[17] = 0.956387200792748594190820022041E-01;
    w[18] = 0.938443990808045656391802376681E-01;
    w[19] = 0.911738786957638847128685771116E-01;
    w[20] = 0.876520930044038111427714627518E-01;
    w[21] = 0.833119242269467552221990746043E-01;
    w[22] = 0.781938957870703064717409188283E-01;
    w[23] = 0.723457941088485062253993564785E-01;
    w[24] = 0.658222227763618468376500637069E-01;
    w[25] = 0.586840934785355471452836373002E-01;
    w[26] = 0.509980592623761761961632446895E-01;
    w[27] = 0.428358980222266806568786466061E-01;
    w[28] = 0.342738629130214331026877322524E-01;
    w[29] = 0.253920653092620594557525897892E-01;
    w[30] = 0.162743947309056706051705622064E-01;
    w[31] = 0.701861000947009660040706373885E-02;
  }
  else if ( n == 33 )
  {
    x[ 0] =  -0.9974246942464552;    
    x[ 1] =  -0.9864557262306425;
    x[ 2] =  -0.9668229096899927;
    x[ 3] =  -0.9386943726111684;    
    x[ 4] =  -0.9023167677434336;    
    x[ 5] =  -0.8580096526765041;    
    x[ 6] =  -0.8061623562741665;    
    x[ 7] =  -0.7472304964495622;    
    x[ 8] =  -0.6817319599697428;    
    x[ 9] =  -0.6102423458363790;    
    x[10] =  -0.5333899047863476;    
    x[11] =  -0.4518500172724507;    
    x[12] =  -0.3663392577480734;    
    x[13] =  -0.2776090971524970;    
    x[14] =  -0.1864392988279916;    
    x[15] =  -0.09363106585473338;
    x[16] =   0.000000000000000;
    x[17] =   0.09363106585473338;
    x[18] =   0.1864392988279916;    
    x[19] =   0.2776090971524970;    
    x[20] =   0.3663392577480734;    
    x[21] =   0.4518500172724507;    
    x[22] =   0.5333899047863476;    
    x[23] =   0.6102423458363790;    
    x[24] =   0.6817319599697428;    
    x[25] =   0.7472304964495622;    
    x[26] =   0.8061623562741665;    
    x[27] =   0.8580096526765041;    
    x[28] =   0.9023167677434336;    
    x[29] =   0.9386943726111684;    
    x[30] =   0.9668229096899927;    
    x[31] =   0.9864557262306425;    
    x[32] =   0.9974246942464552;    
 
    w[ 0] =   0.6606227847587558E-02;
    w[ 1] =   0.1532170151293465E-01;
    w[ 2] =   0.2391554810174960E-01;
    w[ 3] =   0.3230035863232891E-01;
    w[ 4] =   0.4040154133166965E-01;
    w[ 5] =   0.4814774281871162E-01;
    w[ 6] =   0.5547084663166357E-01;
    w[ 7] =   0.6230648253031755E-01;
    w[ 8] =   0.6859457281865676E-01;
    w[ 9] =   0.7427985484395420E-01;
    w[10] =   0.7931236479488685E-01;
    w[11] =   0.8364787606703869E-01;
    w[12] =   0.8724828761884425E-01;
    w[13] =   0.9008195866063859E-01;
    w[14] =   0.9212398664331678E-01;
    w[15] =   0.9335642606559612E-01;
    w[16] =   0.9376844616020999E-01;
    w[17] =   0.9335642606559612E-01;
    w[18] =   0.9212398664331678E-01;
    w[19] =   0.9008195866063859E-01;
    w[20] =   0.8724828761884425E-01;
    w[21] =   0.8364787606703869E-01;
    w[22] =   0.7931236479488685E-01;
    w[23] =   0.7427985484395420E-01;
    w[24] =   0.6859457281865676E-01;
    w[25] =   0.6230648253031755E-01;
    w[26] =   0.5547084663166357E-01;
    w[27] =   0.4814774281871162E-01;
    w[28] =   0.4040154133166965E-01;
    w[29] =   0.3230035863232891E-01;
    w[30] =   0.2391554810174960E-01;
    w[31] =   0.1532170151293465E-01;
    w[32] =   0.6606227847587558E-02;
  }
  else if ( n == 64 )
  {
    x[0] =  - 0.999305041735772139456905624346;
    x[1] =  - 0.996340116771955279346924500676;
    x[2] =  - 0.991013371476744320739382383443;
    x[3] =  - 0.983336253884625956931299302157;
    x[4] =  - 0.973326827789910963741853507352;
    x[5] =  - 0.961008799652053718918614121897;
    x[6] =  - 0.946411374858402816062481491347;
    x[7] =  - 0.929569172131939575821490154559;
    x[8] =  - 0.910522137078502805756380668008;
    x[9] =  - 0.889315445995114105853404038273;
    x[10] = - 0.865999398154092819760783385070;
    x[11] = - 0.840629296252580362751691544696;
    x[12] = - 0.813265315122797559741923338086;
    x[13] = - 0.783972358943341407610220525214;
    x[14] = - 0.752819907260531896611863774886;
    x[15] = - 0.719881850171610826848940217832;
    x[16] = - 0.685236313054233242563558371031;
    x[17] = - 0.648965471254657339857761231993;
    x[18] = - 0.611155355172393250248852971019;
    x[19] = - 0.571895646202634034283878116659;
    x[20] = - 0.531279464019894545658013903544;
    x[21] = - 0.489403145707052957478526307022;
    x[22] = - 0.446366017253464087984947714759;
    x[23] = - 0.402270157963991603695766771260;
    x[24] = - 0.357220158337668115950442615046;
    x[25] = - 0.311322871990210956157512698560;
    x[26] = - 0.264687162208767416373964172510;
    x[27] = - 0.217423643740007084149648748989;
    x[28] = - 0.169644420423992818037313629748;
    x[29] = - 0.121462819296120554470376463492;
    x[30] = - 0.729931217877990394495429419403E-01;
    x[31] = - 0.243502926634244325089558428537E-01;
    x[32] =   0.243502926634244325089558428537E-01;
    x[33] =   0.729931217877990394495429419403E-01;
    x[34] =   0.121462819296120554470376463492;
    x[35] =   0.169644420423992818037313629748;
    x[36] =   0.217423643740007084149648748989;
    x[37] =   0.264687162208767416373964172510;
    x[38] =   0.311322871990210956157512698560;
    x[39] =   0.357220158337668115950442615046;
    x[40] =   0.402270157963991603695766771260;
    x[41] =   0.446366017253464087984947714759;
    x[42] =   0.489403145707052957478526307022;
    x[43] =   0.531279464019894545658013903544;
    x[44] =   0.571895646202634034283878116659;
    x[45] =   0.611155355172393250248852971019;
    x[46] =   0.648965471254657339857761231993;
    x[47] =   0.685236313054233242563558371031;
    x[48] =   0.719881850171610826848940217832;
    x[49] =   0.752819907260531896611863774886;
    x[50] =   0.783972358943341407610220525214;
    x[51] =   0.813265315122797559741923338086;
    x[52] =   0.840629296252580362751691544696;
    x[53] =   0.865999398154092819760783385070;
    x[54] =   0.889315445995114105853404038273;
    x[55] =   0.910522137078502805756380668008;
    x[56] =   0.929569172131939575821490154559;
    x[57] =   0.946411374858402816062481491347;
    x[58] =   0.961008799652053718918614121897;
    x[59] =   0.973326827789910963741853507352;
    x[60] =   0.983336253884625956931299302157;
    x[61] =   0.991013371476744320739382383443;
    x[62] =   0.996340116771955279346924500676;
    x[63] =   0.999305041735772139456905624346;

    w[0] =  0.178328072169643294729607914497E-02;
    w[1] =  0.414703326056246763528753572855E-02;
    w[2] =  0.650445796897836285611736039998E-02;
    w[3] =  0.884675982636394772303091465973E-02;
    w[4] =  0.111681394601311288185904930192E-01;
    w[5] =  0.134630478967186425980607666860E-01;
    w[6] =  0.157260304760247193219659952975E-01;
    w[7] =  0.179517157756973430850453020011E-01;
    w[8] =  0.201348231535302093723403167285E-01;
    w[9] =  0.222701738083832541592983303842E-01;
    w[10] = 0.243527025687108733381775504091E-01;
    w[11] = 0.263774697150546586716917926252E-01;
    w[12] = 0.283396726142594832275113052002E-01;
    w[13] = 0.302346570724024788679740598195E-01;
    w[14] = 0.320579283548515535854675043479E-01;
    w[15] = 0.338051618371416093915654821107E-01;
    w[16] = 0.354722132568823838106931467152E-01;
    w[17] = 0.370551285402400460404151018096E-01;
    w[18] = 0.385501531786156291289624969468E-01;
    w[19] = 0.399537411327203413866569261283E-01;
    w[20] = 0.412625632426235286101562974736E-01;
    w[21] = 0.424735151236535890073397679088E-01;
    w[22] = 0.435837245293234533768278609737E-01;
    w[23] = 0.445905581637565630601347100309E-01;
    w[24] = 0.454916279274181444797709969713E-01;
    w[25] = 0.462847965813144172959532492323E-01;
    w[26] = 0.469681828162100173253262857546E-01;
    w[27] = 0.475401657148303086622822069442E-01;
    w[28] = 0.479993885964583077281261798713E-01;
    w[29] = 0.483447622348029571697695271580E-01;
    w[30] = 0.485754674415034269347990667840E-01;
    w[31] = 0.486909570091397203833653907347E-01;
    w[32] = 0.486909570091397203833653907347E-01;
    w[33] = 0.485754674415034269347990667840E-01;
    w[34] = 0.483447622348029571697695271580E-01;
    w[35] = 0.479993885964583077281261798713E-01;
    w[36] = 0.475401657148303086622822069442E-01;
    w[37] = 0.469681828162100173253262857546E-01;
    w[38] = 0.462847965813144172959532492323E-01;
    w[39] = 0.454916279274181444797709969713E-01;
    w[40] = 0.445905581637565630601347100309E-01;
    w[41] = 0.435837245293234533768278609737E-01;
    w[42] = 0.424735151236535890073397679088E-01;
    w[43] = 0.412625632426235286101562974736E-01;
    w[44] = 0.399537411327203413866569261283E-01;
    w[45] = 0.385501531786156291289624969468E-01;
    w[46] = 0.370551285402400460404151018096E-01;
    w[47] = 0.354722132568823838106931467152E-01;
    w[48] = 0.338051618371416093915654821107E-01;
    w[49] = 0.320579283548515535854675043479E-01;
    w[50] = 0.302346570724024788679740598195E-01;
    w[51] = 0.283396726142594832275113052002E-01;
    w[52] = 0.263774697150546586716917926252E-01;
    w[53] = 0.243527025687108733381775504091E-01;
    w[54] = 0.222701738083832541592983303842E-01;
    w[55] = 0.201348231535302093723403167285E-01;
    w[56] = 0.179517157756973430850453020011E-01;
    w[57] = 0.157260304760247193219659952975E-01;
    w[58] = 0.134630478967186425980607666860E-01;
    w[59] = 0.111681394601311288185904930192E-01;
    w[60] = 0.884675982636394772303091465973E-02;
    w[61] = 0.650445796897836285611736039998E-02;
    w[62] = 0.414703326056246763528753572855E-02;
    w[63] = 0.178328072169643294729607914497E-02;
  }
  else if ( n == 65 )
  {
    x[ 0] =  -0.9993260970754129;    
    x[ 1] =  -0.9964509480618492;    
    x[ 2] =  -0.9912852761768016;    
    x[ 3] =  -0.9838398121870350;    
    x[ 4] =  -0.9741315398335512;    
    x[ 5] =  -0.9621827547180553;    
    x[ 6] =  -0.9480209281684076;    
    x[ 7] =  -0.9316786282287494;    
    x[ 8] =  -0.9131934405428462;    
    x[ 9] =  -0.8926078805047389;    
    x[10] =  -0.8699692949264071;    
    x[11] =  -0.8453297528999303;    
    x[12] =  -0.8187459259226514;    
    x[13] =  -0.7902789574921218;    
    x[14] =  -0.7599943224419998;    
    x[15] =  -0.7279616763294247;    
    x[16] =  -0.6942546952139916;    
    x[17] =  -0.6589509061936252;    
    x[18] =  -0.6221315090854003;    
    x[19] =  -0.5838811896604873;    
    x[20] =  -0.5442879248622271;    
    x[21] =  -0.5034427804550069;    
    x[22] =  -0.4614397015691450;    
    x[23] =  -0.4183752966234090;    
    x[24] =  -0.3743486151220660;    
    x[25] =  -0.3294609198374864;    
    x[26] =  -0.2838154539022487;    
    x[27] =  -0.2375172033464168;    
    x[28] =  -0.1906726556261428;    
    x[29] =  -0.1433895546989752;    
    x[30] =  -0.9577665320919751E-01;
    x[31] =  -0.4794346235317186E-01;
    x[32] =    0.000000000000000;    
    x[33] =   0.4794346235317186E-01;
    x[34] =   0.9577665320919751E-01;
    x[35] =   0.1433895546989752;    
    x[36] =   0.1906726556261428;    
    x[37] =   0.2375172033464168;    
    x[38] =   0.2838154539022487;    
    x[39] =   0.3294609198374864;    
    x[40] =   0.3743486151220660;    
    x[41] =   0.4183752966234090;    
    x[42] =   0.4614397015691450;    
    x[43] =   0.5034427804550069;    
    x[44] =   0.5442879248622271;    
    x[45] =   0.5838811896604873;    
    x[46] =   0.6221315090854003;    
    x[47] =   0.6589509061936252;    
    x[48] =   0.6942546952139916;    
    x[49] =   0.7279616763294247;    
    x[50] =   0.7599943224419998;    
    x[51] =   0.7902789574921218;    
    x[52] =   0.8187459259226514;    
    x[53] =   0.8453297528999303;    
    x[54] =   0.8699692949264071;    
    x[55] =   0.8926078805047389;    
    x[56] =   0.9131934405428462;    
    x[57] =   0.9316786282287494;    
    x[58] =   0.9480209281684076;    
    x[59] =   0.9621827547180553;    
    x[60] =   0.9741315398335512;    
    x[61] =   0.9838398121870350;    
    x[62] =   0.9912852761768016;    
    x[63] =   0.9964509480618492;    
    x[64] =   0.9993260970754129;    
 
    w[ 0] =   0.1729258251300218E-02;
    w[ 1] =   0.4021524172003703E-02;
    w[ 2] =   0.6307942578971821E-02;
    w[ 3] =   0.8580148266881443E-02;
    w[ 4] =   0.1083267878959798E-01;
    w[ 5] =   0.1306031163999490E-01;
    w[ 6] =   0.1525791214644825E-01;
    w[ 7] =   0.1742042199767025E-01;
    w[ 8] =   0.1954286583675005E-01;
    w[ 9] =   0.2162036128493408E-01;
    w[10] =   0.2364812969128723E-01;
    w[11] =   0.2562150693803776E-01;
    w[12] =   0.2753595408845034E-01;
    w[13] =   0.2938706778931066E-01;
    w[14] =   0.3117059038018911E-01;
    w[15] =   0.3288241967636860E-01;
    w[16] =   0.3451861839854901E-01;
    w[17] =   0.3607542322556527E-01;
    w[18] =   0.3754925344825770E-01;
    w[19] =   0.3893671920405121E-01;
    w[20] =   0.4023462927300549E-01;
    w[21] =   0.4143999841724028E-01;
    w[22] =   0.4255005424675579E-01;
    w[23] =   0.4356224359580051E-01;
    w[24] =   0.4447423839508296E-01;
    w[25] =   0.4528394102630023E-01;
    w[26] =   0.4598948914665173E-01;
    w[27] =   0.4658925997223349E-01;
    w[28] =   0.4708187401045461E-01;
    w[29] =   0.4746619823288551E-01;
    w[30] =   0.4774134868124067E-01;
    w[31] =   0.4790669250049590E-01;
    w[32] =   0.4796184939446662E-01;
    w[33] =   0.4790669250049590E-01;
    w[34] =   0.4774134868124067E-01;
    w[35] =   0.4746619823288551E-01;
    w[36] =   0.4708187401045461E-01;
    w[37] =   0.4658925997223349E-01;
    w[38] =   0.4598948914665173E-01;
    w[39] =   0.4528394102630023E-01;
    w[40] =   0.4447423839508296E-01;
    w[41] =   0.4356224359580051E-01;
    w[42] =   0.4255005424675579E-01;
    w[43] =   0.4143999841724028E-01;
    w[44] =   0.4023462927300549E-01;
    w[45] =   0.3893671920405121E-01;
    w[46] =   0.3754925344825770E-01;
    w[47] =   0.3607542322556527E-01;
    w[48] =   0.3451861839854901E-01;
    w[49] =   0.3288241967636860E-01;
    w[50] =   0.3117059038018911E-01;
    w[51] =   0.2938706778931066E-01;
    w[52] =   0.2753595408845034E-01;
    w[53] =   0.2562150693803776E-01;
    w[54] =   0.2364812969128723E-01;
    w[55] =   0.2162036128493408E-01;
    w[56] =   0.1954286583675005E-01;
    w[57] =   0.1742042199767025E-01;
    w[58] =   0.1525791214644825E-01;
    w[59] =   0.1306031163999490E-01;
    w[60] =   0.1083267878959798E-01;
    w[61] =   0.8580148266881443E-02;
    w[62] =   0.6307942578971821E-02;
    w[63] =   0.4021524172003703E-02;
    w[64] =   0.1729258251300218E-02;
  }
  else if ( n == 127 ) 
  {
    x[  0] =  -0.99982213041530614629963254927125E+00;
    x[  1] =  -0.99906293435531189513828920479421E+00;    
    x[  2] =  -0.99769756618980462107441703193392E+00;    
    x[  3] =  -0.99572655135202722663543337085008E+00;    
    x[  4] =  -0.99315104925451714736113079489080E+00;  
    x[  5] =  -0.98997261459148415760778669967548E+00;   
    x[  6] =  -0.98619317401693166671043833175407E+00;    
    x[  7] =  -0.98181502080381411003346312451200E+00;    
    x[  8] =  -0.97684081234307032681744391886221E+00;    
    x[  9] =  -0.97127356816152919228894689830512E+00;    
    x[ 10] =  -0.96511666794529212109082507703391E+00;    
    x[ 11] =  -0.95837384942523877114910286998060E+00;    
    x[ 12] =  -0.95104920607788031054790764659636E+00;   
    x[ 13] =  -0.94314718462481482734544963026201E+00;    
    x[ 14] =  -0.93467258232473796857363487794906E+00;    
    x[ 15] =  -0.92563054405623384912746466814259E+00;    
    x[ 16] =  -0.91602655919146580931308861741716E+00;   
    x[ 17] =  -0.90586645826182138280246131760282E+00;    
    x[ 18] =  -0.89515640941708370896904382642451E+00;   
    x[ 19] =  -0.88390291468002656994525794802849E+00;    
    x[ 20] =  -0.87211280599856071141963753428864E+00;    
    x[ 21] =  -0.85979324109774080981203134414483E+00;   
    x[ 22] =  -0.84695169913409759845333931085437E+00;    
    x[ 23] =  -0.83359597615489951437955716480123E+00;    
    x[ 24] =  -0.81973418036507867415511910167470E+00;   
    x[ 25] =  -0.80537472720468021466656079404644E+00;   
    x[ 26] =  -0.79052633423981379994544995252740E+00;   
    x[ 27] =  -0.77519801587020238244496276354566E+00;  
    x[ 28] =  -0.75939907785653667155666366659810E+00;   
    x[ 29] =  -0.74313911167095451292056688997595E+00;   
    x[ 30] =  -0.72642798867407268553569290153270E+00;    
    x[ 31] =  -0.70927585412210456099944463906757E+00;   
    x[ 32] =  -0.69169312100770067015644143286666E+00; 
    x[ 33] =  -0.67369046373825048534668253831602E+00;
    x[ 34] =  -0.65527881165548263027676505156852E+00;
    x[ 35] =  -0.63646934240029724134760815684175E+00;
    x[ 36] =  -0.61727347512685828385763916340822E+00; 
    x[ 37] =  -0.59770286357006522938441201887478E+00; 
    x[ 38] =  -0.57776938897061258000325165713764E+00; 
    x[ 39] =  -0.55748515286193223292186190687872E+00; 
    x[ 40] =  -0.53686246972339756745816636353452E+00;
    x[ 41] =  -0.51591385950424935727727729906662E+00; 
    x[ 42] =  -0.49465204002278211739494017368636E+00;
    x[ 43] =  -0.47308991924540524164509989939699E+00;
    x[ 44] =  -0.45124058745026622733189858020729E+00;
    x[ 45] =  -0.42911730928019337626254405355418E+00;
    x[ 46] =  -0.40673351568978256340867288124339E+00;
    x[ 47] =  -0.38410279579151693577907781452239E+00;
    x[ 48] =  -0.36123888860586970607092484346723E+00;
    x[ 49] =  -0.33815567472039850137600027657095E+00;
    x[ 50] =  -0.31486716786289498148601475374890E+00; 
    x[ 51] =  -0.29138750639370562079451875284568E+00; 
    x[ 52] =  -0.26773094472238862088834352027938E+00;
    x[ 53] =  -0.24391184465391785797071324453138E+00;
    x[ 54] =  -0.21994466666968754245452337866940E+00;
    x[ 55] =  -0.19584396114861085150428162519610E+00;
    x[ 56] =  -0.17162435953364216500834492248954E+00; 
    x[ 57] =  -0.14730056544908566938932929319807E+00;
    x[ 58] =  -0.12288734577408297172603365288567E+00;
    x[ 59] =  -0.98399521677698970751091751509101E-01;
    x[ 60] =  -0.73851959621048545273440409360569E-01;
    x[ 61] =  -0.49259562331926630315379321821927E-01;
    x[ 62] =  -0.24637259757420944614897071846088E-01;
    x[ 63] =   0.00000000000000000000000000000000E+00;
    x[ 64] =   0.24637259757420944614897071846088E-01;
    x[ 65] =   0.49259562331926630315379321821927E-01;
    x[ 66] =   0.73851959621048545273440409360569E-01;
    x[ 67] =   0.98399521677698970751091751509101E-01;
    x[ 68] =   0.12288734577408297172603365288567E+00;
    x[ 69] =   0.14730056544908566938932929319807E+00;
    x[ 70] =   0.17162435953364216500834492248954E+00;
    x[ 71] =   0.19584396114861085150428162519610E+00;
    x[ 72] =   0.21994466666968754245452337866940E+00;    
    x[ 73] =   0.24391184465391785797071324453138E+00;   
    x[ 74] =   0.26773094472238862088834352027938E+00;   
    x[ 75] =   0.29138750639370562079451875284568E+00;   
    x[ 76] =   0.31486716786289498148601475374890E+00;    
    x[ 77] =   0.33815567472039850137600027657095E+00;   
    x[ 78] =   0.36123888860586970607092484346723E+00;    
    x[ 79] =   0.38410279579151693577907781452239E+00;    
    x[ 80] =   0.40673351568978256340867288124339E+00;  
    x[ 81] =   0.42911730928019337626254405355418E+00;    
    x[ 82] =   0.45124058745026622733189858020729E+00;   
    x[ 83] =   0.47308991924540524164509989939699E+00;   
    x[ 84] =   0.49465204002278211739494017368636E+00; 
    x[ 85] =   0.51591385950424935727727729906662E+00; 
    x[ 86] =   0.53686246972339756745816636353452E+00; 
    x[ 87] =   0.55748515286193223292186190687872E+00;   
    x[ 88] =   0.57776938897061258000325165713764E+00;  
    x[ 89] =   0.59770286357006522938441201887478E+00;  
    x[ 90] =   0.61727347512685828385763916340822E+00;  
    x[ 91] =   0.63646934240029724134760815684175E+00;    
    x[ 92] =   0.65527881165548263027676505156852E+00;  
    x[ 93] =   0.67369046373825048534668253831602E+00;   
    x[ 94] =   0.69169312100770067015644143286666E+00;   
    x[ 95] =   0.70927585412210456099944463906757E+00;   
    x[ 96] =   0.72642798867407268553569290153270E+00;   
    x[ 97] =   0.74313911167095451292056688997595E+00;    
    x[ 98] =   0.75939907785653667155666366659810E+00;   
    x[ 99] =   0.77519801587020238244496276354566E+00;    
    x[100] =   0.79052633423981379994544995252740E+00;   
    x[101] =   0.80537472720468021466656079404644E+00;   
    x[102] =   0.81973418036507867415511910167470E+00;  
    x[103] =   0.83359597615489951437955716480123E+00;   
    x[104] =   0.84695169913409759845333931085437E+00;   
    x[105] =   0.85979324109774080981203134414483E+00; 
    x[106] =   0.87211280599856071141963753428864E+00;  
    x[107] =   0.88390291468002656994525794802849E+00;   
    x[108] =   0.89515640941708370896904382642451E+00;    
    x[109] =   0.90586645826182138280246131760282E+00;   
    x[110] =   0.91602655919146580931308861741716E+00;  
    x[111] =   0.92563054405623384912746466814259E+00; 
    x[112] =   0.93467258232473796857363487794906E+00; 
    x[113] =   0.94314718462481482734544963026201E+00;  
    x[114] =   0.95104920607788031054790764659636E+00; 
    x[115] =   0.95837384942523877114910286998060E+00; 
    x[116] =   0.96511666794529212109082507703391E+00;
    x[117] =   0.97127356816152919228894689830512E+00; 
    x[118] =   0.97684081234307032681744391886221E+00; 
    x[119] =   0.98181502080381411003346312451200E+00;  
    x[120] =   0.98619317401693166671043833175407E+00;
    x[121] =   0.98997261459148415760778669967548E+00;
    x[122] =   0.99315104925451714736113079489080E+00; 
    x[123] =   0.99572655135202722663543337085008E+00; 
    x[124] =   0.99769756618980462107441703193392E+00; 
    x[125] =   0.99906293435531189513828920479421E+00;
    x[126] =   0.99982213041530614629963254927125E+00; 

    w[  0] =   0.45645726109586654495731936146574E-03;
    w[  1] =   0.10622766869538486959954760554099E-02;
    w[  2] =   0.16683488125171936761028811985672E-02;
    w[  3] =   0.22734860707492547802810838362671E-02;
    w[  4] =   0.28772587656289004082883197417581E-02;
    w[  5] =   0.34792893810051465908910894094105E-02;
    w[  6] =   0.40792095178254605327114733456293E-02;
    w[  7] =   0.46766539777779034772638165662478E-02;
    w[  8] =   0.52712596565634400891303815906251E-02;
    w[  9] =   0.58626653903523901033648343751367E-02;
    w[ 10] =   0.64505120486899171845442463868748E-02;
    w[ 11] =   0.70344427036681608755685893032552E-02;
    w[ 12] =   0.76141028256526859356393930849227E-02;
    w[ 13] =   0.81891404887415730817235884718726E-02;
    w[ 14] =   0.87592065795403145773316804234385E-02;
    w[ 15] =   0.93239550065309714787536985834029E-02;
    w[ 16] =   0.98830429087554914716648010899606E-02;
    w[ 17] =   0.10436130863141005225673171997668E-01;
    w[ 18] =   0.10982883090068975788799657376065E-01;
    w[ 19] =   0.11522967656921087154811609734510E-01;
    w[ 20] =   0.12056056679400848183529562144697E-01;
    w[ 21] =   0.12581826520465013101514365424172E-01;
    w[ 22] =   0.13099957986718627426172681912499E-01;
    w[ 23] =   0.13610136522139249906034237533759E-01;
    w[ 24] =   0.14112052399003395774044161633613E-01;
    w[ 25] =   0.14605400905893418351737288078952E-01;
    w[ 26] =   0.15089882532666922992635733981431E-01;
    w[ 27] =   0.15565203152273955098532590262975E-01;
    w[ 28] =   0.16031074199309941802254151842763E-01;
    w[ 29] =   0.16487212845194879399346060358146E-01;
    w[ 30] =   0.16933342169871654545878815295200E-01;
    w[ 31] =   0.17369191329918731922164721250350E-01;
    w[ 32] =   0.17794495722974774231027912900351E-01;
    w[ 33] =   0.18208997148375106468721469154479E-01;
    w[ 34] =   0.18612443963902310429440419898958E-01;
    w[ 35] =   0.19004591238555646611148901044533E-01;
    w[ 36] =   0.19385200901246454628112623489471E-01;
    w[ 37] =   0.19754041885329183081815217323169E-01;
    w[ 38] =   0.20110890268880247225644623956287E-01;
    w[ 39] =   0.20455529410639508279497065713301E-01;
    w[ 40] =   0.20787750081531811812652137291250E-01;
    w[ 41] =   0.21107350591688713643523847921658E-01;
    w[ 42] =   0.21414136912893259295449693233545E-01;
    w[ 43] =   0.21707922796373466052301324695331E-01;
    w[ 44] =   0.21988529885872983756478409758807E-01;
    w[ 45] =   0.22255787825930280235631416460158E-01;
    w[ 46] =   0.22509534365300608085694429903050E-01;
    w[ 47] =   0.22749615455457959852242553240982E-01;
    w[ 48] =   0.22975885344117206754377437838947E-01;
    w[ 49] =   0.23188206663719640249922582981729E-01;
    w[ 50] =   0.23386450514828194170722043496950E-01;
    w[ 51] =   0.23570496544381716050033676844306E-01;
    w[ 52] =   0.23740233018760777777714726703424E-01;
    w[ 53] =   0.23895556891620665983864481754172E-01;
    w[ 54] =   0.24036373866450369675132086026456E-01;
    w[ 55] =   0.24162598453819584716522917710986E-01;
    w[ 56] =   0.24274154023278979833195063936748E-01;
    w[ 57] =   0.24370972849882214952813561907241E-01;
    w[ 58] =   0.24452996155301467956140198471529E-01;
    w[ 59] =   0.24520174143511508275183033290175E-01;
    w[ 60] =   0.24572466031020653286354137335186E-01;
    w[ 61] =   0.24609840071630254092545634003360E-01;
    w[ 62] =   0.24632273575707679066033370218017E-01;
    w[ 63] =   0.24639752923961094419579417477503E-01;
    w[ 64] =   0.24632273575707679066033370218017E-01;
    w[ 65] =   0.24609840071630254092545634003360E-01;
    w[ 66] =   0.24572466031020653286354137335186E-01;
    w[ 67] =   0.24520174143511508275183033290175E-01;
    w[ 68] =   0.24452996155301467956140198471529E-01;
    w[ 69] =   0.24370972849882214952813561907241E-01;
    w[ 70] =   0.24274154023278979833195063936748E-01;
    w[ 71] =   0.24162598453819584716522917710986E-01;
    w[ 72] =   0.24036373866450369675132086026456E-01;
    w[ 73] =   0.23895556891620665983864481754172E-01;
    w[ 74] =   0.23740233018760777777714726703424E-01;
    w[ 75] =   0.23570496544381716050033676844306E-01;
    w[ 76] =   0.23386450514828194170722043496950E-01;
    w[ 77] =   0.23188206663719640249922582981729E-01;
    w[ 78] =   0.22975885344117206754377437838947E-01;
    w[ 79] =   0.22749615455457959852242553240982E-01;
    w[ 80] =   0.22509534365300608085694429903050E-01;
    w[ 81] =   0.22255787825930280235631416460158E-01;
    w[ 82] =   0.21988529885872983756478409758807E-01;
    w[ 83] =   0.21707922796373466052301324695331E-01;
    w[ 84] =   0.21414136912893259295449693233545E-01;
    w[ 85] =   0.21107350591688713643523847921658E-01;
    w[ 86] =   0.20787750081531811812652137291250E-01;
    w[ 87] =   0.20455529410639508279497065713301E-01;
    w[ 88] =   0.20110890268880247225644623956287E-01;
    w[ 89] =   0.19754041885329183081815217323169E-01;
    w[ 90] =   0.19385200901246454628112623489471E-01;
    w[ 91] =   0.19004591238555646611148901044533E-01;
    w[ 92] =   0.18612443963902310429440419898958E-01;
    w[ 93] =   0.18208997148375106468721469154479E-01;
    w[ 94] =   0.17794495722974774231027912900351E-01;
    w[ 95] =   0.17369191329918731922164721250350E-01;
    w[ 96] =   0.16933342169871654545878815295200E-01;
    w[ 97] =   0.16487212845194879399346060358146E-01;
    w[ 98] =   0.16031074199309941802254151842763E-01;
    w[ 99] =   0.15565203152273955098532590262975E-01;
    w[100] =   0.15089882532666922992635733981431E-01;
    w[101] =   0.14605400905893418351737288078952E-01;
    w[102] =   0.14112052399003395774044161633613E-01;
    w[103] =   0.13610136522139249906034237533759E-01;
    w[104] =   0.13099957986718627426172681912499E-01;
    w[105] =   0.12581826520465013101514365424172E-01;
    w[106] =   0.12056056679400848183529562144697E-01;
    w[107] =   0.11522967656921087154811609734510E-01;
    w[108] =   0.10982883090068975788799657376065E-01;
    w[109] =   0.10436130863141005225673171997668E-01;
    w[110] =   0.98830429087554914716648010899606E-02;
    w[111] =   0.93239550065309714787536985834029E-02;
    w[112] =   0.87592065795403145773316804234385E-02;
    w[113] =   0.81891404887415730817235884718726E-02;
    w[114] =   0.76141028256526859356393930849227E-02;
    w[115] =   0.70344427036681608755685893032552E-02;
    w[116] =   0.64505120486899171845442463868748E-02;
    w[117] =   0.58626653903523901033648343751367E-02;
    w[118] =   0.52712596565634400891303815906251E-02;
    w[119] =   0.46766539777779034772638165662478E-02;
    w[120] =   0.40792095178254605327114733456293E-02;
    w[121] =   0.34792893810051465908910894094105E-02;
    w[122] =   0.28772587656289004082883197417581E-02;
    w[123] =   0.22734860707492547802810838362671E-02;
    w[124] =   0.16683488125171936761028811985672E-02;
    w[125] =   0.10622766869538486959954760554099E-02;
    w[126] =   0.45645726109586654495731936146574E-03;
  }
  else if ( n == 255 )
  {
    x[ 0] =      -0.9999557053175637;
    x[ 1] =      -0.9997666213120006;
    x[ 2] =        -0.99942647468017;
    x[ 3] =      -0.9989352412846546;
    x[ 4] =      -0.9982929861369679;
    x[ 5] =      -0.9974998041266158;
    x[ 6] =      -0.9965558144351986;
    x[ 7] =      -0.9954611594800263;
    x[ 8] =      -0.9942160046166302;
    x[ 9] =      -0.9928205380219891;
    x[10] =      -0.9912749706303856;
    x[11] =      -0.9895795360859201;
    x[12] =      -0.9877344906997324;
    x[13] =      -0.9857401134074193;
    x[14] =      -0.9835967057247763;
    x[15] =      -0.9813045917010171;
    x[16] =      -0.9788641178690681;
    x[17] =       -0.976275653192736;
    x[18] =      -0.9735395890106436;
    x[19] =      -0.9706563389768804;
    x[20] =      -0.9676263389983388;
    x[21] =      -0.9644500471687263;
    x[22] =      -0.9611279436992478;
    x[23] =       -0.957660530845962;
    x[24] =      -0.9540483328338163;
    x[25] =      -0.9502918957773683;
    x[26] =      -0.9463917875982043;
    x[27] =      -0.9423485979390644;
    x[28] =      -0.9381629380746873;
    x[29] =      -0.9338354408193861;
    x[30] =      -0.9293667604313699;
    x[31] =      -0.9247575725138244;
    x[32] =      -0.9200085739127664;
    x[33] =       -0.915120482611687;
    x[34] =      -0.9100940376230008;
    x[35] =       -0.904929998876315;
    x[36] =      -0.8996291471035368;
    x[37] =      -0.8941922837208367;
    x[38] =      -0.8886202307074841;
    x[39] =      -0.8829138304815741;
    x[40] =      -0.8770739457726654;
    x[41] =      -0.8711014594913465;
    x[42] =      -0.8649972745957512;
    x[43] =       -0.858762313955043;
    x[44] =      -0.8523975202098902;
    x[45] =      -0.8459038556299511;
    x[46] =       -0.839282301968391;
    x[47] =      -0.8325338603134556;
    x[48] =      -0.8256595509371186;
    x[49] =      -0.8186604131408319;
    x[50] =      -0.8115375050983958;
    x[51] =      -0.8042919036959787;
    x[52] =      -0.7969247043693057;
    x[53] =      -0.7894370209380444;
    x[54] =      -0.7818299854374094;
    x[55] =      -0.7741047479470157;
    x[56] =      -0.7662624764170006;
    x[57] =      -0.7583043564914468;
    x[58] =      -0.7502315913291283;
    x[59] =      -0.7420454014216102;
    x[60] =      -0.7337470244087263;
    x[61] =      -0.7253377148914649;
    x[62] =      -0.7168187442422908;
    x[63] =      -0.7081914004129306;
    x[64] =      -0.6994569877396524;
    x[65] =      -0.6906168267460676;
    x[66] =      -0.6816722539434864;
    x[67] =      -0.6726246216288551;
    x[68] =       -0.663475297680307;
    x[69] =      -0.6542256653503588;
    x[70] =      -0.6448771230567811;
    x[71] =      -0.6354310841711771;
    x[72] =      -0.6258889768052999;
    x[73] =      -0.6162522435951415;
    x[74] =      -0.6065223414828266;
    x[75] =      -0.5967007414963417;
    x[76] =      -0.5867889285271373;
    x[77] =      -0.5767884011056313;
    x[78] =      -0.5667006711746527;
    x[79] =      -0.5565272638608558;
    x[80] =      -0.5462697172441424;
    x[81] =      -0.5359295821251249;
    x[82] =      -0.5255084217906666;
    x[83] =      -0.5150078117775342;
    x[84] =      -0.5044293396341982;
    x[85] =       -0.493774604680817;
    x[86] =       -0.483045217767442;
    x[87] =      -0.4722428010304787;
    x[88] =      -0.4613689876474424;
    x[89] =      -0.4504254215900437;
    x[90] =      -0.4394137573756426;
    x[91] =      -0.4283356598171081;
    x[92] =      -0.4171928037711214;
    x[93] =      -0.4059868738849605;
    x[94] =      -0.3947195643418044;
    x[95] =      -0.3833925786045958;
    x[96] =      -0.3720076291585012;
    x[97] =      -0.3605664372520062;
    x[98] =      -0.3490707326366864;
    x[99] =      -0.3375222533056927;
    x[100] =      -0.3259227452309905;
    x[101] =      -0.3142739620993925;
    x[102] =      -0.3025776650474256;
    x[103] =      -0.2908356223950708;
    x[104] =      -0.2790496093784178;
    x[105] =      -0.2672214078812731;
    x[106] =      -0.2553528061657641;
    x[107] =       -0.243445598601978;
    x[108] =      -0.2315015853966777;
    x[109] =      -0.2195225723211354;
    x[110] =      -0.2075103704381242;
    x[111] =      -0.1954667958281108;
    x[112] =      -0.1833936693146885;
    x[113] =      -0.1712928161892939;
    x[114] =      -0.1591660659352477;
    x[115] =       -0.147015251951162;
    x[116] =      -0.1348422112737553;
    x[117] =      -0.1226487843001178;
    x[118] =      -0.1104368145094688;
    x[119] =     -0.09820814818444755;
    x[120] =     -0.08596463413198061;
    x[121] =     -0.07370812340376778;
    x[122] =     -0.06144046901642827;
    x[123] =     -0.04916352567134998;
    x[124] =     -0.03687914947428402;
    x[125] =     -0.02458919765472701;
    x[126] =     -0.01229552828513332;
    x[127] =                        0;
    x[128] =      0.01229552828513332;
    x[129] =      0.02458919765472701;
    x[130] =      0.03687914947428402;
    x[131] =      0.04916352567134998;
    x[132] =      0.06144046901642827;
    x[133] =      0.07370812340376778;
    x[134] =      0.08596463413198061;
    x[135] =      0.09820814818444755;
    x[136] =       0.1104368145094688;
    x[137] =       0.1226487843001178;
    x[138] =       0.1348422112737553;
    x[139] =        0.147015251951162;
    x[140] =       0.1591660659352477;
    x[141] =       0.1712928161892939;
    x[142] =       0.1833936693146885;
    x[143] =       0.1954667958281108;
    x[144] =       0.2075103704381242;
    x[145] =       0.2195225723211354;
    x[146] =       0.2315015853966777;
    x[147] =        0.243445598601978;
    x[148] =       0.2553528061657641;
    x[149] =       0.2672214078812731;
    x[150] =       0.2790496093784178;
    x[151] =       0.2908356223950708;
    x[152] =       0.3025776650474256;
    x[153] =       0.3142739620993925;
    x[154] =       0.3259227452309905;
    x[155] =       0.3375222533056927;
    x[156] =       0.3490707326366864;
    x[157] =       0.3605664372520062;
    x[158] =       0.3720076291585012;
    x[159] =       0.3833925786045958;
    x[160] =       0.3947195643418044;
    x[161] =       0.4059868738849605;
    x[162] =       0.4171928037711214;
    x[163] =       0.4283356598171081;
    x[164] =       0.4394137573756426;
    x[165] =       0.4504254215900437;
    x[166] =       0.4613689876474424;
    x[167] =       0.4722428010304787;
    x[168] =        0.483045217767442;
    x[169] =        0.493774604680817;
    x[170] =       0.5044293396341982;
    x[171] =       0.5150078117775342;
    x[172] =       0.5255084217906666;
    x[173] =       0.5359295821251249;
    x[174] =       0.5462697172441424;
    x[175] =       0.5565272638608558;
    x[176] =       0.5667006711746527;
    x[177] =       0.5767884011056313;
    x[178] =       0.5867889285271373;
    x[179] =       0.5967007414963417;
    x[180] =       0.6065223414828266;
    x[181] =       0.6162522435951415;
    x[182] =       0.6258889768052999;
    x[183] =       0.6354310841711771;
    x[184] =       0.6448771230567811;
    x[185] =       0.6542256653503588;
    x[186] =        0.663475297680307;
    x[187] =       0.6726246216288551;
    x[188] =       0.6816722539434864;
    x[189] =       0.6906168267460676;
    x[190] =       0.6994569877396524;
    x[191] =       0.7081914004129306;
    x[192] =       0.7168187442422908;
    x[193] =       0.7253377148914649;
    x[194] =       0.7337470244087263;
    x[195] =       0.7420454014216102;
    x[196] =       0.7502315913291283;
    x[197] =       0.7583043564914468;
    x[198] =       0.7662624764170006;
    x[199] =       0.7741047479470157;
    x[200] =       0.7818299854374094;
    x[201] =       0.7894370209380444;
    x[202] =       0.7969247043693057;
    x[203] =       0.8042919036959787;
    x[204] =       0.8115375050983958;
    x[205] =       0.8186604131408319;
    x[206] =       0.8256595509371186;
    x[207] =       0.8325338603134556;
    x[208] =        0.839282301968391;
    x[209] =       0.8459038556299511;
    x[210] =       0.8523975202098902;
    x[211] =        0.858762313955043;
    x[212] =       0.8649972745957512;
    x[213] =       0.8711014594913465;
    x[214] =       0.8770739457726654;
    x[215] =       0.8829138304815741;
    x[216] =       0.8886202307074841;
    x[217] =       0.8941922837208367;
    x[218] =       0.8996291471035368;
    x[219] =        0.904929998876315;
    x[220] =       0.9100940376230008;
    x[221] =        0.915120482611687;
    x[222] =       0.9200085739127664;
    x[223] =       0.9247575725138244;
    x[224] =       0.9293667604313699;
    x[225] =       0.9338354408193861;
    x[226] =       0.9381629380746873;
    x[227] =       0.9423485979390644;
    x[228] =       0.9463917875982043;
    x[229] =       0.9502918957773683;
    x[230] =       0.9540483328338163;
    x[231] =        0.957660530845962;
    x[232] =       0.9611279436992478;
    x[233] =       0.9644500471687263;
    x[234] =       0.9676263389983388;
    x[235] =       0.9706563389768804;
    x[236] =       0.9735395890106436;
    x[237] =        0.976275653192736;
    x[238] =       0.9788641178690681;
    x[239] =       0.9813045917010171;
    x[240] =       0.9835967057247763;
    x[241] =       0.9857401134074193;
    x[242] =       0.9877344906997324;
    x[243] =       0.9895795360859201;
    x[244] =       0.9912749706303856;
    x[245] =       0.9928205380219891;
    x[246] =       0.9942160046166302;
    x[247] =       0.9954611594800263;
    x[248] =       0.9965558144351986;
    x[249] =       0.9974998041266158;
    x[250] =       0.9982929861369679;
    x[251] =       0.9989352412846546;
    x[252] =         0.99942647468017;
    x[253] =       0.9997666213120006;
    x[254] =       0.9999557053175637;

    w[ 0] =    0.0001136736199914808;
    w[ 1] =    0.0002645938711908564;
    w[ 2] =    0.0004156976252681932;
    w[ 3] =    0.0005667579456482639;
    w[ 4] =    0.0007177364780061286;
    w[ 5] =    0.0008686076661194581;
    w[ 6] =     0.001019347976427318;
    w[ 7] =       0.0011699343729388;
    w[ 8] =     0.001320343990022177;
    w[ 9] =     0.001470554042778403;
    w[10] =     0.001620541799041545;
    w[11] =     0.001770284570660304;
    w[12] =     0.001919759711713187;
    w[13] =     0.002068944619501569;
    w[14] =     0.002217816736754017;
    w[15] =     0.002366353554396287;
    w[16] =      0.00251453261459971;
    w[17] =     0.002662331513971696;
    w[18] =      0.00280972790682046;
    w[19] =     0.002956699508457498;
    w[20] =     0.003103224098519095;
    w[21] =     0.003249279524294296;
    w[22] =     0.003394843704053401;
    w[23] =     0.003539894630372244;
    w[24] =     0.003684410373449933;
    w[25] =     0.003828369084417135;
    w[26] =     0.003971748998634907;
    w[27] =     0.004114528438981242;
    w[28] =     0.004256685819126112;
    w[29] =     0.004398199646792759;
    w[30] =      0.00453904852700618;
    w[31] =     0.004679211165326077;
    w[32] =     0.004818666371065699;
    w[33] =      0.00495739306049505;
    w[34] =     0.005095370260027839;
    w[35] =     0.005232577109391968;
    w[36] =     0.005368992864783177;
    w[37] =     0.005504596902000804;
    w[38] =     0.005639368719565862;
    w[39] =     0.005773287941820301;
    w[40] =     0.005906334322007422;
    w[41] =     0.006038487745332765;
    w[42] =     0.006169728232005295;
    w[43] =     0.006300035940257733;
    w[44] =     0.006429391169346602;
    w[45] =     0.006557774362530328;
    w[46] =     0.006685166110026254;
    w[47] =     0.006811547151944815;
    w[48] =     0.006936898381201466;
    w[49] =     0.007061200846405536;
    w[50] =     0.007184435754724984;
    w[51] =     0.007306584474728122;
    w[52] =     0.007427628539199977;
    w[53] =     0.007547549647934514;
    w[54] =     0.007666329670501377;
    w[55] =     0.007783950648986801;
    w[56] =     0.007900394800708624;
    w[57] =     0.008015644520904983;
    w[58] =     0.008129682385395602;
    w[59] =     0.008242491153216323;
    w[60] =     0.008354053769225508;
    w[61] =     0.008464353366682819;
    w[62] =     0.008573373269798925;
    w[63] =     0.008681096996256795;
    w[64] =     0.008787508259703609;
    w[65] =     0.008892590972213036;
    w[66] =     0.008996329246717397;
    w[67] =     0.009098707399409718;
    w[68] =     0.009199709952114802;
    w[69] =     0.009299321634629343;
    w[70] =     0.009397527387030594;
    w[71] =     0.009494312361953241;
    w[72] =     0.009589661926834022;
    w[73] =     0.009683561666124043;
    w[74] =     0.009775997383468165;
    w[75] =     0.009866955103851452;
    w[76] =     0.009956421075711706;
    w[77] =      0.01004438177301882;
    w[78] =      0.01013082389731963;
    w[79] =      0.01021573437974821;
    w[80] =       0.0102991003830022;
    w[81] =      0.01038090930328312;
    w[82] =      0.01046114877220228;
    w[83] =      0.01053980665865038;
    w[84] =      0.01061687107063194;
    w[85] =      0.01069233035706287;
    w[86] =      0.01076617310953212;
    w[87] =      0.01083838816402652;
    w[88] =      0.01090896460261843;
    w[89] =      0.01097789175511656;
    w[90] =      0.01104515920067912;
    w[91] =      0.01111075676938929;
    w[92] =      0.01117467454379268;
    w[93] =      0.01123690286039691;
    w[94] =      0.01129743231113249;
    w[95] =      0.01135625374477508;
    w[96] =      0.01141335826832922;
    w[97] =      0.01146873724837283;
    w[98] =      0.01152238231236217;
    w[99] =      0.01157428534989815;
    w[100] =      0.01162443851395193;
    w[101] =      0.01167283422205182;
    w[102] =      0.01171946515742932;
    w[103] =      0.01176432427012535;
    w[104] =      0.01180740477805627;
    w[105] =      0.01184870016803913;
    w[106] =      0.01188820419677619;
    w[107] =      0.01192591089179929;
    w[108] =      0.01196181455237226;
    w[109] =      0.01199590975035326;
    w[110] =      0.01202819133101508;
    w[111] =      0.01205865441382472;
    w[112] =      0.01208729439318107;
    w[113] =      0.01211410693911137;
    w[114] =      0.01213908799792579;
    w[115] =      0.01216223379283022;
    w[116] =      0.01218354082449738;
    w[117] =      0.01220300587159574;
    w[118] =      0.01222062599127671;
    w[119] =      0.01223639851961942;
    w[120] =      0.01225032107203351;
    w[121] =      0.01226239154361966;
    w[122] =      0.01227260810948789;
    w[123] =      0.01228096922503318;
    w[124] =      0.01228747362616942;
    w[125] =      0.01229212032952021;
    w[126] =      0.01229490863256759;
    w[127] =      0.01229583811375833;
    w[128] =      0.01229490863256759;
    w[129] =      0.01229212032952021;
    w[130] =      0.01228747362616942;
    w[131] =      0.01228096922503318;
    w[132] =      0.01227260810948789;
    w[133] =      0.01226239154361966;
    w[134] =      0.01225032107203351;
    w[135] =      0.01223639851961942;
    w[136] =      0.01222062599127671;
    w[137] =      0.01220300587159574;
    w[138] =      0.01218354082449738;
    w[139] =      0.01216223379283022;
    w[140] =      0.01213908799792579;
    w[141] =      0.01211410693911137;
    w[142] =      0.01208729439318107;
    w[143] =      0.01205865441382472;
    w[144] =      0.01202819133101508;
    w[145] =      0.01199590975035326;
    w[146] =      0.01196181455237226;
    w[147] =      0.01192591089179929;
    w[148] =      0.01188820419677619;
    w[149] =      0.01184870016803913;
    w[150] =      0.01180740477805627;
    w[151] =      0.01176432427012535;
    w[152] =      0.01171946515742932;
    w[153] =      0.01167283422205182;
    w[154] =      0.01162443851395193;
    w[155] =      0.01157428534989815;
    w[156] =      0.01152238231236217;
    w[157] =      0.01146873724837283;
    w[158] =      0.01141335826832922;
    w[159] =      0.01135625374477508;
    w[160] =      0.01129743231113249;
    w[161] =      0.01123690286039691;
    w[162] =      0.01117467454379268;
    w[163] =      0.01111075676938929;
    w[164] =      0.01104515920067912;
    w[165] =      0.01097789175511656;
    w[166] =      0.01090896460261843;
    w[167] =      0.01083838816402652;
    w[168] =      0.01076617310953212;
    w[169] =      0.01069233035706287;
    w[170] =      0.01061687107063194;
    w[171] =      0.01053980665865038;
    w[172] =      0.01046114877220228;
    w[173] =      0.01038090930328312;
    w[174] =       0.0102991003830022;
    w[175] =      0.01021573437974821;
    w[176] =      0.01013082389731963;
    w[177] =      0.01004438177301882;
    w[178] =     0.009956421075711706;
    w[179] =     0.009866955103851452;
    w[180] =     0.009775997383468165;
    w[181] =     0.009683561666124043;
    w[182] =     0.009589661926834022;
    w[183] =     0.009494312361953241;
    w[184] =     0.009397527387030594;
    w[185] =     0.009299321634629343;
    w[186] =     0.009199709952114802;
    w[187] =     0.009098707399409718;
    w[188] =     0.008996329246717397;
    w[189] =     0.008892590972213036;
    w[190] =     0.008787508259703609;
    w[191] =     0.008681096996256795;
    w[192] =     0.008573373269798925;
    w[193] =     0.008464353366682819;
    w[194] =     0.008354053769225508;
    w[195] =     0.008242491153216323;
    w[196] =     0.008129682385395602;
    w[197] =     0.008015644520904983;
    w[198] =     0.007900394800708624;
    w[199] =     0.007783950648986801;
    w[200] =     0.007666329670501377;
    w[201] =     0.007547549647934514;
    w[202] =     0.007427628539199977;
    w[203] =     0.007306584474728122;
    w[204] =     0.007184435754724984;
    w[205] =     0.007061200846405536;
    w[206] =     0.006936898381201466;
    w[207] =     0.006811547151944815;
    w[208] =     0.006685166110026254;
    w[209] =     0.006557774362530328;
    w[210] =     0.006429391169346602;
    w[211] =     0.006300035940257733;
    w[212] =     0.006169728232005295;
    w[213] =     0.006038487745332765;
    w[214] =     0.005906334322007422;
    w[215] =     0.005773287941820301;
    w[216] =     0.005639368719565862;
    w[217] =     0.005504596902000804;
    w[218] =     0.005368992864783177;
    w[219] =     0.005232577109391968;
    w[220] =     0.005095370260027839;
    w[221] =      0.00495739306049505;
    w[222] =     0.004818666371065699;
    w[223] =     0.004679211165326077;
    w[224] =      0.00453904852700618;
    w[225] =     0.004398199646792759;
    w[226] =     0.004256685819126112;
    w[227] =     0.004114528438981242;
    w[228] =     0.003971748998634907;
    w[229] =     0.003828369084417135;
    w[230] =     0.003684410373449933;
    w[231] =     0.003539894630372244;
    w[232] =     0.003394843704053401;
    w[233] =     0.003249279524294296;
    w[234] =     0.003103224098519095;
    w[235] =     0.002956699508457498;
    w[236] =      0.00280972790682046;
    w[237] =     0.002662331513971696;
    w[238] =      0.00251453261459971;
    w[239] =     0.002366353554396287;
    w[240] =     0.002217816736754017;
    w[241] =     0.002068944619501569;
    w[242] =     0.001919759711713187;
    w[243] =     0.001770284570660304;
    w[244] =     0.001620541799041545;
    w[245] =     0.001470554042778403;
    w[246] =     0.001320343990022177;
    w[247] =       0.0011699343729388;
    w[248] =     0.001019347976427318;
    w[249] =    0.0008686076661194581;
    w[250] =    0.0007177364780061286;
    w[251] =    0.0005667579456482639;
    w[252] =    0.0004156976252681932;
    w[253] =    0.0002645938711908564;
    w[254] =    0.0001136736199914808;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET - Fatal error!\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1 through 33, 63, 64, 65, 127 or 255.\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

double *map ( string code, int element_order )

//****************************************************************************80
//
//  Purpose: 
//
//    MAP returns the interpolation matrix for any available element.
//
//  Formula:
//
//    For an element of order N, we suppose we are given N items of data 
//    Q associated with the nodes.
//
//   Let PHI(J)(R,S) be the Lagrange basis polynomial associated with 
//   node J.  PHI(J)(R,S) is 1 at node J, and 0 at each of the other nodes.
//
//   Let P(R,S) be the polynomial of N terms which interpolates the
//   data Q, that is,
//
//      P(R(J),S(J)) = Q(J)
//
//   where the coordinates of node J are (R(J),S(J)).  Then we know
//   that we can write
//
//     P(R,S) = sum ( 1 <= J <= N ) Q(J) * PHI(J)(R,S)
//
//   But P(R,S) also has a standard representation as
//
//     P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I)
//
//   where REXP(I) and SEXP(I) are the exponents of R and S and
//   the A(I) are the appropriate coefficients.
//
//   The interpolation matrix W allows us to immediately compute
//   the standard basis coefficients A from the data Q to be interpolated
//   using the formula:
//
//      A(I) = sum ( 1 <= J <= N ) W(I,J) * Q(J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T4', 'T6' and 'T10'.
//
//    Input, int N, the order associated with the code.
//
//    Output, double MAP[N*N], the interpolation matrix.
//
{
  double area;
  int i;
  int info;
  int j;
  int *pivot;
  double *r;
  int *rexp;
  double rfact;
  double *s;
  int *sexp;
  double sfact;
  double *v;
  double *w;

  pivot = new int[element_order];
  r = new double[element_order];
  rexp = new int[element_order];
  s = new double[element_order];
  sexp = new int[element_order];
  v = new double[element_order*element_order];
//
//  Get the (R,S) location of the nodes.
//
  node_reference ( code, r, s, &area );
//
//  Get the associated monomials.
//
  poly ( code, rexp, sexp );
//
//  Set up the Vandermonde matrix.
//  Factors of the form 0**0 are to be understood as 1.
//
  for ( i = 0; i < element_order; i++ )
  {
    for ( j = 0; j < element_order; j++ )
    {
      if ( rexp[j] == 0 )
      {
        rfact = 1.0;
      }
      else
      {
        rfact = pow ( r[i], rexp[j] );
      }

      if ( sexp[j] == 0 )
      {
        sfact = 1.0;
      }
      else
      {
        sfact = pow ( s[i], sexp[j] );
      }
      v[i+j*element_order] = rfact * sfact;
    }
  }
//
//  Factor the Vandermonde matrix.
//
  info = r8ge_fa ( element_order, v, pivot );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "MAP - Fatal error!\n";
    cout << "  The Vandermonde matrix is singular.\n";
    exit ( 1 );
  }
//
//  Invert the Vandermonde matrix.
//
  w = r8ge_inverse ( element_order, v, pivot );

  delete [] pivot;
  delete [] r;
  delete [] rexp;
  delete [] s;
  delete [] sexp;
  delete [] v;

  return w;
}
//****************************************************************************80

void map_test ( string code )

//****************************************************************************80
//
//  Purpose: 
//
//    MAP_TEST tests the map routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, the code for the element.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T4', 'T6' and 'T10'.
//
{
  int element_order;
  double *w;

  cout << "\n";
  cout << "  MAP_TEST: The interpolation matrix for element " << code << "\n";

  element_order = order_code ( code );

  w = map ( code, element_order );

  r8mat_print ( element_order, element_order, w, 
    "  The interpolation matrix:" );

  delete [] w;

  return;
}
//****************************************************************************80

double *mass_matrix_t6 ( int node_num, int element_num, int element_node[], 
  double node_xy[] )

//****************************************************************************80
//
//  Purpose: 
//
//    MASS_MATRIX_T6 computes the mass matrix, using 6-node triangles.
//
//  Discussion:
//
//    The mass matrix to be estimated has the form:
//
//      A(I,J) = integral ( PHI(I)(X,Y) * PHI(J)(X,Y) ) d Region
//
//    where PHI(I) and PHI(J) are the shape functions associated with
//    the I-th and J-th variables.
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 May 2007
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
//    Input, int ELEMENT_NODE[6*ELEMENT_NUM], the nodes that make up each element.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Output, double MASS_MATRIX_T6[NODE_NUM*NODE_NUM], the mass matrix.
//
{
# define ELEMENT_ORDER 6

  double *a;
  double area;
  double dwdr[ELEMENT_ORDER];
  double dwds[ELEMENT_ORDER];
  int element;
  int global1;
  int global2;
  int local1;
  int local2;
  int norder;
  int p1;
  int p2;
  int p3;
  int quad;
  int quad_num;
  double r;
  double *rtab;
  int rule;
  double s;
  double *stab;
  double w[ELEMENT_ORDER];
  double *weight;

  a = new double[node_num*node_num];
//
//  Zero out the matrix.
//
  for ( global2 = 0; global2 < node_num; global2++ )
  {
    for ( global1 = 0; global1 < node_num; global1++ )
    {
      a[global1+global2*node_num] = 0.0;
    }
  }
//
//  Get the weights and abscissas for a unit triangle.
//
  rule = 12;
  quad_num = triangle_unit_size ( rule );

  rtab = new double[quad_num];
  stab = new double[quad_num];
  weight = new double[quad_num];

  triangle_unit_set ( rule, rtab, stab, weight );
//
//  For each element.
//
  for ( element = 0; element < element_num; element++ )
  {
    p1 = element_node[0+element*6] - 1;
    p2 = element_node[1+element*6] - 1;
    p3 = element_node[2+element*6] - 1;

    area = 0.5 * r8_abs ( 
        node_xy[0+p1*2] * ( node_xy[1+p2*2] - node_xy[1+p3*2] ) 
      + node_xy[0+p2*2] * ( node_xy[1+p3*2] - node_xy[1+p1*2] ) 
      + node_xy[0+p3*2] * ( node_xy[1+p1*2] - node_xy[1+p2*2] ) );

    if ( area == 0.0 )
    {
      cout << "\n";
      cout << "MASS_MATRIX_T6 - Fatal error!\n";
      cout << "  Zero area for element " << element << "\n";
      cout << "  Node 1 = " << p1 << "\n";
      cout << "  X = " << node_xy[0+p1*2] << "\n";
      cout << "  Y = " << node_xy[1+p1*2] << "\n";
      cout << "  Node 2 = " << p2 << "\n";
      cout << "  X = " << node_xy[0+p2*2] << "\n";
      cout << "  Y = " << node_xy[1+p2*2] << "\n";
      cout << "  Node 3 = " << p3 << "\n";
      cout << "  X = " << node_xy[0+p3*2] << "\n";
      cout << "  Y = " << node_xy[1+p3*2] << "\n";
      exit ( 1 );
    }
//
//  For each quadrature point in the element...
//
    for ( quad = 0; quad < quad_num; quad++ )
    {
      r = rtab[quad];
      s = stab[quad];

      shape_t6 ( r, s, w, dwdr, dwds );
//
//  For each basis function PHI(I) associated with a node in the element,
//
      for ( local1 = 0; local1 < 6; local1++ )
      {
        global1 = element_node[local1+element*6] - 1;
//
//  For each "neighbor" basis function PHI(J) associated with a node in
//  the element.
//
        for ( local2 = 0; local2 < 6; local2++ )
        {
          global2 = element_node[local2+element*6] - 1;

          a[global1+global2*node_num] = a[global1+global2*node_num] 
            + area * weight[quad] * w[local1] * w[local2];
        }
      }
    }
  }

  delete [] rtab;
  delete [] stab;
  delete [] weight;

  return a;
# undef ELEMENT_ORDER
}
//****************************************************************************80

int next_boundary_node ( int node, string code )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE returns the next boundary node in any element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Input, string CODE, the element type.
//
//    Output, int NEXT_BOUNDARY_NODE, the index of the next edge node.
//
{
  int value;

  if ( code == "Q4" )
  {
    value = next_boundary_node_q4 ( node );
  }
  else if ( code == "Q8" )
  {
    value = next_boundary_node_q8 ( node );
  }
  else if ( code == "Q9" )
  {
    value = next_boundary_node_q9 ( node );
  }
  else if ( code == "Q12" )
  {
    value = next_boundary_node_q12 ( node );
  }
  else if ( code == "Q16" )
  {
    value = next_boundary_node_q16 ( node );
  }
  else if ( code == "QL" )
  {
    value = next_boundary_node_ql ( node );
  }
  else if ( code == "T3" )
  {
    value = next_boundary_node_t3 ( node );
  }
  else if ( code == "T4" )
  {
    value = next_boundary_node_t4 ( node );
  }
  else if ( code == "T6" )
  {
    value = next_boundary_node_t6 ( node );
  }
  else if ( code == "T10" )
  {
    value = next_boundary_node_t10 ( node );
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_q4 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_Q4 returns the next boundary node in a Q4 element.
//
//  Element Q4:
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_Q4, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_q8 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_Q8 returns the next boundary node in a Q8 element.
//
//  Element Q8:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8     6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_Q8, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_q9 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_Q9 returns the next boundary node in a Q9 element.
//
//  Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_Q9, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else
  {
    value = 1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_q12 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_Q12 returns the next boundary node in a Q12 element.
//
//  Element Q12:
//
//    |
//    1  9-10-11-12
//    |  |        |
//    |  7        8
//    S  |        |
//    |  5        6
//    |  |        |
//    0  1--2--3--4
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_Q12, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 12;
  }
  else if ( node == 12 )
  {
    value = 11;
  }
  else if ( node == 11 )
  {
    value = 10;
  }
  else if ( node == 10 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_q16 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_Q16 returns the next boundary node in a Q16 element.
//
//  Element Q16:
//
//    |
//    1 13--14--15--16
//    |  |   :   :   |
//    |  |   :   :   |
//    |  9..10..11..12
//    S  |   :   :   |
//    |  |   :   :   |
//    |  5...6...7...8
//    |  |   :   :   |
//    |  |   :   :   |  
//    0  1---2---3---4
//    |
//    +--0-----R-----1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_Q16, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 12;
  }
  else if ( node == 12 )
  {
    value = 16;
  }
  else if ( node == 16 )
  {
    value = 15;
  }
  else if ( node == 15 )
  {
    value = 14;
  }
  else if ( node == 14 )
  {
    value = 13;
  }
  else if ( node == 13 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 5;
  }
  else if ( node == 5 ) 
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_ql ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_QL returns the next boundary node in a QL element.
//
//  Element QL:
//
//    |
//    1  4---5---6
//    |  |       |
//    |  |       |
//    S  |       |
//    |  |       |
//    |  |       |
//    0  1---2---3
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_QL, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_t3 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_T3 returns the next boundary node in a T3 element.
//
//  Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_T3, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_t4 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_T4 returns the next boundary node in a T4 element.
//
//  Element T4:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  | 4 \
//    |  |    \
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
//    23 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_T4, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_t6 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_T6 returns the next boundary node in a T6 element.
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_T6, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 5;
  }
  else if ( node == 5 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 6;
  }
  else if ( node == 6 )
  {
    value = 1;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

int next_boundary_node_t10 ( int node )

//****************************************************************************80
//
//  Purpose: 
//
//    NEXT_BOUNDARY_NODE_T10 returns the next boundary node in a T10 element.
//
//  Element T10:
//
//    |
//    1  10
//    |  |\
//    |  | \
//    |  8  9
//    |  |   \
//    S  |    \
//    |  5  6  7
//    |  |      \
//    |  |       \
//    0  1--2--3--4
//    |
//    +--0----R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE, the index of the current node.  An input
//    value of 0 (or any "unusual" value") indicates that the
//    first edge node is desired.
//
//    Output, int NEXT_BOUNDARY_NODE_T10, the index of the next edge node.
//
{
  int value;

  if ( node == 1 )
  {
    value = 2;
  }
  else if ( node == 2 )
  {
    value = 3;
  }
  else if ( node == 3 )
  {
    value = 4;
  }
  else if ( node == 4 )
  {
    value = 7;
  }
  else if ( node == 7 )
  {
    value = 9;
  }
  else if ( node == 9 )
  {
    value = 10;
  }
  else if ( node == 10 )
  {
    value = 8;
  }
  else if ( node == 8 )
  {
    value = 5;
  }
  else
  {
    value = 1;
  }

  return value;
}
//****************************************************************************80

void node_reference ( string code, double r[], double s[], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE returns the basis nodes for any available element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T4', 'T6' and 'T10'.
//
//    Output, double R[N], S[N], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  if ( code == "Q4" )
  {
    node_reference_q4 ( r, s, area );
  }
  else if ( code == "Q8" )
  {
    node_reference_q8 ( r, s, area );
  }
  else if ( code == "Q9" )
  {
    node_reference_q9 ( r, s, area );
  }
  else if ( code == "Q12" )
  {
    node_reference_q12 ( r, s, area );
  }
  else if ( code == "Q16" )
  {
    node_reference_q16 ( r, s, area );
  }
  else if ( code == "QL" )
  {
    node_reference_ql ( r, s, area );
  }
  else if ( code == "T3" )
  {
    node_reference_t3 ( r, s, area );
  }
  else if ( code == "T4" )
  {
    node_reference_t4 ( r, s, area );
  }
  else if ( code == "T6" )
  {
    node_reference_t6 ( r, s, area );
  }
  else if ( code == "T10" )
  {
    node_reference_t10 ( r, s, area );
  }
  else
  {
    cout << "\n";
    cout << "NODE_REFERENCE - Fatal error!\n";
    cout << "  Illegal value of CODE = " << code << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void node_reference_q4 ( double r[4], double s[4], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_Q4 returns the basis nodes for a 4 node quadrilateral.
//
//  Element Q4:
//
//    |
//    1  4-------3
//    |  |       |
//    |  |       |
//    S  |       |
//    |  |       |
//    |  |       |
//    0  1-------2
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[4], S[4], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_q8 ( double r[8], double s[8], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_Q8 returns the basis nodes for an 8 node quadrilateral.
//
//  Comment:
//
//    This element is known as the quadratic "serendipity" element.
//
//  Element Q8:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8     6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[8], S[8], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0; 
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;
  r[6] = 0.5;
  r[7] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;
  s[4] = 0.0;
  s[5] = 0.5;
  s[6] = 1.0;
  s[7] = 0.5;

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_q9 ( double r[9], double s[9], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_Q9 returns the basis nodes for a 9 node quadrilateral.
//
//  Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R[9], S[9], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;
  r[6] = 0.5;
  r[7] = 0.0;
  r[8] = 0.5;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0;
  s[4] = 0.0;
  s[5] = 0.5;
  s[6] = 1.0;
  s[7] = 0.5;
  s[8] = 0.5;

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_q12 ( double r[12], double s[12], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_Q12 returns the basis nodes for a 12 node quadrilateral.
//
//  Discussion:
//
//    This element is known as the cubic "serendipity" element.
//
//  Element Q12:
//
//    |
//    1  9-10-11-12
//    |  |        |
//    |  7        8
//    S  |        |
//    |  5        6
//    |  |        |
//    0  1--2--3--4
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[12], S[12], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  double a = 0.0;
  double b = 1.0 / 3.0;
  double c = 2.0 / 3.0;
  double d = 1.0;

  r[0] = a;
  r[1] = b;
  r[2] = c;
  r[3] = d;
  r[4] = a;
  r[5] = d;
  r[6] = a;
  r[7] = d;
  r[8] = a;
  r[9] = b;
  r[10] = c;
  r[11] = d;

  s[0] = a;
  s[1] = a;
  s[2] = a;
  s[3] = a;
  s[4] = b;
  s[5] = b;
  s[6] = c;
  s[7] = c;
  s[8] = d;
  s[9] = d;
  s[10] = d;
  s[11] = d;

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_q16 ( double r[16], double s[16], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_Q16 returns the basis nodes for a 16 node quadrilateral.
//
//  Element Q16:
//
//    |
//    1 13--14--15--16
//    |  |   :   :   |
//    |  |   :   :   |
//    |  9..10..11..12
//    S  |   :   :   |
//    |  |   :   :   |
//    |  5...6...7...8
//    |  |   :   :   |
//    |  |   :   :   |  
//    0  1---2---3---4
//    |
//    +--0-----R-----1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[16], S[16], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  int i;
  int j;
  int k;

  k = 0;
  for ( i = 0; i <= 3; i++ )
  {
    for ( j = 0; j <= 3; j++ )
    {
      r[k] = ( double ) ( j ) / 3.0;
      s[k] = ( double ) ( i ) / 3.0;
      k = k + 1;
    }
  }

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_ql ( double r[6], double s[6], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_QL returns the basis nodes for a quadratic/linear.
//
//  Element QL:
//
//    |
//    1  4---5---6
//    |  |       |
//    |  |       |
//    S  |       |
//    |  |       |
//    |  |       |
//    0  1---2---3
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[6], S[6], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 0.5;
  r[2] = 1.0;
  r[3] = 0.0;
  r[4] = 0.5;
  r[5] = 1.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 0.0;
  s[3] = 1.0;
  s[4] = 1.0;
  s[5] = 1.0;

  *area = 1.0;

  return;
}
//****************************************************************************80

void node_reference_t3 ( double r[3], double s[3], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_T3 returns the basis nodes for the 3 node triangle.
//
//  Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[3], S[3], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;

  *area = 0.5;

  return;
}
//****************************************************************************80

void node_reference_t4 ( double r[4], double s[4], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_T4 returns the basis nodes for the 4 node triangle.
//
//  Element T4:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  | 4 \
//    |  |    \
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
//    23 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[4], S[4], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;
  r[3] = 1.0 / 3.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 1.0 / 3.0;

  *area = 0.5;

  return;
}
//****************************************************************************80

void node_reference_t6 ( double r[6], double s[6], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_T6 returns the basis nodes for a 6 node triangle.
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[6], S[6], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  r[1] = 1.0;
  r[2] = 0.0;
  r[3] = 0.5;
  r[4] = 0.5;
  r[5] = 0.0;

  s[0] = 0.0;
  s[1] = 0.0;
  s[2] = 1.0;
  s[3] = 0.0;
  s[4] = 0.5;
  s[5] = 0.5;

  *area = 0.5;

  return;
}
//****************************************************************************80

void node_reference_t10 ( double r[10], double s[10], double *area )

//****************************************************************************80
//
//  Purpose: 
//
//    NODE_REFERENCE_T10 returns the basis nodes for a 10 node triangle.
//
//  Element T10:
//
//    |
//    1  10
//    |  |\
//    |  | \
//    |  8  9
//    |  |   \
//    S  |    \
//    |  5  6  7
//    |  |      \
//    |  |       \
//    0  1--2--3--4
//    |
//    +--0----R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R[10], S[10], the coordinates of the basis nodes.
//
//    Output, double *AREA, the area of the element.
//
{
  r[0] = 0.0;
  s[0] = 0.0;

  r[1] = 1.0 / 3.0;
  s[1] = 0.0;

  r[2] = 2.0 / 3.0;
  s[2] = 0.0;

  r[3] = 1.0;
  s[3] = 0.0;

  r[4] = 0.0;
  s[4] = 1.0 / 3.0;

  r[5] = 1.0 / 3.0;
  s[5] = 1.0 / 3.0;

  r[6] = 2.0 / 3.0;
  s[6] = 1.0 / 3.0;

  r[7] = 0.0;
  s[7] = 2.0 / 3.0;

  r[8] = 1.0 / 3.0;
  s[8] = 2.0 / 3.0;

  r[9] = 0.0;
  s[9] = 1.0;

  *area = 0.5;

  return;
}
//****************************************************************************80

int ns_t6_var_count ( int element_num, int element_node[], int node_num, 
  int var_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    NS_T6_VAR_COUNT counts the Navier Stokes variables on a T6 grid.
//
//  Discussion:
//
//    We are given a mesh of T6 elements, and asked to count, in advance,
//    the number of Navier-Stokes variables associated with the grid.
//    In particular, every node has two velocity variables associated with
//    it, but only a node that is a vertex of the element will also have
//    an associated pressure variable.
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
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Output, int VAR_NODE[NODE_NUM+1], used to find the variables 
//    associated with a given node, which are in VAR in locations 
//    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
//    this array points to the location just after the last location in VAR.
//
//    Output, int NS_T6_VAR_COUNT, the number of variables.
//
{
  int count;
  int element;
  int element_order = 6;
  int node;
  int num;
  int order;
  int var_num;
//
//  Our job is easy once we determine which nodes are vertices.
//  So to begin with, let VAR_NODE count the number of variables
//  associated with each node.
//
  for ( node = 0; node < node_num; node++ )
  {
    var_node[node] = 2;
  }

  for ( element = 0; element < element_num; element++ )
  {
    for ( order = 0; order < 3; order++ )
    {
      node = element_node[order+element*element_order];
      var_node[node-1] = 3;
    }
  }
//
//  Count them.
//
  var_num = 0;
  for ( node = 0; node < node_num; node++ )
  {
    var_num = var_num + var_node[node];
  }
//
//  Make pointers.
//
  count = 1;

  for ( node = 0; node < node_num; node++ )
  {
    num = var_node[node];
    var_node[node] = count;
    count = count + num;
  }
  var_node[node_num] = count;

  return var_num;
}
//****************************************************************************80

int *ns_t6_var_set ( int element_num, int element_node[], int node_num, 
  int var_node[], int var_num )

//****************************************************************************80
//
//  Purpose:
//
//    NS_T6_VAR_SET sets the Navier Stokes variables on a T6 grid.
//
//  Discussion:
//
//    We are given a mesh of T6 elements, and asked to create the natural
//    list of indices for Navier-Stokes variables associated with each node.
//    In particular, every node has two velocity variables associated with
//    it, but only a node that is a vertex of the element will also have
//    an associated pressure variable.
//
//    The hard work has been done for us alread, because the variables
//    have been counted, and the pointers to the occurrence of the
//    first variable associated with each node have been created.
//
//    The indexing of the nodes can be arbitrary, although a bad
//    indexing will result in a miserably large bandwidth (if band
//    storage is being tried for the stiffness matrix).  Here, we
//    simply try to natural ordering, that is, the variables are
//    numbered in order of the node with which they are associated.
//
//    For the Navier Stokes problem on a T6 grid, we take it as
//    understood that each node has either 2 or 3 variables associated
//    with it, that the first two are always the horizontal and
//    then vertical velocity coefficients, and that the third, if
//    present, is a pressure coefficient.
//
//    In other settings, it might be necessary not merely to assign
//    the variables an index, but also to identify them as to type.
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
//  Parameters:
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM]; 
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int VAR_NODE[NODE_NUM+1], used to find the variables 
//    associated with a given node, which are in VAR in locations 
//    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
//    this array points to the location just after the last location in VAR.
//
//    Input, int VAR_NUM, the number of variables.
//
//    Output, int NS_T6_VAR_SET[VAR_NUM], the indexes of the variables, which
//    are simply 1, 2, 3, ..., VAR_NUM.
//
{
  int element_order = 6;
  int i;
  int *var;

  var = new int[var_num];

  for ( i = 0; i < var_num; i++ )
  {
    var[i] = i + 1;
  }
  
  return var;
}
//****************************************************************************80

int order_code ( string code )

//****************************************************************************80
//
//  Purpose: 
//
//    ORDER_CODE returns the order for each element.
//
//  List:
//
//    CODE  Order  Definition
//    ----  -----  ----------
//    Q4     4     4 node linear Lagrange/serendipity quadrilateral;
//    Q8     8     8 node quadratic serendipity quadrilateral;
//    Q9     9     9 node quadratic Lagrange quadrilateral;
//    Q12   12     12 node cubic serendipity quadrilateral;
//    Q16   16     16 node cubic Lagrange quadrilateral;
//    QL     6     6 node linear/quadratic quadrilateral;
//    T3     3     3 node linear triangle;
//    T4     4     4 node bubble triangle;
//    T6     6     6 node quadratic triangle;
//    T10   10     10 node cubic triangle.
// 
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, the code for the element.
//
//    Output, int ORDER_CODE, the order of the element.
//
{
  int value;

  if ( code == "Q4" )
  {
    value = 4;
  }
  else if ( code == "Q8" )
  {
    value = 8;
  }
  else if ( code == "Q9" )
  {
    value = 9;
  }
  else if ( code == "Q12" )
  {
    value = 12;
  }
  else if ( code == "Q16" )
  {
    value = 16;
  }
  else if ( code == "QL" )
  {
    value = 6;
  }
  else if ( code == "T3" )
  {
    value = 3;
  }
  else if ( code == "T4" )
  {
    value = 4;
  }
  else if ( code == "T6" )
  {
    value = 6;
  }
  else if ( code == "T10" )
  {
    value = 10;
  }
  else
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

void physical_to_reference_t3 ( double t[], int n, double phy[], double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHYSICAL_TO_REFERENCE_T3 maps physical points to reference points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (X,Y) in the physical triangle, the routine computes the value
//    of the corresponding image point (XSI,ETA) in reference space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the X and Y coordinates
//    of the vertices.  The vertices are assumed to be the images of
//    (0,0), (1,0) and (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double PHY[2*N], the coordinates of physical points
//    to be transformed.
//
//    Output, double REF[2*N], the coordinates of the corresponding
//    points in the reference space.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {

    ref[0+j*2] = ( ( t[1+2*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] ) ) 
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );

    ref[1+j*2] = ( ( t[0+1*2] - t[0+0*2] ) * ( phy[1+j*2] - t[1+0*2] )   
                 - ( t[1+1*2] - t[1+0*2] ) * ( phy[0+j*2] - t[0+0*2] ) ) 
               / ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2]   - t[0+0*2] )   
                 - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2]   - t[1+0*2] ) );
  }
  return;
}
//****************************************************************************80

void points_plot ( string file_name, int node_num, double node_xy[], 
  bool node_label )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_PLOT plots a pointset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_NAME, the name of the file to create.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the nodes.
//
//    Input, bool NODE_LABEL, is TRUE if the nodes are to be labeled.
//
//  Local parameters:
//
//    int CIRCLE_SIZE, controls the size of the circles depicting
//    the nodes.  Currently set to 5.  3 is pretty small, and 1 is
//    barely visible.
//
{
  int circle_size = 3;
  int delta;
  ofstream file_unit;
  int i;
  int node;
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

  file_unit.open ( file_name.c_str ( ) );

  if ( !file_unit )
  {
    cout << "\n";
    cout << "POINTS_PLOT - Fatal error!\n";
    cout << "  Could not open the output EPS file.\n";
    exit ( 1 );
  }

  file_unit << "%!PS-Adobe-3.0 EPSF-3.0\n";
  file_unit << "%%Creator: points_plot.C\n";
  file_unit << "%%Title: " << file_name << "\n";

  file_unit << "%%Pages: 1\n";
  file_unit << "%%BoundingBox:  "
    << x_ps_min << "  "
    << y_ps_min << "  "
    << x_ps_max << "  "
    << y_ps_max << "\n";
  file_unit << "%%Document-Fonts: Times-Roman\n";
  file_unit << "%%LanguageLevel: 1\n";
  file_unit << "%%EndComments\n";
  file_unit << "%%BeginProlog\n";
  file_unit << "/inch {72 mul} def\n";
  file_unit << "%%EndProlog\n";
  file_unit << "%%Page:      1     1\n";
  file_unit << "save\n";
  file_unit << "%\n";
  file_unit << "% Set the RGB line color to very light gray.\n";
  file_unit << "%\n";
  file_unit << " 0.9000 0.9000 0.9000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "% Draw a gray border around the page.\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  moveto\n";
  file_unit << x_ps_max << "  "
            << y_ps_min << "  lineto\n";
  file_unit << x_ps_max << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_max << "  lineto\n";
  file_unit << x_ps_min << "  "
            << y_ps_min << "  lineto\n";
  file_unit << "stroke\n";
  file_unit << "%\n";
  file_unit << "% Set RGB line color to black.\n";
  file_unit << "%\n";
  file_unit << " 0.0000 0.0000 0.0000 setrgbcolor\n";
  file_unit << "%\n";
  file_unit << "%  Set the font and its size:\n";
  file_unit << "%\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.50 inch scalefont\n";
  file_unit << "setfont\n";
  file_unit << "%\n";
  file_unit << "%  Print a title:\n";
  file_unit << "%\n";
  file_unit << "%  210  702 moveto\n";
  file_unit << "%(Pointset) show\n";
  file_unit << "%\n";
  file_unit << "% Define a clipping polygon\n";
  file_unit << "%\n";
  file_unit << "newpath\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  moveto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << x_ps_max_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_max_clip << "  lineto\n";
  file_unit << x_ps_min_clip << "  "
            << y_ps_min_clip << "  lineto\n";
  file_unit << "clip newpath\n";
//
//  Draw the nodes.
//
  file_unit << "%\n";
  file_unit << "%  Draw filled dots at each node:\n";
  file_unit << "%\n";
  file_unit << "%  Set the color to blue:\n";
  file_unit << "%\n";
  file_unit << "0.000  0.150  0.750  setrgbcolor\n";
  file_unit << "%\n";

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

    file_unit << "newpath  " 
      << x_ps << "  " 
      << y_ps << "  "
      << circle_size << " 0 360 arc closepath fill\n";
  }
//
//  Label the nodes.
//
  file_unit << "%\n";
  file_unit << "%  Label the nodes:\n";
  file_unit << "%\n";
  file_unit << "%  Set the color to darker blue:\n";
  file_unit << "%\n";
  file_unit << "0.000  0.250  0.850  setrgbcolor\n";
  file_unit << "/Times-Roman findfont\n";
  file_unit << "0.20 inch scalefont\n";
  file_unit << "setfont\n";

  file_unit << "%\n";

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

    file_unit << "newpath  " 
      << x_ps     << "  " 
      << y_ps + 5 << "  moveto ("
      << node+1   << ") show\n";
  }

  file_unit << "%\n";
  file_unit << "restore showpage\n";
  file_unit << "%\n";
  file_unit << "% End of page\n";
  file_unit << "%\n";
  file_unit << "%%Trailer\n";
  file_unit << "%%EOF\n";

  file_unit.close ( );

  return;
}
//****************************************************************************80

void poly ( string code, int rexp[], int sexp[] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY returns the polynomial terms associated with any available element.
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
//    'T4', 'T6' and 'T10'.
//
//    Output, int REXP(N), SEXP(N), the powers of R and S associated
//    with each polynomial.
//
{
  if ( code == "Q4" )
  {
    poly_q4 ( rexp, sexp );
  }
  else if ( code == "Q8" )
  {
    poly_q8 ( rexp, sexp );
  }
  else if ( code == "Q9" )
  {
    poly_q9 ( rexp, sexp );
  }
  else if ( code == "Q12" )
  {
    poly_q12 ( rexp, sexp );
  }
  else if ( code == "Q16" )
  {
    poly_q16 ( rexp, sexp );
  }
  else if ( code == "QL" )
  {
    poly_ql ( rexp, sexp );
  }
  else if ( code == "T3" )
  {
    poly_t3 ( rexp, sexp );
  }
  else if ( code == "T4" )
  {
    cout << "\n";
    cout << "POLY - Fatal error!\n";
    cout << "  The T4 element does not follow the pattern!\n";
    exit ( 1 );
  }
  else if ( code == "T6" )
  {
    poly_t6 ( rexp, sexp );
  }
  else if ( code == "T10" )
  {
    poly_t10 ( rexp, sexp );
  }
  else
  {
    cout << "\n";
    cout << "POLY - Fatal error!\n";
    cout << "  Illegal value of CODE = " << code << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void poly_q4 ( int rexp[4], int sexp[4] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_Q4 returns the monomials associated with a 4 node quadrilateral.
//
//  Element Q4:
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
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[4], SEXP[4], the powers of R and S associated
//    with each polynomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 1;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 1;

  return;
}
//****************************************************************************80

void poly_q8 ( int rexp[8], int sexp[8] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_Q8 returns the monomials associated with an 8 node quadrilateral.
//
//  Element Q8:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8     6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[8], SEXP[8], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 1;
  rexp[7] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 2;
  sexp[7] = 1;

  return;
}
//****************************************************************************80

void poly_q9 ( int rexp[9], int sexp[9] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_Q9 returns the monomials associated with a 9 node quadrilateral.
//
//  Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[9], SEXP[9], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 1;
  rexp[7] = 2;
  rexp[8] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 2;
  sexp[7] = 1;
  sexp[8] = 2;

  return;
}
//****************************************************************************80

void poly_q12 ( int rexp[12], int sexp[12] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_Q12 returns the monomials associated with a 12 node quadrilateral.
//
//  Element Q12:
//
//    |
//    1  9-10-11-12
//    |  |        |
//    |  7        8
//    S  |        |
//    |  5        6
//    |  |        |
//    0  1--2--3--4
//    |
//    +--0---R---1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[12], SEXP[12], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;
  rexp[10] = 1;
  rexp[11] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;
  sexp[10] = 3;
  sexp[11] = 1;

  return;
}
//****************************************************************************80

void poly_q16 ( int rexp[16], int sexp[16] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_Q16 returns the monomials associated with a 16 node quadrilateral.
//
//  Element Q16:
//
//    |
//    1 13--14--15--16
//    |  |   :   :   |
//    |  |   :   :   |
//    |  9..10..11..12
//    S  |   :   :   |
//    |  |   :   :   |
//    |  5...6...7...8
//    |  |   :   :   |
//    |  |   :   :   |  
//    0  1---2---3---4
//    |
//    +--0-----R-----1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[16], SEXP[16], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;
  rexp[10] = 1;
  rexp[11] = 2;
  rexp[12] = 3;
  rexp[13] = 2;
  rexp[14] = 3;
  rexp[15] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;
  sexp[10] = 3;
  sexp[11] = 2;
  sexp[12] = 1;
  sexp[13] = 3;
  sexp[14] = 2;
  sexp[15] = 3;

  return;
}
//****************************************************************************80

void poly_ql ( int rexp[6], int sexp[6] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_QL returns the monomials for a quadratic/linear quadrilateral.
//
//  Element QL:
//
//    |
//    1  4---5---6
//    |  |       |
//    |  |       |
//    S  |       |
//    |  |       |
//    |  |       |
//    0  1---2---3
//    |
//    +--0---R---1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[6], SEXP[6], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 1;
  rexp[4] = 2;
  rexp[5] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 1;
  sexp[4] = 0;
  sexp[5] = 1;

  return;
}
//****************************************************************************80

void poly_t3 ( int rexp[3], int sexp[3] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_T3 returns the monomials associated with a 3 node triangle.
//
//  Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
//    0  1-----2
//    |
//    +--0--R--1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[3], SEXP[3], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;

  return;
}
//****************************************************************************80

void poly_t6 ( int rexp[6], int sexp[6] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_T6 returns the monomials associated with a 6 node triangle.
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[6], SEXP[6], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;

  return;
}
//****************************************************************************80

void poly_t10 ( int rexp[10], int sexp[10] )

//****************************************************************************80
//
//  Purpose: 
//
//    POLY_T10 returns the monomials associated with a 10 node triangle.
//
//  Element T10:
//
//    |
//    1  10
//    |  |\
//    |  | \
//    |  8  9
//    |  |   \
//    S  |    \
//    |  5  6  7
//    |  |      \
//    |  |       \
//    0  1--2--3--4
//    |
//    +--0----R---1-->
//
//  Formula:
//
//    Given coefficients A(I), the polynomial interpolant at (R,S) is
//
//      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int REXP[10], SEXP[10], the powers of R and S associated
//    with each monomial.
//
{
  rexp[0] = 0;
  rexp[1] = 0;
  rexp[2] = 1;
  rexp[3] = 0;
  rexp[4] = 1;
  rexp[5] = 2;
  rexp[6] = 0;
  rexp[7] = 1;
  rexp[8] = 2;
  rexp[9] = 3;

  sexp[0] = 0;
  sexp[1] = 1;
  sexp[2] = 0;
  sexp[3] = 2;
  sexp[4] = 1;
  sexp[5] = 0;
  sexp[6] = 3;
  sexp[7] = 2;
  sexp[8] = 1;
  sexp[9] = 0;

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
//    02 April 2005
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
  if ( 0.0 <= x )
  {
    return x;
  } 
  else
  {
    return ( -x );
  }
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the roundoff unit for R8's.
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
//    01 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the double precision round-off unit.
//
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" double precision value.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number, 
//    and is usually defined in math.h, or sometimes in stdlib.h.
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
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  return HUGE_VAL;
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
//        X         R8_NINT
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
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//****************************************************************************80

double r8_power ( double r, int p )

//****************************************************************************80
//
//  Purpose:
//
//    R8_POWER computes the P-th power of R.
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
//  Parameters:
//
//    Input, double R, the base.
//
//    Input, int P, the power, which may be negative.
//
//    Output, double R8_POWER, the value of the P-th power of R.
// 
{
  double value;
//
//  Special case.  R^0 = 1.
//
  if ( p == 0 )
  {
    value = 1.0;
  }
//
//  Special case.  All positive powers of 0 are 0.
//
  else if ( r == 0.0 )
  {
    value = 0.0;
  }
  else if ( 1 <= p )
  {
    value = pow ( r, p );
  }
  else
  {
    value = 1.0 / pow ( r, -p );
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
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
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
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
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
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

int r8ge_fa ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.
//
//  Discussion:
//
//    The R8GE storage format is used for a general M by N matrix.  A physical 
//    storage space is made for each logical entry.  The two dimensional logical
//    array is mapped to a vector, in which storage is by columns.
//
//    R8GE_FA is a simplified version of the LINPACK routine DGEFA.
//
//    The two dimensional array is stored by columns in a one dimensional
//    array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N], the matrix to be factored.
//    On output, A contains an upper triangular matrix and the multipliers
//    which were used to obtain it.  The factorization can be written
//    A = L * U, where L is a product of permutation and unit lower
//    triangular matrices and U is upper triangular.
//
//    Output, int PIVOT[N], a vector of pivot indices.
//
//    Output, int R8GE_FA, singularity flag.
//    0, no singularity detected.
//    nonzero, the factorization failed on the INFO-th step.
//
{
  int i;
  int j;
  int k;
  int l;
  double t;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;

    for ( i = k+1; i <= n; i++ )
    {
      if ( fabs ( a[l-1+(k-1)*n] ) < fabs ( a[i-1+(k-1)*n] ) )
      {
        l = i;
      }
    }

    pivot[k-1] = l;
//
//  If the pivot index is zero, the algorithm has failed.
//
    if ( a[l-1+(k-1)*n] == 0.0 )
    {
      cout << "\n";
      cout << "R8GE_FA - Fatal error!\n";
      cout << "  Zero pivot on step " << k << "\n";
      return k;
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      t              = a[l-1+(k-1)*n];
      a[l-1+(k-1)*n] = a[k-1+(k-1)*n];
      a[k-1+(k-1)*n] = t;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      a[i-1+(k-1)*n] = -a[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t              = a[l-1+(j-1)*n];
        a[l-1+(j-1)*n] = a[k-1+(j-1)*n];
        a[k-1+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + a[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }

    }

  }

  pivot[n-1] = n;

  if ( a[n-1+(n-1)*n] == 0.0 )
  {
    cout << "\n";
    cout << "R8GE_FA - Fatal error!\n";
    cout << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
}
//****************************************************************************80

double *r8ge_inverse ( int n, double a[], int pivot[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.
//
//  Discussion:
//
//    The R8GE storage format is used for a general M by N matrix.  A physical 
//    storage space is made for each logical entry.  The two dimensional logical
//    array is mapped to a vector, in which storage is by columns.
//
//    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
//    DGEDI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix A.
//
//    Input, double A[N*N], the factor information computed by R8GE_FA.
//
//    Input, int PIVOT(N), the pivot vector from R8GE_FA.
//
//    Output, double R8GE_INVERSE[N*N], the inverse matrix.
//
{
  double *b;
  int i;
  int j;
  int k;
  double temp;

  b = new double[n*n];
//
//  Compute Inverse(U).
//
  for ( k = 1; k <= n; k++ )
  {
    for ( i = 1; i <= k-1; i++ )
    {
      b[i-1+(k-1)*n] = -b[i-1+(k-1)*n] / a[k-1+(k-1)*n];
    }
    b[k-1+(k-1)*n] = 1.0 / a[k-1+(k-1)*n];

    for ( j = k+1; j <= n; j++ )
    {
      b[k-1+(j-1)*n] = 0.0;
      for ( i = 1; i <= k; i++ )
      {
        b[i-1+(j-1)*n] = b[i-1+(j-1)*n] + b[i-1+(k-1)*n] * a[k-1+(j-1)*n];
      }
    }
  }
//
//  Multiply Inverse(U) by Inverse(L).
//
  for ( k = n-1; 1 <= k; k-- )
  {
    for ( i = k+1; i <= n; i++ )
    {
      b[i-1+(k-1)*n] = 0.0;
    }

    for ( j = k+1; j <= n; j++ )
    {
      for ( i = 1; i <= n; i++ )
      {
        b[i-1+(k-1)*n] = b[i-1+(k-1)*n] + b[i-1+(j-1)*n] * a[j-1+(k-1)*n];
      }
    }

    if ( pivot[k-1] != k )
    {
      for ( i = 1; i <= n; i++ )
      {
        temp = b[i-1+(k-1)*n];
        b[i-1+(k-1)*n] = b[i-1+(pivot[k-1]-1)*n];
        b[i-1+(pivot[k-1]-1)*n] = temp;
      }

    }

  }

  return b;
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

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    The doubly dimensioned array A is treated as a one dimensional vector,
//    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 August 2003
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
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title to be printed.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2004
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
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title for the matrix.
{
# define INCX 5

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
//  Print the columns of the matrix, in strips of 5.
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
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "  ---\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
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

void reference_sample ( string code, int *seed, double *r, double *s )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_SAMPLE samples a reference element.
//
//  Discussion:
//
//    The routine either samples the unit triangle or the unit square.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include "Q4", "Q8", "Q9", "Q12", "Q16", "QL", "T3", 
//    "T4", "T6" and "T10".
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double *R, *S, a random point in the reference element.
//
{
  if ( code == "Q4" || 
       code == "Q8" || 
       code == "Q9" ||
       code == "Q12" ||
       code == "Q16" ||
       code == "QL" )
  {
    *r = r8_uniform_01 ( seed );
    *s = r8_uniform_01 ( seed );
  }
  else if ( code == "T3" || code == "T4" || code == "T6" || code == "T10" )
  {
    *r = r8_uniform_01 ( seed );
    *s = r8_uniform_01 ( seed );

    if ( 1.0 < *r + *s )
    {
      *r = 1.0 - *r;
      *s = 1.0 - *s;
    }
  }
  else
  {
    cout << "\n";
    cout << "REFERENCE_SAMPLE - Fatal error!\n";
    cout << "  Illegal code = \"" << code << "\".\n";
    exit ( 1 );
  }

  return;
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

void reference_to_physical_t3 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 3 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    Note that this routine may also be appropriate for an order 6
//    triangle, if the mapping between reference and physical space
//    is linear.  This implies, in particular, that the sides of the
//    image triangle are straight and that the "midside" nodes in the
//    physical triangle are halfway along the sides of
//    the physical triangle.
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    24 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0) and
//    (0,1) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = t[i+0*2] * ( 1.0 - ref[0+j*2] - ref[1+j*2] ) 
                 + t[i+1*2] *       + ref[0+j*2]                
                 + t[i+2*2] *                    + ref[1+j*2];
    }
  }

  return;
}
//****************************************************************************80

void reference_to_physical_t6 ( double t[], int n, double ref[], double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 6 physical triangle and a point
//    (XSI,ETA) in the reference triangle, the routine computes the value
//    of the corresponding image point (X,Y) in physical space.
//
//    The mapping from (XSI,ETA) to (X,Y) has the form:
//
//      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
//                 + D1 * XSI    + E1 * ETA     + F1
//
//      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
//                 + D2 * XSI    + E2 * ETA     + F2
//
//  Reference Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*6], the coordinates of the vertices.
//    The vertices are assumed to be the images of (0,0), (1,0),
//    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[2*N], points in the reference triangle.
//
//    Output, double PHY[2*N], corresponding points in the
//    physical triangle.
//
{
  double a[2];
  double b[2];
  double c[2];
  double d[2];
  double e[2];
  double f[2];
  int i;
  int j;

  for ( i = 0; i < 2; i++ )
  {
    a[i] =   2.0 * t[i+0*2] + 2.0 * t[i+1*2]
           - 4.0 * t[i+3*2];

    b[i] =   4.0 * t[i+0*2] 
           - 4.0 * t[i+3*2] + 4.0 * t[i+4*2] - 4.0 * t[i+5*2];

    c[i] =   2.0 * t[i+0*2]                  + 2.0 * t[i+2*2] 
                                             - 4.0 * t[i+5*2];

    d[i] = - 3.0 * t[i+0*2] -       t[i+1*2]
           + 4.0 * t[i+3*2];

    e[i] = - 3.0 * t[i+0*2]                  -       t[i+2*2]
                                             + 4.0 * t[i+5*2];
    f[i] =         t[i+0*2];

  }

  for ( i = 0; i < 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      phy[i+j*2] = a[i] * ref[0+j*2] * ref[0+j*2] 
                 + b[i] * ref[0+j*2] * ref[1+j*2] 
                 + c[i] * ref[1+j*2] * ref[1+j*2]
                 + d[i] * ref[0+j*2]
                 + e[i] * ref[1+j*2] 
                 + f[i];
    }
  }

  return;
}
//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
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
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal. 
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ ) 
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) ) 
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length ) 
  {
    for ( i = nchar; i < s1_length; i++ ) 
    {
      if ( s1[i] != ' ' ) 
      {
        return false;
      }
    } 
  }
  else if ( nchar < s2_length ) 
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' ) 
      {
        return false;
      }
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
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

double serene ( string type, double ve, double vn, double vne, double vnw, 
  double vs, double vse, double vsw, double vw )

//****************************************************************************80
//
//  Purpose: 
//
//    SERENE interpolates data using a Q8 element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string TYPE, tells SERENE the geometry of the
//    finite element that surrounds the point of interest.  The options
//    are displayed in the following table, which suggests the meaning
//    of each option by its position:
//
//        |   |
//     NW * N * NE
//        |   |
//     -*-*-*-*-*-
//        |   |
//      W * C * E
//        |   |
//     -*-*-*-*-*-
//        |   |
//     SW * S * SE
//        |   |
//
//    Input, double VE, VN, VNE, VNW, VS, VSE, VSW, VW,
//    are the values of the function at the nodes to the east,
//    north, northeast, northwest, south, southeast, southwest and
//    west of the point of interest.  If the finite element is of
//    type 'C', then all 8 values are needed.  However, if the
//    finite element is of type 'SE', for instance, then only three
//    values are needed, namely VE, VN, and VNW, since these are
//    the only node positions defined in such a finite element.
//
//    Output, double SERENE, the interpolated value of the
//    function at the point of interest.
//
{
  double eta;
  double pe;
  double pn;
  double pne;
  double pnw;
  double ps;
  double pse;
  double psw;
  double pw;
  double vterp;
  double xsi;
//
//  To make this routine more general, simply pass in the values of XSI
//  and ETA at which the interpolated value is desired.
//
//  By setting XSI = ETA = 0, we are asking for the interpolated value
//  at the center of the finite element.
//
  xsi = 0.0;
  eta = 0.0;
//
//  8 node center
//
//  Polynomial space is spanned by:
//
//         1
//       x    y
//    x^2  xy  y^2
//      x^2y xy^2
//
//
//    ^   1    4--7--3
//    |        !     !
//    E        !     !
//    T   0    8  X  6
//    A        !     !
//    |        !     !
//    V  -1    1--5--2
//
//            -1  0  1
//
//           <---XSI--->
//
  if ( type == "C" )
  {
    psw = - 0.25 * ( 1.0 - xsi ) * ( 1.0 - eta )
      * ( 1.0 + xsi + eta );
    pse = - 0.25 * ( 1.0 + xsi ) * ( 1.0 - eta )
      * ( 1.0 - xsi + eta );
    pne = - 0.25 * ( 1.0 + xsi ) * ( 1.0 + eta )
      * ( 1.0 - xsi - eta );
    pnw = - 0.25 * ( 1.0 - xsi ) * ( 1.0 + eta )
      * ( 1.0 + xsi - eta );
    ps =    0.50 * ( 1.0 - xsi ) * ( 1.0 + xsi )
      * ( 1.0 - eta );
    pe =    0.50 * ( 1.0 + xsi ) * ( 1.0 + eta ) 
      * ( 1.0 - eta );
    pn =    0.50 * ( 1.0 - xsi ) * ( 1.0 + xsi ) 
      * ( 1.0 + eta );
    pw =    0.50 * ( 1.0 - xsi ) * ( 1.0 + eta ) 
      * ( 1.0 - eta );

    vterp = vsw * psw + vse * pse + vne * pne + vnw * pnw 
      + vs * ps + ve * pe + vn * pn + vw * pw;
  }
//
//  5 node side
//
//    ^   1
//    |
//    E
//    T   0    8  X  6
//    A        !     !
//    |        !     !
//    V  -1    1--5--2
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "N" )
  {
    psw =  0.5 * ( xsi - 1.0 ) * ( 1.0 + xsi + eta );
    pse = -0.5 * ( xsi + 1.0 ) * ( 1.0 - xsi + eta );
    ps =  -          ( xsi + 1.0 ) * ( xsi - 1.0 );
    pe =   0.5 * ( xsi + 1.0 ) * ( eta + 1.0 );
    pw =  -0.5 * ( xsi - 1.0 ) * ( eta + 1.0 );

    vterp = vsw * psw + vse * pse + vs * ps + ve * pe + vw * pw;
  }
//
//    ^   1    4--7
//    |        !
//    E        !
//    T   0    8  X
//    A        !
//    |        !
//    V  -1    1--5
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "E" )
  {
    pse =  0.5 * ( eta - 1.0 ) * ( 1.0 + xsi + eta );
    pne = -0.5 * ( eta + 1.0 ) * ( 1.0 + xsi - eta );
    ps =  -0.5 * ( xsi + 1.0 ) * ( eta - 1.0 );
    pn =   0.5 * ( xsi + 1.0 ) * ( eta + 1.0 );
    pw =  -          ( eta + 1.0 ) * ( eta - 1.0 );

    vterp = vse * pse + vne * pne + vs * ps + vn * pn + vw * pw;
  }
//
//  5 node side
//
//    ^   1       7--3
//    |              !
//    E              !
//    T   0       X  6
//    A              !
//    |              !
//    V  -1       5--2
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "W" )
  {
    pse =  0.5 * ( eta - 1.0 ) * ( 1.0 - xsi + eta );
    pne = -0.5 * ( eta + 1.0 ) * ( 1.0 - xsi - eta );
    ps =   0.5 * ( xsi - 1.0 ) * ( eta - 1.0 );
    pe =  -          ( eta - 1.0 ) * ( eta + 1.0 );
    pn =  -0.5 * ( xsi - 1.0 ) * ( eta + 1.0 );

    vterp = vse * pse + vne * pne + vs * ps + ve * pe + vn * pn;
  }
//
//  5 node side
//
//    ^   1    4--7--3
//    |        !     !
//    E        !     !
//    T   0    8  X  6
//    A
//    |
//    V  -1
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "S" )
  {
    pne = -0.5 * ( xsi + 1.0 ) * ( 1.0 - xsi - eta );
    pnw =  0.5 * ( xsi - 1.0 ) * ( 1.0 + xsi - eta );
    pe =  -0.5 * ( eta - 1.0 ) * ( xsi + 1.0 );
    pn =  -          ( xsi + 1.0 ) * ( xsi - 1.0 );
    pw =   0.5 * ( eta - 1.0 ) * ( xsi - 1.0 );

    vterp = vne * pne + vnw * pnw + ve * pe + vn * pn + vw * pw;
  }
//
//  3 node corner
//
//  Polynomial space is spanned by:
//
//         1
//       x    y
//
//
//    ^   1
//    |
//    E
//    T   0    8  X
//    A        !
//    |        !
//    V  -1    1--5
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "NE" )
  {
    psw = -1.0 - xsi - eta;
    ps =   1.0 + xsi;
    pw =   1.0       + eta;

    vterp = vsw * psw + vs * ps + vw * pw;
  }
//
//  3 node corner
//
//  Polynomial space is spanned by:
//
//         1
//       x    y
//
//    ^   1
//    |
//    E
//    T   0       X  6
//    A              !
//    |              !
//    V  -1       5--2
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "NW" )
  {
    pse = 1.0 + xsi - eta;
    ps =  1.0 - xsi;
    pe =  1.0       + eta;

    vterp = vse * pse + vs * ps + ve * pe;
  }
//
//  3 node corner
//
//  Polynomial space is spanned by:
//         1
//       x    y
//
//
//    ^   1    4--7
//    |        !
//    E        !
//    T   0    8  X
//    A
//    |
//    V  -1
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "SE" )
  {
    pnw = - 1.0 - xsi + eta;
    pn =    1.0 + xsi;
    pw =    1.0       - eta;

    vterp = vnw * pnw + vn * pn + vw * pw;
  }
//
//  3 node corner
//
//  Polynomial space is spanned by:
//
//         1
//       x    y
//
//    ^   1       7--3
//    |              !
//    E              !
//    T   0       X  6
//    A
//    |
//    V  -1
//
//            -1  0  1
//
//           <---XSI--->
//
  else if ( type == "SW" )
  {
    pne = - 1.0 + xsi + eta;
    pe =    1.0       - eta;
    pn =    1.0 - xsi;

    vterp = vne * pne + ve * pe + vn * pn;
  }

  return vterp;
}
//****************************************************************************80

void shape ( string code, double r, double s, double t[], 
  double dtdr[], double dtds[] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE evaluates shape functions for any available element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T4', 'T6' and 'T10'.
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[N], the basis functions at the point.
//
//    Output, double DTDR[N], the R basis derivatives at the point.
//
//    Output, double DTDS[N], the S basis derivatives at the point.
//
{
  if ( s_eqi ( code, "Q4" ) )
  {
    shape_q4 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q8" ) )
  {
    shape_q8 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q9" ) )
  {
    shape_q9 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q12" ) )
  {
    shape_q12 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "Q16" ) )
  {
    shape_q16 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "QL" ) )
  {
    shape_ql ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T3" ) )
  {
    shape_t3 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T4" ) )
  {
    shape_t4 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T6" ) )
  {
    shape_t6 ( r, s, t, dtdr, dtds );
  }
  else if ( s_eqi ( code, "T10" ) )
  {
    shape_t10 ( r, s, t, dtdr, dtds );
  }
  else
  {
    cout << "\n";
    cout << "SHAPE - Fatal error!\n";
    cout << "  Unrecognized code = " << code << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void shape_q4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_Q4 evaluates shape functions for a 4 node quadrilateral.
//
//  Element Q4:
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[4], the basis functions at the point.
//
//    Output, double DTDR[4], the R basis derivatives at the point.
//
//    Output, double DTDS[4], the S basis derivatives at the point.
//
{
  t[0] = ( 1.0 - r ) * ( 1.0 - s );
  t[1] =             r   * ( 1.0 - s );
  t[2] =             r   *             s;
  t[3] = ( 1.0 - r ) *             s;

  dtdr[0] = - 1.0 + s;
  dtdr[1] =   1.0 - s;   
  dtdr[2] =             s;
  dtdr[3] =           - s;

  dtds[0] = - 1.0 + r;
  dtds[1] =           - r;
  dtds[2] =             r;
  dtds[3] =   1.0 - r;

  return;
}
//****************************************************************************80

void shape_q8 ( double r, double s, double t[8], double dtdr[8], 
  double dtds[8] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_Q8 evaluates shape functions for an 8 node quadrilateral.
//
//  Comment:
//
//    This element is known as the "serendipity" element.
//
//  Element Q8:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8     6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[8], the basis functions at the point.
//
//    Output, double DTDR[8], the R basis derivatives at the point.
//
//    Output, double DTDS[8], the S basis derivatives at the point.
//
{
  t[0] =                 ( r - 1.0 )     * ( s - 1.0 ) 
    * ( 1.0 - 2.0 * r - 2.0 * s );
  t[1] =             r                       * ( s - 1.0 ) 
    * ( 1.0 - 2.0 * r + 2.0 * s );
  t[2] =             r                   * s                   
    * ( 2.0 * r + 2.0 * s - 3.0 );
  t[3] =                 ( r - 1.0 ) * s                   
    * ( 2.0 * r - 2.0 * s + 1.0 );
  t[4] =   4.0 * r * ( r - 1.0 )     * ( s - 1.0 );
  t[5] = - 4.0 * r                   * s * ( s - 1.0 );
  t[6] = - 4.0 * r * ( r - 1.0 ) * s;   
  t[7] =   4.0 *     ( r - 1.0 ) * s * ( s - 1.0 );

  dtdr[0] = ( s - 1.0 ) * ( - 4.0 * r - 2.0 * s + 3.0 );
  dtdr[1] = ( s - 1.0 ) * ( - 4.0 * r + 2.0 * s + 1.0 );
  dtdr[2] =   s         * (   4.0 * r + 2.0 * s - 3.0 );
  dtdr[3] =   s         * (   4.0 * r - 2.0 * s - 1.0 );
  dtdr[4] =   4.0 * ( 2.0 * r - 1.0 )     * ( s - 1.0 );
  dtdr[5] = - 4.0 *                     s * ( s - 1.0 );
  dtdr[6] = - 4.0 * ( 2.0 * r - 1.0 ) * s;
  dtdr[7] =   4.0 *                     s * ( s - 1.0 );

  dtds[0] = ( r - 1.0 ) * ( - 4.0 * s - 2.0 * r + 3.0 );
  dtds[1] =   r *       (   4.0 * s - 2.0 * r - 1.0 );
  dtds[2] =   r *       (   4.0 * s + 2.0 * r - 3.0 );
  dtds[3] = ( r - 1.0 ) * ( - 4.0 * s + 2.0 * r + 1.0 );
  dtds[4] =   4.0 * r * ( r - 1.0 );
  dtds[5] = - 4.0 * r               * ( 2.0 * s - 1.0 );
  dtds[6] = - 4.0 * r * ( r - 1.0 );
  dtds[7] =   4.0 *     ( r - 1.0 ) * ( 2.0 * s - 1.0 );

  return;
}
//****************************************************************************80

void shape_q9 ( double r, double s, double t[9], double dtdr[9], 
  double dtds[9] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_Q9 evaluates shape functions for a 9 node quadrilateral.
//
//  Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[9], the basis functions at the point.
//
//    Output, double DTDR[9], the R basis derivatives at the point.
//
//    Output, double DTDS[9], the S basis derivatives at the point.
//
{
  t[0] =    4.0 * ( r - 1.0 ) * ( r - 0.5 ) * ( s - 1.0 ) 
            * ( s - 0.5 );
  t[1] =    4.0 * r * ( r - 0.5 ) * ( s - 1.0 ) * ( s - 0.5 );
  t[2] =    4.0 * r * ( r - 0.5 ) * s * ( s - 0.5 );
  t[3] =    4.0 * ( r - 1.0 ) * ( r - 0.5 ) * s * ( s - 0.5 );
  t[4] = -  8.0 * r * ( r - 1.0 ) * ( s - 1.0 ) * ( s - 0.5 );
  t[5] = -  8.0 * r * ( r - 0.5 ) * s * ( s - 1.0 );
  t[6] = -  8.0 * r * ( r - 1.0 ) * s * ( s - 0.5 );
  t[7] = -  8.0 * ( r - 1.0 ) * ( r - 0.5 ) * s * ( s - 1.0 );
  t[8] =   16.0 * r * ( r - 1.0 ) * s * ( s - 1.0 );

  dtdr[0] =   4.0 * ( 2.0 * r - 1.5 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[1] =   4.0 * ( 2.0 * r - 0.5 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[2] =   4.0 * ( 2.0 * r - 0.5 ) * s * ( s - 0.5 );
  dtdr[3] =   4.0 * ( 2.0 * r - 1.5 ) * s * ( s - 0.5 );

  dtdr[4] = - 8.0 * ( 2.0 * r - 1.0 ) * ( s - 1.0 ) 
              * ( s - 0.5 );
  dtdr[5] = - 8.0 * ( 2.0 * r - 0.5 ) * s * ( s - 1.0 );
  dtdr[6] = - 8.0 * ( 2.0 * r - 1.0 ) * s * ( s - 0.5 );
  dtdr[7] = - 8.0 * ( 2.0 * r - 1.5 ) * s * ( s - 1.0 );
  dtdr[8] =  16.0 * ( 2.0 * r - 1.0 ) * s * ( s - 1.0 );

  dtds[0] =   4.0 * ( r - 1.0 ) * ( r - 0.5 ) 
              * ( 2.0 * s - 1.5 );
  dtds[1] =   4.0 * r * ( r - 0.5 ) * ( 2.0 * s - 1.5 );
  dtds[2] =   4.0 * r * ( r - 0.5 ) * ( 2.0 * s - 0.5 );
  dtds[3] =   4.0 * ( r - 1.0 ) * ( r - 0.5 ) 
            * ( 2.0 * s - 0.5 );
  dtds[4] = - 8.0 * r * ( r - 1.0 ) * ( 2.0 * s - 1.5 );
  dtds[5] = - 8.0 * r * ( r - 0.5 ) * ( 2.0 * s - 1.0 );
  dtds[6] = - 8.0 * r * ( r - 1.0 ) * ( 2.0 * s - 0.5 );
  dtds[7] = - 8.0 * ( r - 1.0 ) * ( r - 0.5 ) 
            * ( 2.0 * s - 1.0 );
  dtds[8] =  16.0 * r * ( r - 1.0 ) * ( 2.0 * s - 1.0 );

  return;
}
//****************************************************************************80

void shape_q12 ( double r, double s, double t[12], double dtdr[12], 
  double dtds[12] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_Q12 evaluates shape functions for a 12 node quadrilateral.
//
//  Element Q12:
//
//    |
//    1  9-10-11-12
//    |  |        |
//    |  7        8
//    S  |        |
//    |  5        6
//    |  |        |
//    0  1--2--3--4
//    |
//    +--0---R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[12], the basis functions at the point.
//
//    Output, double DTDR[12], the R basis derivatives at the point.
//
//    Output, double DTDS[12], the S basis derivatives at the point.
//
{
  double a;
  double b;
  double c;
  double corner;
  double d;
  double dcdr;
  double dcds;

  a = 0.0;
  b = 1.0 / 3.0;
  c = 2.0 / 3.0;
  d = 1.0;

  corner = 9.0 * ( ( 2.0 * r - 1.0 ) * ( 2.0 * r - 1.0 )
                     + ( 2.0 * s - 1.0 ) * ( 2.0 * s - 1.0 ) ) 
    - 10.0;

  t[0] =     0.125  * ( r - d ) * ( s - d ) * corner;
  t[1] =  - 13.5    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - d );
  t[2] =    13.5    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - d );
  t[3] =   - 0.125  * ( r - a ) * ( s - d ) * corner;
  t[4] =  - 13.5    * ( r - d ) * ( s - a ) * ( s - c ) * ( s - d );
  t[5] =    13.5    * ( r - a ) * ( s - a ) * ( s - c ) * ( s - d );
  t[6] =    13.5    * ( r - d ) * ( s - a ) * ( s - b ) * ( s - d );
  t[7] =  - 13.5    * ( r - a ) * ( s - a ) * ( s - b ) * ( s - d );
  t[8] =   - 0.125  * ( r - d ) * ( s - a ) * corner;
  t[9] =   13.5    * ( r - a ) * ( r - c ) * ( r - d ) * ( s - a );
  t[10] = - 13.5    * ( r - a ) * ( r - b ) * ( r - d ) * ( s - a );
  t[11] =    0.125  * ( r - a ) * ( s - a ) * corner;
 
  dcdr = 36.0 * ( 2.0 * r - 1.0 );

  dtdr[0] =  0.125 * ( s - d ) * ( ( r - d ) * dcdr + corner );
  dtdr[1] =  - 13.5 * ( s - d ) * ( 3.0 * r * r 
    - 2.0 * ( a + c + d ) * r + a * c + c * d + d * a ); 
  dtdr[2] =    13.5 * ( s - d ) * ( 3.0 * r * r 
    - 2.0 * ( a + b + d ) * r + a * b + b * d + d * a );
  dtdr[3] = - 0.125 * ( s - d ) * ( ( r - a ) * dcdr + corner );
  dtdr[4] = - 13.5 * ( s - a ) * ( s - c ) * ( s - d );
  dtdr[5] =   13.5 * ( s - a ) * ( s - c ) * ( s - d );
  dtdr[6] =   13.5 * ( s - a ) * ( s - b ) * ( s - d );
  dtdr[7] = - 13.5 * ( s - a ) * ( s - b ) * ( s - d );
  dtdr[8] = - 0.125 * ( s - a ) * ( ( r - d ) * dcdr + corner );
  dtdr[9] =   13.5 * ( s - a ) * ( 3.0 * r * r 
    - 2.0 * ( a + c + d ) * r + a * c + c * d + d * a );
  dtdr[10] = - 13.5 * ( s - a ) * ( 3.0 * r * r
    - 2.0 * ( a + b + d ) * r + a * b + b * d + d * a );
  dtdr[11] = 0.125 * ( s - a ) * ( ( r - a ) * dcdr + corner );

  dcds = 36.0 * ( 2.0 * s - 1.0 );

  dtds[0] =  0.125 * ( r - d ) * ( corner + ( s - d ) * dcds );
  dtds[1] =  - 13.5 * ( r - a ) * ( r - c ) * ( r - d );
  dtds[2] =  13.5 * ( r - a ) * ( r - b ) * ( r - d );
  dtds[3] = - 0.125  * ( r - a ) * ( corner + ( s - d ) * dcds );
  dtds[4] =  - 13.5 * ( r - d ) * ( 3.0 * s * s 
    - 2.0 * ( a + c + d ) * s + a * c + c * d + d * a );
  dtds[5] =  13.5 * ( r - a ) * ( 3.0 * s * s 
    - 2.0 * ( a + c + d ) * s + a * c + c * d + d * a );
  dtds[6] =  13.5 * ( r - d ) * ( 3.0 * s * s 
    - 2.0 * ( a + b + d ) * s + a * b + b * d + d * a );
  dtds[7] =  - 13.5 * ( r - a ) * ( 3.0 * s * s 
    - 2.0 * ( a + b + d ) * s + a * b + b * d + d * a );
  dtds[8] =  - 0.125 * ( r - d ) * ( corner + ( s - a ) * dcds );
  dtds[9] = 13.5 * ( r - a ) * ( r - c ) * ( r - d ) ;
  dtds[10] = - 13.5 * ( r - a ) * ( r - b ) * ( r - d ) ;
  dtds[11] = 0.125 * ( r - a ) * ( corner + ( s - a ) * dcds );

  return;
}
//****************************************************************************80

void shape_q16 ( double r, double s, double t[16], double dtdr[16], 
  double dtds[16] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_Q16 evaluates shape functions for a 16 node quadrilateral.
//
//  Diagram:
//
//    |
//    1 13--14--15--16
//    |  |   :   :   |
//    |  |   :   :   |
//    |  9..10..11..12
//    S  |   :   :   |
//    |  |   :   :   |
//    |  5...6...7...8
//    |  |   :   :   |
//    |  |   :   :   |  
//    0  1---2---3---4
//    |
//    +--0-----R-----1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[16], the basis functions at the point.
//
//    Output, double DTDR[16], the R basis derivatives at the point.
//
//    Output, double DTDS[16], the S basis derivatives at the point.
//
{
  double dabc;
  double dabd;
  double dacd;
  double dbcd;
  double ra;
  double rb;
  double rc;
  double rd;
  double sa;
  double sb;
  double sc;
  double sd;

  ra = r - 0.0;
  rb = r - 1.0 / 3.0;
  rc = r - 2.0 / 3.0;
  rd = r - 1.0;

  sa = s - 0.0;
  sb = s - 1.0 / 3.0;
  sc = s - 2.0 / 3.0;
  sd = s - 1.0;

  t[0]  =   (  81.0 / 4.0 ) * rb * rc * rd * sb * sc * sd;
  t[1]  = - ( 243.0 / 4.0 ) * ra * rc * rd * sb * sc * sd;
  t[2]  =   ( 243.0 / 4.0 ) * ra * rb * rd * sb * sc * sd;
  t[3]  = - (  81.0 / 4.0 ) * ra * rb * rc * sb * sc * sd;

  t[4]  = - ( 243.0 / 4.0 ) * rb * rc * rd * sa * sc * sd;
  t[5]  =   ( 729.0 / 4.0 ) * ra * rc * rd * sa * sc * sd;
  t[6]  = - ( 729.0 / 4.0 ) * ra * rb * rd * sa * sc * sd;
  t[7]  =   ( 243.0 / 4.0 ) * ra * rb * rc * sa * sc * sd;

  t[8]  =   ( 243.0 / 4.0 ) * rb * rc * rd * sa * sb * sd;
  t[9]  = - ( 729.0 / 4.0 ) * ra * rc * rd * sa * sb * sd;
  t[10] =   ( 729.0 / 4.0 ) * ra * rb * rd * sa * sb * sd;
  t[11] = - ( 243.0 / 4.0 ) * ra * rb * rc * sa * sb * sd;

  t[12] = - (  81.0 / 4.0 ) * rb * rc * rd * sa * sb * sc;
  t[13] =   ( 243.0 / 4.0 ) * ra * rc * rd * sa * sb * sc;
  t[14] = - ( 243.0 / 4.0 ) * ra * rb * rd * sa * sb * sc;
  t[15] =   (  81.0 / 4.0 ) * ra * rb * rc * sa * sb * sc;

  dbcd = 3.0 * r * r -  4.0 * r       + 11.0 / 9.0;
  dacd = 3.0 * r * r - 10.0 * r / 3.0 +  2.0 / 3.0;
  dabd = 3.0 * r * r -  8.0 * r / 3.0 +  1.0 / 3.0;
  dabc = 3.0 * r * r -  2.0 * r       +  2.0 / 9.0;

  dtdr[0]  =   (  81.0 / 4.0 ) * dbcd * sb * sc * sd;
  dtdr[1]  = - ( 243.0 / 4.0 ) * dacd * sb * sc * sd;
  dtdr[2]  =   ( 243.0 / 4.0 ) * dabd * sb * sc * sd;
  dtdr[3]  = - (  81.0 / 4.0 ) * dabc * sb * sc * sd;
  dtdr[4]  = - ( 243.0 / 4.0 ) * dbcd * sa * sc * sd;
  dtdr[5]  =   ( 729.0 / 4.0 ) * dacd * sa * sc * sd;
  dtdr[6]  = - ( 729.0 / 4.0 ) * dabd * sa * sc * sd;
  dtdr[7]  =   ( 243.0 / 4.0 ) * dabc * sa * sc * sd;
  dtdr[8]  =   ( 243.0 / 4.0 ) * dbcd * sa * sb * sd;
  dtdr[9]  = - ( 729.0 / 4.0 ) * dacd * sa * sb * sd;
  dtdr[10] =   ( 729.0 / 4.0 ) * dabd * sa * sb * sd;
  dtdr[11] = - ( 243.0 / 4.0 ) * dabc * sa * sb * sd;
  dtdr[12] = - (  81.0 / 4.0 ) * dbcd * sa * sb * sc;
  dtdr[13] =   ( 243.0 / 4.0 ) * dacd * sa * sb * sc;
  dtdr[14] = - ( 243.0 / 4.0 ) * dabd * sa * sb * sc;
  dtdr[15] =   (  81.0 / 4.0 ) * dabc * sa * sb * sc;

  dbcd = 3.0 * s * s -  4.0 * s       + 11.0 / 9.0;
  dacd = 3.0 * s * s - 10.0 * s / 3.0 +  2.0 / 3.0;
  dabd = 3.0 * s * s -  8.0 * s / 3.0 +  1.0 / 3.0;
  dabc = 3.0 * s * s -  2.0 * s       +  2.0 / 9.0;

  dtds[0]  =   (  81.0 / 4.0 ) * rb * rc * rd * dbcd;
  dtds[1]  = - ( 243.0 / 4.0 ) * ra * rc * rd * dbcd;
  dtds[2]  =   ( 243.0 / 4.0 ) * ra * rb * rd * dbcd;
  dtds[3]  = - (  81.0 / 4.0 ) * ra * rb * rc * dbcd;
  dtds[4]  = - ( 243.0 / 4.0 ) * rb * rc * rd * dacd;
  dtds[5]  =   ( 729.0 / 4.0 ) * ra * rc * rd * dacd;
  dtds[6]  = - ( 729.0 / 4.0 ) * ra * rb * rd * dacd;
  dtds[7]  =   ( 243.0 / 4.0 ) * ra * rb * rc * dacd;
  dtds[8]  =   ( 243.0 / 4.0 ) * rb * rc * rd * dabd;
  dtds[9]  = - ( 729.0 / 4.0 ) * ra * rc * rd * dabd;
  dtds[10] =   ( 729.0 / 4.0 ) * ra * rb * rd * dabd;
  dtds[11] = - ( 243.0 / 4.0 ) * ra * rb * rc * dabd;
  dtds[12] = - (  81.0 / 4.0 ) * rb * rc * rd * dabc;
  dtds[13] =   ( 243.0 / 4.0 ) * ra * rc * rd * dabc;
  dtds[14] = - ( 243.0 / 4.0 ) * ra * rb * rd * dabc;
  dtds[15] =   (  81.0 / 4.0 ) * ra * rb * rc * dabc;
  
  return;
}
//****************************************************************************80

void shape_ql ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_QL evaluates shape functions for a 6 node quadratic/linear.
//
//  Diagram:
//
//    |
//    1  4--5--6
//    |  |     |
//    |  |     |
//    S  |     |
//    |  |     |
//    |  |     |
//    0  1--2--3
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[6], the basis functions at the point.
//
//    Output, double DTDR[6], the R basis derivatives at the point.
//
//    Output, double DTDS[6], the S basis derivatives at the point.
//
{
  t[0] = - 2.0 *     ( r - 0.5 ) * ( r - 1.0 )     * ( s - 1.0 );
  t[1] =   4.0 * r                   * ( r - 1.0 )     * ( s - 1.0 );
  t[2] = - 2.0 * r * ( r - 0.5 )                       * ( s - 1.0 );
  t[3] =   2.0 *     ( r - 0.5 ) * ( r - 1.0 ) * s;
  t[4] = - 4.0 * r                   * ( r - 1.0 ) * s;
  t[5] =   2.0 * r * ( r - 0.5 )                   * s;

  dtdr[0] = 2.0 * ( - 2.0 * r + 1.5 )     * ( s - 1.0 );
  dtdr[1] = 4.0 * (   2.0 * r - 1.0 )     * ( s - 1.0 );
  dtdr[2] = 2.0 * ( - 2.0 * r + 0.5 )     * ( s - 1.0 );
  dtdr[3] = 2.0 * (   2.0 * r - 1.5 ) * s;
  dtdr[4] = 4.0 * ( - 2.0 * r + 1.0 ) * s;
  dtdr[5] = 2.0 * (   2.0 * r - 0.5 ) * s;

  dtds[0] = - 2.0 *     ( r - 0.5 ) * ( r - 1.0 );
  dtds[1] =   4.0 * r                   * ( r - 1.0 );
  dtds[2] = - 2.0 * r * ( r - 0.5 );
  dtds[3] =   2.0 *     ( r - 0.5 ) * ( r - 1.0 );
  dtds[4] = - 4.0 * r                   * ( r - 1.0 );
  dtds[5] =   2.0 * r * ( r - 0.5 );

  return;
}
//****************************************************************************80

void shape_t3 ( double r, double s, double t[3], double dtdr[3], 
  double dtds[3] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_T3 evaluates shape functions for a 3 node triangle.
//
//  Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[3], the basis functions at the point.
//
//    Output, double DTDR[3], the R basis derivatives at the point.
//
//    Output, double DTDS[3], the S basis derivatives at the point.
//
{
  t[0] = 1.0 - r - s;
  t[1] =           r;
  t[2] =               s;

  dtdr[0] = -1.0;
  dtdr[1] =  1.0;
  dtdr[2] =  0.0;

  dtds[0] = -1.0;
  dtds[1] =  0.0;
  dtds[2] =  1.0;

  return;
}
//****************************************************************************80

void shape_t4 ( double r, double s, double t[4], double dtdr[4], 
  double dtds[4] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_T4 evaluates shape functions for a T4 triangle.
//
//  Element T4:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  | 4 \
//    |  |    \
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
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[4], the basis functions at the point.
//
//    Output, double DTDR[3], the R basis derivatives at the point.
//
//    Output, double DTDS[3], the S basis derivatives at the point.
//
{
  t[0] = ( 1.0 - 9.0 * r * s ) * ( 1.0 - r - s );
  t[1] = r * ( 1.0 - 9.0 * ( 1.0 - r - s ) * s );
  t[2] = s * ( 1.0 - 9.0 * ( 1.0 - r - s ) * r );
  t[3] = 27.0            * ( 1.0 - r - s ) * r * s;

  dtdr[0] = -1.0 +  9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[1] =  1.0 +  9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[2] =         9.0 * ( - s + 2.0 * r * s + s * s );
  dtdr[3] =      - 27.0 * ( - s + 2.0 * r * s + s * s );

  dtds[0] = -1.0 +  9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[1] =         9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[2] =  1.0 +  9.0 * ( - r + r * r + 2.0 * r * s );
  dtds[3] =      - 27.0 * ( - r + r * r + 2.0 * r * s );

  return;
}
//****************************************************************************80

void shape_t6 ( double r, double s, double t[6], double dtdr[6], 
  double dtds[6] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_T6 evaluates shape functions for a 6 node triangle.
//
//  Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[6], the basis functions at the point.
//
//    Output, double DTDR[6], the R basis derivatives at the point.
//
//    Output, double DTDS[6], the S basis derivatives at the point.
//
{
  t[0] = 2.0 *     ( 1.0 - r - s ) * ( 0.5 - r - s );
  t[1] = 2.0 * r * ( r - 0.5 );
  t[2] = 2.0 * s * ( s - 0.5 );
  t[3] = 4.0 * r * ( 1.0 - r - s );
  t[4] = 4.0 * r * s;
  t[5] = 4.0 * s * ( 1.0 - r - s );

  dtdr[0] = - 3.0 + 4.0 * r + 4.0 * s;
  dtdr[1] = - 1.0 + 4.0 * r;
  dtdr[2] =   0.0;
  dtdr[3] =   4.0 - 8.0 * r - 4.0 * s;
  dtdr[4] =                           4.0 * s;
  dtdr[5] =                         - 4.0 * s;

  dtds[0] = - 3.0 + 4.0 * r + 4.0 * s;
  dtds[1] =   0.0;
  dtds[2] = - 1.0               + 4.0 * s;
  dtds[3] =           - 4.0 * r;
  dtds[4] =             4.0 * r;
  dtds[5] =   4.0 - 4.0 * r - 8.0 * s;

  return;
}
//****************************************************************************80

void shape_t10 ( double r, double s, double t[10], double dtdr[10], 
  double dtds[10] )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_T10 evaluates shape functions for a 10 node triangle.
//
//  Diagram:
//
//    |
//    1  10
//    |  |\
//    |  | \
//    |  8  9
//    |  |   \
//    S  |    \
//    |  5  6  7
//    |  |      \
//    |  |       \
//    0  1--2--3--4
//    |
//    +--0----R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the reference coordinates of a point.
//
//    Output, double T[10], the basis functions at the point.
//
//    Output, double DTDR[10], the R basis derivatives at the point.
//
//    Output, double DTDS[10], the S basis derivatives at the point.
//
{
  double a;
  double b;
  double c;

  a = 1.0 / 3.0;
  b = 2.0 / 3.0;
  c = 1.0;

  t[0] = 4.5 * ( a - r - s ) * ( b - r - s ) * ( c - r - s );
  t[1] = 13.5 * r * ( b - r - s ) * ( c - r - s );
  t[2] = - 13.5 * r * ( a - r ) * ( c - r - s );
  t[3] = 4.5 * r * ( a - r ) * ( b - r );
  t[4] = 13.5 * s * ( b - r - s ) * ( c - r - s );
  t[5] = 27.0 * r * s * ( c - r - s );
  t[6] = 13.5 * r * s * ( r - a );
  t[7] = 13.5 * s * ( s - a ) * ( c - r - s );
  t[8] = 13.5 * r * s * ( s - a );
  t[9] = 4.5 * s * ( a - s ) * ( b - s );

  dtdr[0] = 4.5 * ( ( a - s ) * ( 2.0 * r - c - b + 2.0 * s ) 
    - ( s - b ) * ( s - c ) - 2.0 * ( 2.0 * s - b - c ) * r 
    - 3.0 * r * r );
  dtdr[1] = 13.5 * ( 
    ( s - b ) * ( s - c ) + 2.0 * ( 2.0 * s - b - c ) * r 
    + 3.0 * r * r );
  dtdr[2] = - 13.5 * ( a * ( c - s ) + 2.0 * ( s - a - c ) * r 
    + 3.0 * r * r );
  dtdr[3] = 4.5 * ( a * b - 2.0 * ( a + b ) * r + 3.0 * r * r );
  dtdr[4] = 13.5 * s * ( 2.0 * s - b - c + 2.0 * r );
  dtdr[5] = 27.0 * s * ( c - s - 2.0 * r );
  dtdr[6] = 13.5 * s * ( 2.0 * r - a );
  dtdr[7] = - 13.5 * s * ( s - a );
  dtdr[8] = 13.5 * s * ( s - a );
  dtdr[9] = 0.0;

  dtds[0] = 4.5 * ( ( a - r ) * ( 2.0 * s - c - b + 2.0 * r ) 
    - ( r - b ) * ( r - c ) - 2.0 * ( 2.0 * r - b - c ) * s 
    - 3.0 * s * s );
  dtds[1] = 13.5 * r * ( 2.0 * s + 2.0 * r - b - c );
  dtds[2] = 13.5 * r * ( a - r );
  dtds[3] = 0.0;
  dtds[4] = 13.5 * ( ( r - b ) * ( r - c ) + 
    2.0 * ( 2.0 * r - b - c ) * s + 3.0 * s * s );
  dtds[5] = 27.0 * r * ( c - r - 2.0 * s );
  dtds[6] = 13.5 * r * ( r - a );
  dtds[7] = - 13.5 * ( a * ( c - r ) + 2.0 * ( r - c - a ) * s 
    + 3.0 * s * s );
  dtds[8] = 13.5 * r * ( 2.0 * s - a );
  dtds[9] = 4.5 * ( a * b - 2.0 * ( a + b ) * s + 3.0 * s * s );

  return;
}
//****************************************************************************80

void shape_test ( string code )

//****************************************************************************80
//
//  Purpose: 
//
//    SHAPE_TEST verifies the shape function values at the basis nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string CODE, identifies the element to be used.
//    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
//    'T3', 'T6' and 'T10'.
//
{
  double area;
  double *dtdr;
  double *dtds;
  int element_order;
  int i;
  int j;
  double *r;
  double rsum;
  double *s;
  double ssum;
  double *t;

  cout << "\n";
  cout << "  SHAPE_TEST: Verify shape functions of type " << code << "\n";

  element_order = order_code ( code );

  dtdr = new double[element_order];
  dtds = new double[element_order];
  r = new double[element_order];
  s = new double[element_order];
  t = new double[element_order];

  node_reference ( code, r, s, &area );

  cout << "\n";
  cout << "  Element order = " << element_order << "\n";
  cout << "  Basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  for ( i = 0; i < element_order; i++ )
  {
    shape ( code, r[i], s[i], t, dtdr, dtds );
    for ( j = 0; j < element_order; j++ )
    {
      cout << "  " << setw(7) << t[j];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The R and S derivatives should sum to 0.\n";
  cout << "\n";
  cout << "  dTdR sum, dTdS sum:\n";
  cout << "\n";
  for ( i = 0; i < element_order; i++ )
  {
    shape ( code, r[i], s[i], t, dtdr, dtds );
    rsum = 0.0;
    for ( j = 0; j < element_order; j++ )
    {
      rsum = rsum + dtdr[j];
    }
    ssum = 0.0;
    for ( j = 0; j < element_order; j++ )
    {
      ssum = ssum + dtds[j];
    }
    cout << "  " << setw(14) << rsum
         << "  " << setw(14) << ssum << "\n";
  }

  delete [] dtdr;
  delete [] dtds;
  delete [] r;
  delete [] s;
  delete [] t;

  return;
}
//****************************************************************************80

int sphere_grid_element_num ( string code, int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_ELEMENT_NUM returns the number of elements in a sphere grid.
//
//  Discussion:
//
//    The number of elements generated will be NELEMX * NELEMY for
//    quadrilaterals, or 2 * NELEMX * ( NELEMY - 1 ) for triangles.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
//
//    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
//    X and Y directions.  
//
//    Output, int SPHERE_GRID_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  if ( code == "Q4" )
  {
    element_num = sphere_grid_q4_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q9" )
  {
    element_num = sphere_grid_q9_element_num ( nelemx, nelemy );
  }
  else if ( code == "Q16" )
  {
    element_num = sphere_grid_q16_element_num ( nelemx, nelemy );
  }
  else if ( code == "T3" )
  {
    element_num = sphere_grid_t3_element_num ( nelemx, nelemy );
  }
  else if ( code == "T6" )
  {
    element_num = sphere_grid_t6_element_num ( nelemx, nelemy );
  }
  else
  {
    cout << "\n";
    cout << "SPHERE_GRID_ELEMENT_NUM - Fatal error!\n";
    cout << "  Illegal value of CODE = \"" << code << "\".\n";
    element_num = -1;
    exit ( 1 );
  }

  return element_num;
}
//****************************************************************************80

int sphere_grid_node_num ( string code, int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_NODE_NUM returns the number of nodes in a sphere grid.
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
//  Parameters:
//
//    Input, string CODE, identifies the element desired.
//    Legal values include 'Q4', 'Q9', 'Q16', 'T3', 'T6'.
//
//    Input, int NELEMX, NELEMY, the number of quadrilaterals along the
//    X and Y directions.  
//
//    Output, int SPHERE_GRID_NODE_NUM, the number of elements in the grid.
//
{
  int node_num;

  if ( code == "Q" )
  {
    node_num = sphere_grid_q4_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q9" )
  {
    node_num = sphere_grid_q9_node_num ( nelemx, nelemy );
  }
  else if ( code == "Q16" )
  {
    node_num = sphere_grid_q16_node_num ( nelemx, nelemy );
  }
  else if ( code == "T3" )
  {
    node_num = sphere_grid_t3_node_num ( nelemx, nelemy );
  }
  else if ( code == "T6" )
  {
    node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  }
  else
  {
    cout << "\n";
    cout << "SPHERE_GRID_NODE_NUM - Fatal error!\n";
    cout << "  Illegal value of CODE = \"" << code << "\".\n";
    node_num = -1;
    exit ( 1 );
  }

  return node_num;
}
//****************************************************************************80

int *sphere_grid_q4_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q4_ELEMENT produces a Q4 sphere grid.
//
//  Discussion:
//
//    This would be the same as the grid in a plane, except that all the
//    nodes along the bottom edge are identified (replaced by a single node
//    that is the south pole) and similarly for the top edge, and the
//    nodes on the extreme right edge are identified pairwise with those 
//    on the extreme left edge.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  1,  3,  2;
//         1,  1,  4,  3;
//         1,  1,  2,  4;
//         2,  3,  6,  5;
//         3,  4,  7,  6;
//         4,  2,  5,  7;
//         5,  6,  9,  8;
//         6,  7, 10,  9;
//         7,  5,  8, 10;
//         8,  9, 11, 11;
//         9, 10, 11, 11;
//        10,  8, 11, 11;
//
//  Grid:
//
//   11----11----11----11
//    |     |     |     |
//    | E10 | E11 | E12 |
//    |     |     |     |
//    8-----9----10-----8
//    |     |     |     |
//    | E7  | E8  | E9  |
//    |     |     |     |
//    5-----6-----7-----5
//    |     |     |     |
//    | E4  | E5  | E6  |
//    |     |     |     |
//    2-----3-----4-----2
//    |     |     |     |
//    | E1  | E2  | E3  |
//    |     |     |     |
//    1-----1-----1-----1
//
//  Reference Element Q4:
//
//    |
//    1  4------3
//    |  |      |
//    S  |      |
//    |  |      |
//    |  |      |
//    0  1------2
//    |
//    +--0--R---1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int SPHERE_GRID_Q4_ELEMENT[4*NELEMX*NELEMY], the nodes 
//    that form each element.
//
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 4;
  int i;
  int j;

  element_num = sphere_grid_q4_element_num ( nelemx, nelemy );

  element_node = new int[element_order*element_num];
  element = 0;

  for( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * nelemx + 2 - nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + i - 1;

      element_node[0+element*element_order] = base2;
      if ( i < nelemx )
      {
        element_node[1+element*element_order] = base2 + 1;
      }
      else
      {
        element_node[1+element*element_order] = base1;
      }
      element_node[2+element*element_order] = 
        element_node[1+element*element_order] + nelemx;
      element_node[3+element*element_order] = 
        element_node[0+element*element_order] + nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[2+element*element_order] = base1 + nelemx;
        element_node[3+element*element_order] = base1 + nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
//****************************************************************************80

int sphere_grid_q4_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q4_ELEMENT_NUM counts the elements in a Q4 sphere grid.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, int SPHERE_GRID_Q4_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int sphere_grid_q4_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
// 
//    SPHERE_GRID_Q4_NODE_NUM counts nodes in a Q4 sphere grid.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_Q4_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = nelemx * ( nelemy - 1 ) + 2;

  return node_num;
}
//****************************************************************************80

double *sphere_grid_q4_node_xyz ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q4_NODE_XYZ produces node coordinates for a Q4 sphere grid.
//
//  Discussion:
//
//    The number of nodes to be generated is
//
//      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, double SPHERE_GRID_Q4_NODE_XYZ[3*NODE_NUM], 
//    the node coordinates.
//
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  node_xyz = new double[3*node_num];
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;    

  for ( j = nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( nelemy );

    for ( i = 1; i <= nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;  
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
//****************************************************************************80

int *sphere_grid_q9_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q9_ELEMENT produces a Q9 sphere grid.
//
//  Discussion:
//
//    This would be the same as the grid in a plane, except that all the
//    nodes along the bottom edge are identified (replaced by a single node
//    that is the south pole) and similarly for the top edge, and the
//    nodes on the extreme right edge are identified pairwise with those 
//    on the extreme left edge.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  1, 10,  8,  1,  4,  9,  2,  3;
//         1,  1, 12, 10,  1,  6, 11,  4,  5;
//         1,  1,  8, 12,  1,  2, 13,  6,  7;
//         8, 10, 22, 20,  9, 16, 21, 14, 15;
//        10, 12, 24, 22, 11, 18, 23, 16, 17;
//        12,  8, 20, 24, 13, 14, 25, 18, 19;
//        20, 22, 34, 32, 21, 28, 33, 26, 27;
//        22, 24, 36, 34, 23, 30, 35, 28, 29;
//        24, 20, 32, 36, 25, 26, 37, 30, 31;
//        32, 34, 44, 44, 33, 40, 44, 38, 39;
//        34, 36, 44, 44, 35, 42, 44, 40, 41;
//        36, 32, 44, 44, 37, 38, 44, 42, 43;
//
//  Grid:
//
//   44-44-44-44-44-44-44
//    |     |     |     |
//   38 39 40 41 42 43 38
//    |     |     |     |
//   32-33-34-35-36-37-32
//    |     |     |     |
//   26 27 28 29 30 31 26
//    |     |     |     |
//   20-21-22-23-24-25-20
//    |     |     |     |
//   14 15 16 17 18 19 14
//    |     |     |     |
//    8--9-10-11-12-13--8
//    |     |     |     |
//    2  3  4  5  6  7  2
//    |     |     |     |
//    1--1--1--1--1--1--1
//
//  Reference Element Q9:
//
//    |
//    1  4--7--3
//    |  |     |
//    |  |     |
//    S  8  9  6
//    |  |     |
//    |  |     |
//    0  1--5--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_Q9_ELEMENT[9*NELEMX*NELEMY], 
//    the nodes that form each element.
//
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 9;
  int i;
  int j;

  element_num = sphere_grid_q9_element_num ( nelemx, nelemy );
  element_node = new int[element_order*element_num];
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 2 * ( i - 1 );

      element_node[0+element*element_order] = base2;
      element_node[4+element*element_order] = base2 + 1;

      if ( i < nelemx )
      {
        element_node[1+element*element_order] = base2 + 2;
      }
      else
      {
        element_node[1+element*element_order] = base1;
      }

      element_node[7+element*element_order] = 
        element_node[0+element*element_order] + 2 * nelemx;
      element_node[8+element*element_order] = 
        element_node[4+element*element_order] + 2 * nelemx;
      element_node[5+element*element_order] = 
        element_node[1+element*element_order] + 2 * nelemx;

      element_node[3+element*element_order] = 
        element_node[7+element*element_order] + 2 * nelemx;
      element_node[6+element*element_order] = 
        element_node[8+element*element_order] + 2 * nelemx;
      element_node[2+element*element_order] = 
        element_node[5+element*element_order] + 2 * nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[4+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[3+element*element_order] = base1 + 4 * nelemx;
        element_node[6+element*element_order] = base1 + 4 * nelemx;
        element_node[2+element*element_order] = base1 + 4 * nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
//****************************************************************************80

int sphere_grid_q9_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q9_ELEMENT_NUM counts the elements in a Q9 sphere grid.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, int SPHERE_GRID_Q9_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;
 
  return element_num;
}
//****************************************************************************80

int sphere_grid_q9_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q9_NODE_NUM counts nodes in a Q9 sphere grid.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_Q9_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

  return node_num;
}
//****************************************************************************80

double *sphere_grid_q9_node_xyz ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q9_NODE_XYZ produces node coordinates for a Q9 sphere grid.
//
//  Discussion:
//
//    The number of nodes to be generated is
//
//      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, double SPHERE_GRID_Q9_NODE_XYZ[3*NODE_NUM], 
//    the node coordinates.
//
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_q9_node_num ( nelemx, nelemy );
  node_xyz = new double[3*node_num];
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = 2 * nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 2 * nelemy );

    for ( i = 1; i <= 2 * nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 2 * nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
//****************************************************************************80

int *sphere_grid_q16_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q16_ELEMENT produces a Q16 sphere grid.
//
//  Discussion:
//
//    This would be the same as the grid in a plane, except that all the
//    nodes along the bottom edge are identified (replaced by a single node
//    that is the south pole) and similarly for the top edge, and the
//    nodes on the extreme right edge are identified pairwise with those 
//    on the extreme left edge.
//
//  Example:
//
//    Input:
//
//      NELEMX = 2, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  1,  1,  1,  2,  3,  4,  5,  8,  9, 10, 11, 14, 15, 16, 17;
//         1,  1,  1,  1,  5,  6,  7,  2, 11, 12, 13,  8, 17, 18, 19, 14;
//        14, 15, 16, 17, 20, 21, 22, 23, 26, 27, 28, 29, 32, 32, 32, 32;
//        17, 18, 19, 14, 23, 24, 25, 20, 29, 30, 31, 26, 32, 32, 32, 32.
//
//  Grid:
//
//   32-32-32-32-32-32-32
//    |        |        |
//    |        |        |
//   26 27 28 29 30 31 26
//    |        |        |
//    |        |        |
//   20 21 22 23 24 25 20
//    |        |        |
//    | E3     | E4     |
//   14-15-16-17-18-19-14
//    |        |        |
//    |        |        |
//    8  9 10 11 12 13  8
//    |        |        |
//    |        |        |
//    2  3  4  5  6  7  2
//    |        |        |
//    | E1     | E2     |
//    1--1--1--1--1--1--1
//
//  Reference Element Q16:
//
//    |
//    1 13--14--15--16
//    |  |   :   :   |
//    |  |   :   :   |
//    |  9..10..11..12
//    S  |   :   :   |
//    |  |   :   :   |
//    |  5...6...7...8
//    |  |   :   :   |
//    |  |   :   :   |
//    0  1---2---3---4
//    |
//    +--0-----R-----1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  The number of elements generated will be
//    NELEMX * NELEMY.
//
//    Output, int SPHERE_GRID_Q16_ELEMENT[16*NELEMX*NELEMY], 
//    the nodes that form each element.
//
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 16;
  int i;
  int j;

  element_num = sphere_grid_q16_element_num ( nelemx, nelemy );
  element_node = new int[element_order*nelemx*nelemy];

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 3 * ( 3 * nelemx ) + 2 - 3 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 3 * ( i - 1 );

      element_node[0+element*element_order] = base2;
      element_node[1+element*element_order] = base2 + 1;
      element_node[2+element*element_order] = base2 + 2;

      if ( i < nelemx )
      {
        element_node[3+element*element_order] = base2 + 3;
      }
      else
      {
        element_node[3+element*element_order] = base1;
      }

      element_node[4+element*element_order] = 
        element_node[0+element*element_order] + 3 * nelemx;
      element_node[5+element*element_order] = 
        element_node[1+element*element_order] + 3 * nelemx;
      element_node[6+element*element_order] = 
        element_node[2+element*element_order] + 3 * nelemx;
      element_node[7+element*element_order] = 
        element_node[3+element*element_order] + 3 * nelemx;

      element_node[8+element*element_order] = 
        element_node[4+element*element_order] + 3 * nelemx;
      element_node[9+element*element_order] = 
        element_node[5+element*element_order] + 3 * nelemx;
      element_node[10+element*element_order] = 
        element_node[6+element*element_order] + 3 * nelemx;
      element_node[11+element*element_order] = 
        element_node[7+element*element_order] + 3 * nelemx;

      element_node[12+element*element_order] = 
        element_node[8+element*element_order] + 3 * nelemx;
      element_node[13+element*element_order] = 
        element_node[9+element*element_order] + 3 * nelemx;
      element_node[14+element*element_order] = 
        element_node[10+element*element_order] + 3 * nelemx;
      element_node[15+element*element_order] = 
        element_node[11+element*element_order] + 3 * nelemx;

      if ( j == 1 )
      {
        element_node[0+element*element_order] = 1;
        element_node[1+element*element_order] = 1;
        element_node[2+element*element_order] = 1;
        element_node[3+element*element_order] = 1;
      }
      else if ( j == nelemy )
      {
        element_node[12+element*element_order] = base1 + 9 * nelemx;
        element_node[13+element*element_order] = base1 + 9 * nelemx;
        element_node[14+element*element_order] = base1 + 9 * nelemx;
        element_node[15+element*element_order] = base1 + 9 * nelemx;
      }
      element = element + 1;
    }
  }
  return element_node;
}
//****************************************************************************80

int sphere_grid_q16_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q16_ELEMENT_NUM counts the elements in a Q16 sphere grid.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 2
//
//    Output:
//
//      ELEMENT_NUM = NELEMX * NELEMY = 6
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, int SPHERE_GRID_Q16_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = nelemx * nelemy;

  return element_num;
}
//****************************************************************************80

int sphere_grid_q16_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q16_NODE_NUM counts nodes in a Q16 sphere grid.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_Q16_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 9 * nelemx * nelemy - 3 * nelemx + 2;

  return node_num;
}
//****************************************************************************80

double *sphere_grid_q16_node_xyz ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_Q16_NODE_XYZ produces node coordinates for a Q16 sphere grid.
//
//  Discussion:
//
//    The number of nodes to be generated is
//
//      NODE_NUM = 9 * NELEMX * NELEMY - 3 * NELEMX + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, double SPHERE_GRID_Q16_NODE_XYZ[3*NODE_NUM], the node coordinates.
//
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_q16_node_num ( nelemx, nelemy );
  node_xyz = new double[3*node_num];
  node = 0;

  for ( j = 3 * nelemy + 1; 1 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 3 * nelemy );

    if ( j == 1 )
    {
      node_xyz[0+node*3] =  0.0;
      node_xyz[1+node*3] =  0.0;
      node_xyz[2+node*3] =  1.0;
      node = node + 1;
    }
    else if ( j < 3 * nelemy + 1 )
    {
      for ( i = 1; i <= 3 * nelemx; i++ )
      {
        theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 3 * nelemx );

        node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
        node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
        node_xyz[2+node*3] =                 cos ( phi );
        node = node + 1;      
      }
    }
    else if ( j == 3 * nelemy + 1 )
    {
      node_xyz[0+node*3] =  0.0;
      node_xyz[1+node*3] =  0.0;
      node_xyz[2+node*3] = -1.0;
      node = node + 1;
    }
  }
  return node_xyz;
}
//****************************************************************************80

int *sphere_grid_t3_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T3_ELEMENT produces a T3 sphere grid.
//
//  Discussion:
//
//    This would be the same as the grid in a plane, except that all the
//    nodes along the bottom edge are identified (replaced by a single node
//    that is the south pole) and similarly for the top edge, and the
//    nodes on the extreme right edge are identified pairwise with those 
//    on the extreme left edge.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NODE =
//         1,  3,  2;
//         1,  4,  3;
//         1,  2,  4;
//         2,  3,  5
//         6,  5,  3
//         3,  4,  6
//         7,  6,  4;
//         4,  2,  7;
//         5,  7,  2;
//         5,  6,  8;
//         9,  8,  6;
//         6,  7,  9;
//        10,  9,  7;
//         7,  5, 10;
//         8, 10,  5;
//         8,  9, 11;
//         9, 10, 11;
//        10,  8, 11;
//
//  Grid:
//
//   11    11    11    11
//    | \   | \   | \   |
//    |  \  |  \  |  \  |
//    |E16\ |E17 \|E18\ |
//    8-----9----10-----8
//    | \E11| \E13| \E15|
//    |  \  |  \  |  \  |
//    |E10\ |E12\ |E14\ |
//    5-----6-----7-----5
//    | \E5 | \E7 | \E9 |
//    |  \  |  \  |  \  |
//    |E4 \ |E6 \ |E8 \ |
//    2-----3-----4-----2
//      \E1 | \E2 | \E3 |
//       \  |  \  |  \  |
//        \ |   \ |   \ |
//          1     1     1
//
//  Reference Element T3:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  |  \
//    |  |   \
//    |  |    \
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
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_T3_ELEMENT[3*(2*NELEMX*(NELEMY-1)], 
//    the nodes that form each element.
//
{
  int base1;
  int base2;
  int element;
  int *element_node;
  int element_num;
  int element_order = 3;
  int i;
  int j;
  int ne;
  int nw;
  int se;
  int sw;

  element_num = sphere_grid_t3_element_num ( nelemx, nelemy );
  element_node = new int[element_order*element_num];

  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * nelemx + 2 - nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + i - 1;

      sw = base2;
      if ( i < nelemx )
      {
        se = base2 + 1;
      }
      else
      {
        se = base1;
      }
      nw = sw + nelemx;
      ne = se + nelemx;

      if ( j == 1 )
      {
        sw = 1;
        se = 1;
      }
      else if ( j == nelemx )
      {
        nw = base1 + nelemx;
        ne = base1 + nelemx;
      }

      if ( 1 < j )
      {
        element_node[0+element*element_order] = sw;
        element_node[1+element*element_order] = se;
        element_node[2+element*element_order] = nw;
        element = element + 1;
      }

      if ( j < nelemy )
      {
        element_node[0+element*element_order] = ne;
        element_node[1+element*element_order] = nw;
        element_node[2+element*element_order] = se;
        element = element + 1;
      }
    }
  }
  return element_node;
}
//****************************************************************************80

int sphere_grid_t3_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T3_ELEMENT_NUM counts the elements in a T3 sphere grid.
//
//  Example:
//
//    Input:
//
//      NELEMX = 6, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, int SPHERE_GRID_T3_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * ( nelemy - 1 );

  return element_num;
}
//****************************************************************************80

int sphere_grid_t3_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T3_NODE_NUM counts nodes in a T3 sphere grid.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_T3_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = nelemx * ( nelemy - 1 ) + 2;

  return node_num;
}
//****************************************************************************80

double *sphere_grid_t3_node_xyz ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T3_NODE_XYZ produces node coordinates for a T3 sphere grid.
//
//  Discussion:
//
//    The number of nodes to be generated is
//
//      NODE_NUM = NELEMX * ( NELEMY - 1 ) + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, double SPHERE_GRID_T3_NODE_XYZ[3*NODE_NUM], 
//    the node coordinates.
//
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t3_node_num ( nelemx, nelemy );
  node_xyz = new double[3*node_num];

  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( nelemy );

    for ( i = 1; i <= nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }
  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
}
//****************************************************************************80

int *sphere_grid_t6_element ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T6_ELEMENT produces a T6 sphere grid.
//
//  Discussion:
//
//    This would be the same as the grid in a plane, except that all the
//    nodes along the bottom edge are identified (replaced by a single node
//    that is the south pole) and similarly for the top edge, and the
//    nodes on the extreme right edge are identified pairwise with those 
//    on the extreme left edge.
//
//  Example:
//
//    Input:
//
//      NELEMX = 3, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NODE =
//        10,  8,  1,  9,  3,  4;
//        12, 10,  1, 11,  5,  6;
//         8, 12,  1, 13,  7,  2;
//         8, 10, 20,  9, 15, 14;
//        22, 20, 10, 21, 15, 16;
//        10, 12, 22, 11, 17, 16;
//        24, 22, 12, 23, 17, 18;
//        12,  8, 24, 13, 19, 18;
//        20, 24,  8, 25, 19, 14;
//        20, 22, 32, 21, 27, 26;
//        34, 32, 22, 33, 27, 28;
//        22, 24, 34, 23, 29, 28;
//        36, 34, 24, 35, 29, 30;
//        24, 20, 36, 25, 31, 30;
//        32, 36, 20, 37, 31, 26;
//        32, 34, 44, 33, 39, 38;
//        34, 36, 44, 35, 41, 40;
//        36, 32, 44, 37, 43, 42;
//
//  Grid:
//
//   44    44    44
//    |\    |\    |\
//   38 39 40 41 42 43 
//    |    \|    \|    \
//   32-33-34-35-36-37-32
//    |\    |\    |\    |
//   26 27 28 29 30 31 26
//    |    \|    \|    \|
//   20-21-22-23-24-25-20
//    |\    |\    |\    |
//   14 15 16 17 18 19 14
//    |    \|    \|    \|
//    8--9-10-11-12-13--8
//     \    |\    |\    |
//       3  4  5  6  7  2
//         \|    \|    \|
//          1     1     1
//
//  Reference Element T6:
//
//    |
//    1  3
//    |  |\
//    |  | \
//    S  6  5
//    |  |   \
//    |  |    \
//    0  1--4--2
//    |
//    +--0--R--1-->
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_T6_ELEMENT[6*2*NELEMX*(NELEMY-1)], 
//    the nodes that form each element.
//
{
  int base1;
  int base2;
  int c;
  int e;
  int element;
  int *element_node;
  int element_num;
  int element_order = 6;
  int i;
  int j;
  int n;
  int ne;
  int nw;
  int s;
  int se;
  int sw;
  int w;

  element_num = sphere_grid_t6_element_num ( nelemx, nelemy );
  element_node = new int[element_order*element_num];
  element = 0;

  for ( j = 1; j <= nelemy; j++ )
  {
    base1 = ( j - 1 ) * 2 * ( 2 * nelemx ) + 2 - 2 * nelemx;

    for ( i = 1; i <= nelemx; i++ )
    {
      base2 = base1 + 2 * ( i - 1 );

      sw = base2;
      s = base2 + 1;
      if ( i < nelemx )
      {
        se = base2 + 2;
      }
      else
      {
        se = base1;
      }

      w = sw + 2 * nelemx;
      c = s  + 2 * nelemx;
      e = se + 2 * nelemx;

      nw = w + 2 * nelemx;
      n  = c + 2 * nelemx;
      ne = e + 2 * nelemx;

      if ( j == 1 )
      {
        sw = 1;
        s  = 1;
        se = 1;
      }
      else if ( j == nelemy )
      {
        nw = base1 + 4 * nelemx;
        n  = base1 + 4 * nelemx;
        ne = base1 + 4 * nelemx;
      }

      if ( 1 < j )
      {
        element_node[0+element*element_order] = sw;
        element_node[1+element*element_order] = se;
        element_node[2+element*element_order] = nw;
        element_node[3+element*element_order] = s;
        element_node[4+element*element_order] = c;
        element_node[5+element*element_order] = w;
        element = element + 1;
      }

      if ( j < nelemy )
      {
        element_node[0+element*element_order] = ne;
        element_node[1+element*element_order] = nw;
        element_node[2+element*element_order] = se;
        element_node[3+element*element_order] = n;
        element_node[4+element*element_order] = c;
        element_node[5+element*element_order] = e;
        element = element + 1;
      }
    }
  }
  return element_node;
}
//****************************************************************************80

int sphere_grid_t6_element_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T6_ELEMENT_NUM counts the elements in a T6 sphere grid.
//
//  Example:
//
//    Input:
//
//      NELEMX = 6, NELEMY = 4
//
//    Output:
//
//      ELEMENT_NUM = 2 * NELEMX * ( NELEMY - 1 ) = 36
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions. 
//
//    Output, int SPHERE_GRID_T6_ELEMENT_NUM, the number of elements in the grid.
//
{
  int element_num;

  element_num = 2 * nelemx * ( nelemy - 1 );

  return element_num;
}
//****************************************************************************80

int sphere_grid_t6_node_num ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T6_NODE_NUM counts nodes in a T6 sphere grid.
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
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.
//
//    Output, int SPHERE_GRID_T6_NODE_NUM, the number of nodes in the grid.
//
{
  int node_num;

  node_num = 4 * nelemx * nelemy - 2 * nelemx + 2;

  return node_num;
}
//****************************************************************************80

double *sphere_grid_t6_node_xyz ( int nelemx, int nelemy )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_GRID_T6_NODE_XYZ produces node coordinates for a T6 sphere grid.
//
//  Discussion:
//
//    The number of nodes to be generated is
//
//      NODE_NUM = 4 * NELEMX * NELEMY - 2 * NELEMX + 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NELEMX, NELEMY, the number of elements along the
//    X and Y directions.  
//
//    Output, double SPHERE_GRID_T6_NODE_XYZ[3*NODE_NUM], 
//    the node coordinates.
//
{
  int i;
  int j;
  int node;
  int node_num;
  double *node_xyz;
  double phi;
  double pi = 3.141592653589793;
  double theta;

  node_num = sphere_grid_t6_node_num ( nelemx, nelemy );
  node_xyz = new double[3*node_num];
  node = 0;

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] = -1.0;
  node = node + 1;

  for ( j = 2 * nelemy; 2 <= j; j-- )
  {
    phi = ( double ) ( j - 1 ) * pi / ( double ) ( 2 * nelemy );

    for ( i = 1; i <= 2 * nelemx; i++ )
    {
      theta = ( double ) ( i - 1 ) * 2.0 * pi / ( double ) ( 2 * nelemx );

      node_xyz[0+node*3] = cos ( theta ) * sin ( phi );
      node_xyz[1+node*3] = sin ( theta ) * sin ( phi );
      node_xyz[2+node*3] =                 cos ( phi );
      node = node + 1;
    }
  }

  node_xyz[0+node*3] =  0.0;
  node_xyz[1+node*3] =  0.0;
  node_xyz[2+node*3] =  1.0;
  node = node + 1;

  return node_xyz;
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
//****************************************************************************80

void triangle_unit_set ( int rule, double xtab[], double ytab[], 
  double weight[] )

//****************************************************************************80
//
//  Purpose: 
//
//    TRIANGLE_UNIT_SET sets a quadrature rule in a unit triangle.
//
//  Integration region:
//
//    Points (X,Y) such that
//
//      0 <= X,
//      0 <= Y, and
//      X + Y <= 1.
//
//  Graph:
//
//      ^
//    1 | *
//      | |\
//    Y | | \
//      | |  \
//    0 | *---*
//      +------->
//        0 X 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  References:
//
//    H R Schwarz,
//    Methode der Finiten Elemente,
//    Teubner Studienbuecher, 1980.
//
//    Strang and Fix,
//    An Analysis of the Finite Element Method,
//    Prentice Hall, 1973, page 184.
//
//    Arthur H Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971.
//
//    O C Zienkiewicz,
//    The Finite Element Method,
//    McGraw Hill, Third Edition, 1977, page 201.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//     1, NORDER =  1, precision 1, Zienkiewicz #1.
//     2, NORDER =  3, precision 1, the "vertex rule".
//     3, NORDER =  3, precision 2, Strang and Fix formula #1.
//     4, NORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
//     5, NORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
//     6, NORDER =  6, precision 3, Strang and Fix formula #4.
//     7, NORDER =  6, precision 3, Stroud formula T2:3-1.
//     8, NORDER =  6, precision 4, Strang and Fix formula #5.
//     9, NORDER =  7, precision 4, Strang and Fix formula #6.
//    10, NORDER =  7, precision 5, Strang and Fix formula #7,
//        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
//    11, NORDER =  9, precision 6, Strang and Fix formula #8.
//    12, NORDER = 12, precision 6, Strang and Fix formula #9.
//    13, ORDER = 13, precision 7, Strang and Fix formula #10.
//        Note that there is a typographical error in Strang and Fix
//        which lists the value of the XSI(3) component of the
//        last generator point as 0.4869... when it should be 0.04869...
//    14, ORDER =  7, precision ?.
//    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//    16, ORDER = 64, precision 15, triangular product Gauss rule.
//    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
//    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
//    20, ORDER = 37, precision 13, from ACM TOMS #706.
//
//    Output, double XTAB[NORDER], YTAB[NORDER], the abscissas.
//
//    Output, double WEIGHT[NORDER], the weights of the rule.
//
{
  double a;
  double b;
  double c;
  double d;
  double e;
  double f;
  double g;
  double h;
  int i;
  int j;
  int k;
  int order2;
  double p;
  double q;
  double r;
  double s;
  double t;
  double u;
  double v;
  double w;
  double w1;
  double w2;
  double w3;
  double w4;
  double w5;
  double w6;
  double w7;
  double w8;
  double w9;
  double weight1[8];
  double weight2[8];
  double wx;
  double x;
  double xtab1[8];
  double xtab2[8];
  double y;
  double z;
//
//  1 point, precision 1.
//
  if ( rule == 1 )
  {
    xtab[0] = 1.0 / 3.0;
    ytab[0] = 1.0 / 3.0;
    weight[0] = 1.0;
  }
//
//  3 points, precision 1, the "vertex rule".
//
  else if ( rule == 2 )
  {
    xtab[0] = 1.0;
    xtab[1] = 0.0;
    xtab[2] = 0.0;

    ytab[0] = 0.0;
    ytab[1] = 1.0;
    ytab[2] = 0.0;

    weight[0] = 1.0 / 3.0;
    weight[1] = 1.0 / 3.0;
    weight[2] = 1.0 / 3.0;
  }
//
//  3 points, precision 2, Strang and Fix formula #1.
//
  else if ( rule == 3 )
  {
    xtab[0] = 4.0 / 6.0;
    xtab[1] = 1.0 / 6.0;
    xtab[2] = 1.0 / 6.0;

    ytab[0] = 1.0 / 6.0;
    ytab[1] = 4.0 / 6.0;
    ytab[2] = 1.0 / 6.0;

    weight[0] = 1.0 / 3.0;
    weight[1] = 1.0 / 3.0;
    weight[2] = 1.0 / 3.0;
  }
//
//  3 points, precision 2, Strang and Fix formula #2.
//
  else if ( rule == 4 )
  {
    xtab[0] = 0.0;
    xtab[1] = 1.0 / 2.0;
    xtab[2] = 1.0 / 2.0;

    ytab[0] = 1.0 / 2.0;
    ytab[1] = 0.0;
    ytab[2] = 1.0 / 2.0;

    weight[0] = 1.0 / 3.0;
    weight[1] = 1.0 / 3.0;
    weight[2] = 1.0 / 3.0;
  }
//
//  4 points, precision 3, Strang and Fix formula #3.
//
  else if ( rule == 5 )
  {
    xtab[0] = 10.0 / 30.0;
    xtab[1] = 18.0 / 30.0;
    xtab[2] =  6.0 / 30.0;
    xtab[3] =  6.0 / 30.0;

    ytab[0] = 10.0 / 30.0;
    ytab[1] =  6.0 / 30.0;
    ytab[2] = 18.0 / 30.0;
    ytab[3] =  6.0 / 30.0;

    weight[0] = -27.0 / 48.0;
    weight[1] =  25.0 / 48.0;
    weight[2] =  25.0 / 48.0;
    weight[3] =  25.0 / 48.0;
  }
//
//  6 points, precision 3, Strang and Fix formula #4.
//
  else if ( rule == 6 )
  {
    xtab[0] = 0.659027622374092;
    xtab[1] = 0.659027622374092;
    xtab[2] = 0.231933368553031;
    xtab[3] = 0.231933368553031;
    xtab[4] = 0.109039009072877;
    xtab[5] = 0.109039009072877;

    ytab[0] = 0.231933368553031;
    ytab[1] = 0.109039009072877;
    ytab[2] = 0.659027622374092;
    ytab[3] = 0.109039009072877;
    ytab[4] = 0.659027622374092;
    ytab[5] = 0.231933368553031;

    weight[0] = 1.0 / 6.0;
    weight[1] = 1.0 / 6.0;
    weight[2] = 1.0 / 6.0;
    weight[3] = 1.0 / 6.0;
    weight[4] = 1.0 / 6.0;
    weight[5] = 1.0 / 6.0;
  }
//
//  6 points, precision 3, Stroud T2:3-1.
//
  else if ( rule == 7 )
  {
    xtab[0] = 0.0;
    xtab[1] = 0.5;
    xtab[2] = 0.5;
    xtab[3] = 2.0 /  3.0;
    xtab[4] = 1.0 /  6.0;
    xtab[5] = 1.0 /  6.0;

    ytab[0] = 0.5;
    ytab[1] = 0.0;
    ytab[2] = 0.5;
    ytab[3] = 1.0 /  6.0;
    ytab[4] = 2.0 /  3.0;
    ytab[5] = 1.0 /  6.0;

    weight[0] = 1.0 / 30.0;
    weight[1] = 1.0 / 30.0;
    weight[2] = 1.0 / 30.0;
    weight[3] = 3.0 / 10.0;
    weight[4] = 3.0 / 10.0;
    weight[5] = 3.0 / 10.0;
  }
//
//  6 points, precision 4, Strang and Fix, formula #5.
//
  else if ( rule == 8 )
  {
    xtab[0] = 0.816847572980459;
    xtab[1] = 0.091576213509771;
    xtab[2] = 0.091576213509771;
    xtab[3] = 0.108103018168070;
    xtab[4] = 0.445948490915965;
    xtab[5] = 0.445948490915965;

    ytab[0] = 0.091576213509771;
    ytab[1] = 0.816847572980459;
    ytab[2] = 0.091576213509771;
    ytab[3] = 0.445948490915965;
    ytab[4] = 0.108103018168070;
    ytab[5] = 0.445948490915965;

    weight[0] = 0.109951743655322;
    weight[1] = 0.109951743655322;
    weight[2] = 0.109951743655322;
    weight[3] = 0.223381589678011;
    weight[4] = 0.223381589678011;
    weight[5] = 0.223381589678011;
  }
//
//  7 points, precision 4, Strang and Fix formula #6.
//
  else if ( rule == 9 )
  {
    xtab[0] = 1.0 / 3.0;
    xtab[1] = 0.736712498968435;
    xtab[2] = 0.736712498968435;
    xtab[3] = 0.237932366472434;
    xtab[4] = 0.237932366472434;
    xtab[5] = 0.025355134551932;
    xtab[6] = 0.025355134551932;

    ytab[0] = 1.0 / 3.0;
    ytab[1] = 0.237932366472434;
    ytab[2] = 0.025355134551932;
    ytab[3] = 0.736712498968435;
    ytab[4] = 0.025355134551932;
    ytab[5] = 0.736712498968435;
    ytab[6] = 0.237932366472434;

    weight[0] = 0.375;
    weight[1] = 0.1041666666666667;
    weight[2] = 0.1041666666666667;
    weight[3] = 0.1041666666666667;
    weight[4] = 0.1041666666666667;
    weight[5] = 0.1041666666666667;
    weight[6] = 0.1041666666666667;
  }
//
//  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1
//
  else if ( rule == 10 )
  {
    xtab[0] = 1.0 / 3.0;
    xtab[1] = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    xtab[2] = ( 6.0 -           sqrt ( 15.0 ) ) / 21.0;
    xtab[3] = ( 6.0 -           sqrt ( 15.0 ) ) / 21.0;
    xtab[4] = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    xtab[5] = ( 6.0 +           sqrt ( 15.0 ) ) / 21.0;
    xtab[6] = ( 6.0 +           sqrt ( 15.0 ) ) / 21.0;

    ytab[0] = 1.0 / 3.0;
    ytab[1] = ( 6.0 -           sqrt ( 15.0 ) ) / 21.0;
    ytab[2] = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    ytab[3] = ( 6.0 -           sqrt ( 15.0 ) ) / 21.0;
    ytab[4] = ( 6.0 +           sqrt ( 15.0 ) ) / 21.0;
    ytab[5] = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    ytab[6] = ( 6.0 +           sqrt ( 15.0 ) ) / 21.0;

    weight[0] = 0.225;
    weight[1] = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    weight[2] = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    weight[3] = ( 155.0 - sqrt ( 15.0 ) ) / 1200.0;
    weight[4] = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;
    weight[5] = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;
    weight[6] = ( 155.0 + sqrt ( 15.0 ) ) / 1200.0;
  }
//
//  9 points, precision 6, Strang and Fix formula #8.
//
  else if ( rule == 11 )
  {
    xtab[0] = 0.124949503233232;
    xtab[1] = 0.437525248383384;
    xtab[2] = 0.437525248383384;
    xtab[3] = 0.797112651860071;
    xtab[4] = 0.797112651860071;
    xtab[5] = 0.165409927389841;
    xtab[6] = 0.165409927389841;
    xtab[7] = 0.037477420750088;
    xtab[8] = 0.037477420750088;

    ytab[0] = 0.437525248383384;
    ytab[1] = 0.124949503233232;
    ytab[2] = 0.437525248383384;
    ytab[3] = 0.165409927389841;
    ytab[4] = 0.037477420750088;
    ytab[5] = 0.797112651860071;
    ytab[6] = 0.037477420750088;
    ytab[7] = 0.797112651860071;
    ytab[8] = 0.165409927389841;

    weight[0] = 0.205950504760887;
    weight[1] = 0.205950504760887;
    weight[2] = 0.205950504760887;
    weight[3] = 0.063691414286223;
    weight[4] = 0.063691414286223;
    weight[5] = 0.063691414286223;
    weight[6] = 0.063691414286223;
    weight[7] = 0.063691414286223;
    weight[8] = 0.063691414286223;
  }
//
//  12 points, precision 6, Strang and Fix, formula #9.
//
  else if ( rule == 12 )
  {
    xtab[0] = 0.873821971016996;
    xtab[1] = 0.063089014491502;
    xtab[2] = 0.063089014491502;
    xtab[3] = 0.249286745170910;
    xtab[4] = 0.501426509658179;
    xtab[5] = 0.249286745170910;
    xtab[6] = 0.636502499121399;
    xtab[7] = 0.636502499121399;
    xtab[8] = 0.310352451033785;
    xtab[9] = 0.310352451033785;
    xtab[10] = 0.053145049844816;
    xtab[11] = 0.053145049844816;

    ytab[0] = 0.063089014491502;
    ytab[1] = 0.873821971016996;
    ytab[2] = 0.063089014491502;
    ytab[3] = 0.501426509658179;
    ytab[4] = 0.249286745170910;
    ytab[5] = 0.249286745170910;
    ytab[6] = 0.310352451033785;
    ytab[7] = 0.053145049844816;
    ytab[8] = 0.636502499121399;
    ytab[9] = 0.053145049844816;
    ytab[10] = 0.636502499121399;
    ytab[11] = 0.310352451033785;

    weight[0] = 0.050844906370207;
    weight[1] = 0.050844906370207;
    weight[2] = 0.050844906370207;
    weight[3] = 0.116786275726379;
    weight[4] = 0.116786275726379;
    weight[5] = 0.116786275726379;
    weight[6] = 0.082851075618374;
    weight[7] = 0.082851075618374;
    weight[8] = 0.082851075618374;
    weight[9] = 0.082851075618374;
    weight[10] = 0.082851075618374;
    weight[11] = 0.082851075618374;
  }
//
//  13 points, precision 7, Strang and Fix, formula #10.
//
  else if ( rule == 13 )
  {
    xtab[0] = 0.479308067841923;
    xtab[1] = 0.260345966079038;
    xtab[2] = 0.260345966079038;
    xtab[3] = 0.869739794195568;
    xtab[4] = 0.065130102902216;
    xtab[5] = 0.065130102902216;
    xtab[6] = 0.638444188569809;
    xtab[7] = 0.638444188569809;
    xtab[8] = 0.312865496004875;
    xtab[9] = 0.312865496004875;
    xtab[10] = 0.048690315425316;
    xtab[11] = 0.048690315425316;
    xtab[12] = 1.0 / 3.0;

    ytab[0] = 0.260345966079038;
    ytab[1] = 0.479308067841923;
    ytab[2] = 0.260345966079038;
    ytab[3] = 0.065130102902216;
    ytab[4] = 0.869739794195568;
    ytab[5] = 0.065130102902216;
    ytab[6] = 0.312865496004875;
    ytab[7] = 0.048690315425316;
    ytab[8] = 0.638444188569809;
    ytab[9] = 0.048690315425316;
    ytab[10] = 0.638444188569809;
    ytab[11] = 0.312865496004875;
    ytab[12] = 1.0 / 3.0;

    weight[0] = 0.175615257433204;
    weight[1] = 0.175615257433204;
    weight[2] = 0.175615257433204;
    weight[3] = 0.053347235608839;
    weight[4] = 0.053347235608839;
    weight[5] = 0.053347235608839;
    weight[6] = 0.077113760890257;
    weight[7] = 0.077113760890257;
    weight[8] = 0.077113760890257;
    weight[9] = 0.077113760890257;
    weight[10] = 0.077113760890257;
    weight[11] = 0.077113760890257;
    weight[12] = -0.149570044467670;
  }
//
//  7 points, precision ?.
//
  else if ( rule == 14 )
  {
    a = 1.0 / 3.0;
    b = 1.0;
    c = 0.5;
    z = 0.0;

    u = 27.0 / 60.0;
    v =  3.0 / 60.0;
    w =  8.0 / 60.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = z;
    xtab[3] = z;
    xtab[4] = z;
    xtab[5] = c;
    xtab[6] = c;

    ytab[0] = a;
    ytab[1] = z;
    ytab[2] = b;
    ytab[3] = z;
    ytab[4] = c;
    ytab[5] = z;
    ytab[6] = c;

    weight[0] = u;
    weight[1] = v;
    weight[2] = v;
    weight[3] = v;
    weight[4] = w;
    weight[5] = w;
    weight[6] = w;
    weight[7] = w;
  }
//
//  16 points.
//
  else if ( rule == 15 )
  {
//
//  Legendre rule of order 4.
//
    order2 = 4;

    xtab1[0] = -0.861136311594052575223946488893;
    xtab1[1] = -0.339981043584856264802665759103;
    xtab1[2] =  0.339981043584856264802665759103;
    xtab1[3] =  0.861136311594052575223946488893;

    weight1[0] = 0.347854845137453857373063949222;
    weight1[1] = 0.652145154862546142626936050778;
    weight1[2] = 0.652145154862546142626936050778;
    weight1[3] = 0.347854845137453857373063949222;

    for ( i = 0; i < order2; i++ )
    {
      xtab1[i] = 0.5 * ( xtab1[i] + 1.0 );
    }

    weight2[0] = 0.1355069134;
    weight2[1] = 0.2034645680;
    weight2[2] = 0.1298475476;
    weight2[3] = 0.0311809709;

    xtab2[0] = 0.0571041961;
    xtab2[1] = 0.2768430136;
    xtab2[2] = 0.5835904324;
    xtab2[3] = 0.8602401357;

    k = 0;
    for ( i = 0; i < order2; i++ )
    {
      for ( j = 0; j < order2; j++ )
      {
        xtab[k] = xtab2[j];
        ytab[k] = xtab1[i] * ( 1.0 - xtab2[j] );
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
//
//  64 points, precision 15.
//
  else if ( rule == 16 )
  {
//
//  Legendre rule of order 8.
//
    order2 = 8;

    xtab1[0] = -0.960289856497536231683560868569;
    xtab1[1] = -0.796666477413626739591553936476;
    xtab1[2] = -0.525532409916328985817739049189;
    xtab1[3] = -0.183434642495649804939476142360;
    xtab1[4] =  0.183434642495649804939476142360;
    xtab1[5] =  0.525532409916328985817739049189;
    xtab1[6] =  0.796666477413626739591553936476;
    xtab1[7] =  0.960289856497536231683560868569;

    weight1[0] = 0.101228536290376259152531354310;
    weight1[1] = 0.222381034453374470544355994426;
    weight1[2] = 0.313706645877887287337962201987;
    weight1[3] = 0.362683783378361982965150449277;
    weight1[4] = 0.362683783378361982965150449277;
    weight1[5] = 0.313706645877887287337962201987;
    weight1[6] = 0.222381034453374470544355994426;
    weight1[7] = 0.101228536290376259152531354310;

    weight2[0] = 0.00329519144;
    weight2[1] = 0.01784290266;
    weight2[2] = 0.04543931950;
    weight2[3] = 0.07919959949;
    weight2[4] = 0.10604735944;
    weight2[5] = 0.11250579947;
    weight2[6] = 0.09111902364;
    weight2[7] = 0.04455080436;

    xtab2[0] = 0.04463395529;
    xtab2[1] = 0.14436625704;
    xtab2[2] = 0.28682475714;
    xtab2[3] = 0.45481331520;
    xtab2[4] = 0.62806783542;
    xtab2[5] = 0.78569152060;
    xtab2[6] = 0.90867639210;
    xtab2[7] = 0.98222008485;

    k = 0;
    for ( j = 0; j < order2; j++ )
    {
      for ( i = 0; i < order2; i++ )
      {
        xtab[k] = 1.0 - xtab2[j];
        ytab[k] = 0.5 * ( 1.0 + xtab1[i] ) * xtab2[j];
        weight[k] = weight1[i] * weight2[j];
        k = k + 1;
      }
    }
  }
//
//  19 points, precision 8.
//
  else if ( rule == 17 )
  {
    a = 1.0 / 3.0;
    b = ( 9.0 + 2.0 * sqrt ( 15.0 ) ) / 21.0;
    c = ( 6.0 -       sqrt ( 15.0 ) ) / 21.0;
    d = ( 9.0 - 2.0 * sqrt ( 15.0 ) ) / 21.0;
    e = ( 6.0 +       sqrt ( 15.0 ) ) / 21.0;
    f = ( 40.0 - 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) + 2.0 * sqrt ( 105.0 ) ) / 90.0;
    g = ( 25.0 +  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) - sqrt ( 105.0 ) ) / 90.0;
    p = ( 40.0 + 10.0 * sqrt ( 15.0 ) 
      + 10.0 * sqrt ( 7.0 ) - 2.0 * sqrt ( 105.0 ) ) / 90.0;
    q = ( 25.0 -  5.0 * sqrt ( 15.0 ) 
      -  5.0 * sqrt ( 7.0 ) + sqrt ( 105.0 ) ) / 90.0;
    r = ( 40.0 + 10.0 * sqrt ( 7.0 ) ) / 90.0;
    s = ( 25.0 +  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      - sqrt ( 105.0 ) ) / 90.0;
    t = ( 25.0 -  5.0 * sqrt ( 15.0 ) - 5.0 * sqrt ( 7.0 ) 
      + sqrt ( 105.0 ) ) / 90.0;

    w1 = ( 7137.0 - 1800.0 * sqrt ( 7.0 ) ) / 62720.0;
    w2 = -9301697.0 / 4695040.0 - 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      + 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w2 = w2 / 3.0;
    w3 = -9301697.0 / 4695040.0 + 13517313.0 * sqrt ( 15.0 ) 
      / 23475200.0 
      + 764885.0 * sqrt ( 7.0 ) / 939008.0 
      - 198763.0 * sqrt ( 105.0 ) / 939008.0;
    w3 = w3 / 3.0;
    w4 = ( 102791225.0 - 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      + 9914825.0 * sqrt ( 105.0 ) ) / 59157504.0;
    w4 = w4 / 3.0;
    w5 = ( 102791225.0 + 23876225.0 * sqrt ( 15.0 ) 
      - 34500875.0 * sqrt ( 7.0 ) 
      - 9914825 * sqrt ( 105.0 ) ) / 59157504.0;
    w5 = w5 / 3.0;
    w6 = ( 11075.0 - 3500.0 * sqrt ( 7.0 ) ) / 8064.0;
    w6 = w6 / 6.0;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = r;
    xtab[15] = s;
    xtab[16] = s;
    xtab[17] = t;
    xtab[18] = t;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = t;
    ytab[15] = r;
    ytab[16] = t;
    ytab[17] = r;
    ytab[18] = s;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
//
//  19 points, precision 9.
//
  else if ( rule == 18 )
  {
    a = 1.0 / 3.0;
    b = 0.02063496160252593;
    c = 0.4896825191987370;
    d = 0.1258208170141290;
    e = 0.4370895914929355;
    f = 0.6235929287619356;
    g = 0.1882035356190322;
    r = 0.9105409732110941;
    s = 0.04472951339445297;
    t = 0.7411985987844980;
    u = 0.03683841205473626;
    v = 0.22196288916076574;

    w1 = 0.09713579628279610;
    w2 = 0.03133470022713983;
    w3 = 0.07782754100477543;
    w4 = 0.07964773892720910;
    w5 = 0.02557767565869810;
    w6 = 0.04328353937728940;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = r;
    xtab[11] = s;
    xtab[12] = s;
    xtab[13] = t;
    xtab[14] = t;
    xtab[15] = u;
    xtab[16] = u;
    xtab[17] = v;
    xtab[18] = v;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = s;
    ytab[11] = r;
    ytab[12] = s;
    ytab[13] = u;
    ytab[14] = v;
    ytab[15] = t;
    ytab[16] = v;
    ytab[17] = t;
    ytab[18] = u;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w6;
    weight[17] = w6;
    weight[18] = w6;
  }
//
//  28 points, precision 11.
//
  else if ( rule == 19 )
  {
    a = 1.0 / 3.0;
    b = 0.9480217181434233;
    c = 0.02598914092828833;
    d = 0.8114249947041546;
    e = 0.09428750264792270;
    f = 0.01072644996557060;
    g = 0.4946367750172147;
    p = 0.5853132347709715;
    q = 0.2073433826145142;
    r = 0.1221843885990187;
    s = 0.4389078057004907;
    t = 0.6779376548825902;
    u = 0.04484167758913055;
    v = 0.27722066752827925;
    w = 0.8588702812826364;
    x = 0.0;
    y = 0.1411297187173636;

    w1 = 0.08797730116222190;
    w2 = 0.008744311553736190;
    w3 = 0.03808157199393533;
    w4 = 0.01885544805613125;
    w5 = 0.07215969754474100;
    w6 = 0.06932913870553720;
    w7 = 0.04105631542928860;
    w8 = 0.007362383783300573;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;
    xtab[10] = p;
    xtab[11] = q;
    xtab[12] = q;
    xtab[13] = r;
    xtab[14] = s;
    xtab[15] = s;
    xtab[16] = t;
    xtab[17] = t;
    xtab[18] = u;
    xtab[19] = u;
    xtab[20] = v;
    xtab[21] = v;
    xtab[22] = w;
    xtab[23] = w;
    xtab[24] = x;
    xtab[25] = x;
    xtab[26] = y;
    xtab[27] = y;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;
    ytab[10] = q;
    ytab[11] = p;
    ytab[12] = q;
    ytab[13] = s;
    ytab[14] = r;
    ytab[15] = s;
    ytab[16] = u;
    ytab[17] = v;
    ytab[18] = t;
    ytab[19] = v;
    ytab[20] = t;
    ytab[21] = u;
    ytab[22] = x;
    ytab[23] = y;
    ytab[24] = w;
    ytab[25] = y;
    ytab[26] = w;
    ytab[27] = x;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w7;
    weight[20] = w7;
    weight[21] = w7;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w8;
    weight[26] = w8;
    weight[27] = w8;
  }
//
//  37 points, precision 13.
//
  else if ( rule == 20 )
  {
    a = 1.0 / 3.0;
    b = 0.950275662924105565450352089520;
    c = 0.024862168537947217274823955239;
    d = 0.171614914923835347556304795551;
    e = 0.414192542538082326221847602214;
    f = 0.539412243677190440263092985511;
    g = 0.230293878161404779868453507244;

    w1 = 0.051739766065744133555179145422;
    w2 = 0.008007799555564801597804123460;
    w3 = 0.046868898981821644823226732071;
    w4 = 0.046590940183976487960361770070;
    w5 = 0.031016943313796381407646220131;
    w6 = 0.010791612736631273623178240136;
    w7 = 0.032195534242431618819414482205;
    w8 = 0.015445834210701583817692900053;
    w9 = 0.017822989923178661888748319485;
    wx = 0.037038683681384627918546472190;

    xtab[0] = a;
    xtab[1] = b;
    xtab[2] = c;
    xtab[3] = c;
    xtab[4] = d;
    xtab[5] = e;
    xtab[6] = e;
    xtab[7] = f;
    xtab[8] = g;
    xtab[9] = g;

    ytab[0] = a;
    ytab[1] = c;
    ytab[2] = b;
    ytab[3] = c;
    ytab[4] = e;
    ytab[5] = d;
    ytab[6] = e;
    ytab[7] = g;
    ytab[8] = f;
    ytab[9] = g;

    weight[0] = w1;
    weight[1] = w2;
    weight[2] = w2;
    weight[3] = w2;
    weight[4] = w3;
    weight[5] = w3;
    weight[6] = w3;
    weight[7] = w4;
    weight[8] = w4;
    weight[9] = w4;
    weight[10] = w5;
    weight[11] = w5;
    weight[12] = w5;
    weight[13] = w6;
    weight[14] = w6;
    weight[15] = w6;
    weight[16] = w7;
    weight[17] = w7;
    weight[18] = w7;
    weight[19] = w8;
    weight[20] = w8;
    weight[21] = w8;
    weight[22] = w8;
    weight[23] = w8;
    weight[24] = w8;
    weight[25] = w9;
    weight[26] = w9;
    weight[27] = w9;
    weight[28] = w9;
    weight[29] = w9;
    weight[30] = w9;
    weight[31] = wx;
    weight[32] = wx;
    weight[33] = wx;
    weight[34] = wx;
    weight[35] = wx;
    weight[36] = wx;

    a = 0.772160036676532561750285570113;
    b = 0.113919981661733719124857214943;

    xtab[10] = a;
    ytab[10] = b;

    xtab[11] = b;
    ytab[11] = a;

    xtab[12] = b;
    ytab[12] = b;

    a = 0.009085399949835353883572964740;
    b = 0.495457300025082323058213517632;

    xtab[13] = a;
    ytab[13] = b;

    xtab[14] = b;
    ytab[14] = a;

    xtab[15] = b;
    ytab[15] = b;

    a = 0.062277290305886993497083640527;
    b = 0.468861354847056503251458179727;

    xtab[16] = a;
    ytab[16] = b;

    xtab[17] = b;
    ytab[17] = a;

    xtab[18] = b;
    ytab[18] = b;

    a = 0.022076289653624405142446876931;
    b = 0.851306504174348550389457672223;
    c = 1.0 - a - b;

    xtab[19] = a;
    ytab[19] = b;

    xtab[20] = a;
    ytab[20] = c;

    xtab[21] = b;
    ytab[21] = a;

    xtab[22] = b;
    ytab[22] = c;

    xtab[23] = c;
    ytab[23] = a;

    xtab[24] = c;
    ytab[24] = b;

    a = 0.018620522802520968955913511549;
    b = 0.689441970728591295496647976487;
    c = 1.0 - a - b;

    xtab[25] = a;
    ytab[25] = b;

    xtab[26] = a;
    ytab[26] = c;

    xtab[27] = b;
    ytab[27] = a;

    xtab[28] = b;
    ytab[28] = c;

    xtab[29] = c;
    ytab[29] = a;

    xtab[30] = c;
    ytab[30] = b;

    a = 0.096506481292159228736516560903;
    b = 0.635867859433872768286976979827;
    c = 1.0 - a - b;

    xtab[31] = a;
    ytab[31] = b;

    xtab[32] = a;
    ytab[32] = c;

    xtab[33]  = b;
    ytab[33] = a;

    xtab[34] = b;
    ytab[34] = c;

    xtab[35] = c;
    ytab[35] = a;

    xtab[36] = c;
    ytab[36] = b;
  }
  else
  {
    cout << "\n";
    cout << "TRIANGLE_UNIT_SET - Fatal error!\n";
    cout << "  Illegal value of RULE = " << rule << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int triangle_unit_size ( int rule )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_SIZE returns the "size" of a unit triangle quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gilbert Strang, George Fix,
//    An Analysis of the Finite Element Method,
//    Prentice Hall, 1973,
//    TA335.S77.
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    McGraw Hill, Third Edition, 1977, page 202.
//
//  Parameters:
//
//    Input, int RULE, the index of the rule.
//
//     1, ORDER =  1, precision 1, Zienkiewicz #1.
//     2, ORDER =  3, precision 1, the "vertex rule".
//     3, ORDER =  3, precision 2, Strang and Fix formula #1.
//     4, ORDER =  3, precision 2, Strang and Fix formula #2, Zienkiewicz #2.
//     5, ORDER =  4, precision 3, Strang and Fix formula #3, Zienkiewicz #3.
//     6, ORDER =  6, precision 3, Strang and Fix formula #4.
//     7, ORDER =  6, precision 3, Stroud formula T2:3-1.
//     8, ORDER =  6, precision 4, Strang and Fix formula #5.
//     9, ORDER =  7, precision 4, Strang and Fix formula #6.
//    10, ORDER =  7, precision 5, Strang and Fix formula #7,
//        Stroud formula T2:5-1, Zienkiewicz #4, Schwarz Table 2.2.
//    11, ORDER =  9, precision 6, Strang and Fix formula #8.
//    12, ORDER = 12, precision 6, Strang and Fix formula #9.
//    13, ORDER = 13, precision 7, Strang and Fix formula #10.
//    14, ORDER =  7, precision ?.
//    15, ORDER = 16, precision 7, conical product Gauss, Stroud formula T2:7-1.
//    16, ORDER = 64, precision 15, triangular product Gauss rule.
//    17, ORDER = 19, precision 8, from CUBTRI, ACM TOMS #584.
//    18, ORDER = 19, precision 9, from TRIEX, Lyness and Jespersen.
//    19, ORDER = 28, precision 11, from TRIEX, Lyness and Jespersen.
//    20, ORDER = 37, precision 13, from ACM TOMS #706.
//
//    Output, int TRIANGLE_UNIT_SIZE, the order of the rule.
//
{
  int value;

  if ( rule == 1 )
  {
    value = 1;
  }
  else if ( rule == 2 )
  {
    value = 3;
  }
  else if ( rule == 3 )
  {
    value = 3;
  }
  else if ( rule == 4 )
  {
    value = 3;
  }
  else if ( rule == 5 )
  {
    value = 4;
  }
  else if ( rule == 6 )
  {
    value = 6;
  }
  else if ( rule == 7 )
  {
    value = 6;
  }
  else if ( rule == 8 )
  {
    value = 6;
  }
  else if ( rule == 9 )
  {
    value = 7;
  }
  else if ( rule == 10 )
  {
    value = 7;
  }
  else if ( rule == 11 )
  {
    value = 9;
  }
  else if ( rule == 12 )
  {
    value = 12;
  }
  else if ( rule == 13 )
  {
    value = 13;
  }
  else if ( rule == 14 )
  {
    value = 7;
  }
  else if ( rule == 15 )
  {
    value = 16;
  }
  else if ( rule == 16 )
  {
    value = 64;
  }
  else if ( rule == 17 )
  {
    value = 19;
  }
  else if ( rule == 18 )
  {
    value = 19;
  }
  else if ( rule == 19 )
  {
    value = 28;
  }
  else if ( rule == 20 )
  {
    value = 37;
  }
  else
  {
    value = -1;
  }

  return value;
}
