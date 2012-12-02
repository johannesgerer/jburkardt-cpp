# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <fstream>
# include <cstring>

using namespace std;

# include "fem3d_pack.hpp"

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

double *basis_brick8 ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK8: BRICK8 basis functions at natural coordinates.
//
//  Discussion:
//
//      8------7        t  s
//     /|     /|        | /
//    5------6 |        |/
//    | |    | |        0-------r
//    | 4----|-3        
//    |/     |/        
//    1------2        
//                   
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], natural coordinates of evaluation
//    points.
//
//    Output, double BASIS_BRICK8[8*N], the basis function values.
//
{
  int j;
  double *phi;

  phi = new double[8*n];

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*8] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 8.0;
    phi[1+j*8] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 8.0;
    phi[2+j*8] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 8.0;
    phi[3+j*8] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 8.0;
    phi[4+j*8] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 8.0;
    phi[5+j*8] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 8.0;
    phi[6+j*8] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 8.0;
    phi[7+j*8] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 8.0;
  }

  return phi;
}
//****************************************************************************80

void basis_brick8_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK8_TEST verifies BASIS_BRICK8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
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
  int i;
  int j;
  int n;
  int node_num = 8;
  double *p;
  double *phi;
  double phi_sum;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "BASIS_BRICK8_TEST:\n";
  cout << "  Verify basis functions for element BRICK8.\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  n = node_num;

  p = nodes_brick8 ( );

  phi = basis_brick8 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    cout << "  ";
    for ( i = 0; i < node_num; i++ )
    {
      cout << setw(7) << phi[i+j*node_num];
    }
    cout << "\n";
  }

  delete [] p;
  delete [] phi;

  cout << "\n";
  cout << "  The basis function values at ANY point P\n";
  cout << "  should sum to 1:\n";
  cout << "\n";
  cout << "    ------------P-------------     PHI_SUM\n";
  cout << "\n";

  n = test_num;
  seed = 123456789;

  p = r8mat_uniform_01_new ( 3, n, &seed );

  phi = basis_brick8 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    phi_sum = r8vec_sum ( node_num, phi+j*node_num );
    cout << "  " << setw(8) << p[0+j*3]
         << "  " << setw(8) << p[1+j*3]
         << "  " << setw(8) << p[2+j*3]
         << "  " << setw(8) << phi_sum << "\n";
  }

  delete [] p;
  delete [] phi;

  return;
}
//****************************************************************************80

double *basis_brick20 ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
//
//  Discussion:
//
//        8----19---7
//       /|        /|
//     20 |      18 |        t   s
//     /  16     /  15       |  /
//    5----17---6   |        | /
//    |   |     |   |        |/
//    |   4--11-|---3        0---------r
//   13  /     14  /        
//    | 12      | 10       
//    |/        |/        
//    1----9----2
//                   
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], natural coordinates of evaluation
//    points.
//
//    Output, double BASIS_BRICK20[20*N], the basis function values.
//
{
  int j;
  double *phi;

  phi = new double[20*n];

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*20] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) 
      * ( - p[0+j*3] - p[1+j*3] - p[2+j*3] - 2.0 ) / 8.0;
    phi[1+j*20] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) 
      * ( + p[0+j*3] - p[1+j*3] - p[2+j*3] - 2.0 ) / 8.0;
    phi[2+j*20] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 - p[2+j*3] ) 
      * ( + p[0+j*3] + p[1+j*3] - p[2+j*3] - 2.0 ) / 8.0;
    phi[3+j*20] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 - p[2+j*3] ) 
      * ( - p[0+j*3] + p[1+j*3] - p[2+j*3] - 2.0 ) / 8.0;
    phi[4+j*20] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) 
      * ( - p[0+j*3] - p[1+j*3] + p[2+j*3] - 2.0 ) / 8.0;
    phi[5+j*20] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) 
      * ( + p[0+j*3] - p[1+j*3] + p[2+j*3] - 2.0 ) / 8.0;
    phi[6+j*20] = 
      ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 + p[2+j*3] ) 
      * ( + p[0+j*3] + p[1+j*3] + p[2+j*3] - 2.0 ) / 8.0;
    phi[7+j*20] = 
      ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) * ( 1.0 + p[2+j*3] ) 
      * ( - p[0+j*3] + p[1+j*3] + p[2+j*3] - 2.0 ) / 8.0;

    phi[8+j*20] =  ( 1.0 + p[0+j*3] ) * ( 1.0 - p[0+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[9+j*20] =  ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[10+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 - p[0+j*3] ) 
                 * ( 1.0 + p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[11+j*20] = ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[12+j*20] = ( 1.0 - p[0+j*3] ) * ( 1.0 - p[1+j*3] ) 
                 * ( 1.0 + p[2+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[13+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 - p[1+j*3] ) 
                 * ( 1.0 + p[2+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[14+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 + p[2+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[15+j*20] = ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 + p[2+j*3] ) * ( 1.0 - p[2+j*3] ) / 4.0;
    phi[16+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 - p[0+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 4.0;
    phi[17+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 4.0;
    phi[18+j*20] = ( 1.0 + p[0+j*3] ) * ( 1.0 - p[0+j*3] ) 
                 * ( 1.0 + p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 4.0;
    phi[19+j*20] = ( 1.0 - p[0+j*3] ) * ( 1.0 + p[1+j*3] ) 
                 * ( 1.0 - p[1+j*3] ) * ( 1.0 + p[2+j*3] ) / 4.0;
  }
  return phi;
}
//****************************************************************************80

void basis_brick20_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK20_TEST verifies BASIS_BRICK20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
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
  int i;
  int j;
  int n;
  int node_num = 20;
  double *p;
  double *phi;
  double phi_sum;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "BASIS_BRICK20_TEST:\n";
  cout << "  Verify basis functions for element BRICK20.\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  n = node_num;

  p = nodes_brick20 ( );

  phi = basis_brick20 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    cout << "  ";
    for ( i = 0; i < node_num; i++ )
    {
      cout << setw(7) << phi[i+j*node_num];
    }
    cout << "\n";
  }

  delete [] p;
  delete [] phi;

  cout << "\n";
  cout << "  The basis function values at ANY point P\n";
  cout << "  should sum to 1:\n";
  cout << "\n";
  cout << "    ------------P-------------     PHI_SUM\n";
  cout << "\n";

  n = test_num;
  seed = 123456789;

  p = r8mat_uniform_01_new ( 3, n, &seed );

  phi = basis_brick20 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    phi_sum = r8vec_sum ( node_num, phi+j*node_num );
    cout << "  " << setw(8) << p[0+j*3]
         << "  " << setw(8) << p[1+j*3]
         << "  " << setw(8) << p[2+j*3]
         << "  " << setw(8) << phi_sum << "\n";
  }

  delete [] p;
  delete [] phi;

  return;
}
//****************************************************************************80

double *basis_brick27 ( int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK20: BRICK20 basis functions at natural coordinates.
//
//  Discussion:
//
//        8----19---7
//       /|         /
//     20 |   26   /|
//     /          / |
//    5----17----6  |
//    |   |      |  |
//    |  16---24-|-15
//    |  /|      | /|
//    |25 |  27  |23|        t
//    |/         |/ |        |   s
//   13----22---14  |        |  /
//    |   |      |  |        | /
//    |   |      |  |        |/
//    |   4--11--|--3        0---------r
//    |  /       | /        
//    | 12   21  |10       
//    |/         |/        
//    1----9-----2
//                   
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], natural coordinates of evaluation
//    points.
//
//    Output, double BASIS_BRICK27[27*N)], the basis function values.
//
{
  int j;
  double *phi;
  double rm;
  double rp;
  double rz;
  double sm;
  double sp;
  double sz;
  double tm;
  double tp;
  double tz;

  phi = new double[27*n];

  for ( j = 0; j < n; j++ )
  {
    rm = p[0+j*3] + 1.0;
    rz = p[0+j*3];
    rp = p[0+j*3] - 1.0;

    sm = p[1+j*3] + 1.0;
    sz = p[1+j*3];
    sp = p[1+j*3] - 1.0;

    tm = p[2+j*3] + 1.0;
    tz = p[2+j*3];
    tp = p[2+j*3] - 1.0;

    phi[0+j*27]  =        rz * rp      * sz * sp      * tz * tp / 8.0;
    phi[1+j*27]  =   rm * rz           * sz * sp      * tz * tp / 8.0;
    phi[2+j*27]  =   rm * rz      * sm * sz           * tz * tp / 8.0;
    phi[3+j*27]  =        rz * rp * sm * sz           * tz * tp / 8.0;
    phi[4+j*27]  =        rz * rp      * sz * sp * tm * tz      / 8.0;
    phi[5+j*27]  =   rm * rz           * sz * sp * tm * tz      / 8.0;
    phi[6+j*27]  =   rm * rz      * sm * sz      * tm * tz      / 8.0;
    phi[7+j*27]  =        rz * rp * sm * sz      * tm * tz      / 8.0;

    phi[8+j*27]  = - rm      * rp      * sz * sp      * tz * tp / 4.0;
    phi[9+j*27]  = - rm * rz      * sm      * sp      * tz * tp / 4.0;
    phi[10+j*27] = - rm      * rp * sm * sz           * tz * tp / 4.0;
    phi[11+j*27] = -      rz * rp * sm      * sp      * tz * tp / 4.0;
    phi[12+j*27] = -      rz * rp      * sz * sp * tm      * tp / 4.0;
    phi[13+j*27] = - rm * rz           * sz * sp * tm      * tp / 4.0;
    phi[14+j*27] = - rm * rz      * sm * sz      * tm      * tp / 4.0;
    phi[15+j*27] = -      rz * rp * sm * sz      * tm      * tp / 4.0;
    phi[16+j*27] = - rm      * rp      * sz * sp * tm * tz      / 4.0;
    phi[17+j*27] = - rm * rz      * sm      * sp * tm * tz      / 4.0;
    phi[18+j*27] = - rm      * rp * sm * sz      * tm * tz      / 4.0;
    phi[19+j*27] = -      rz * rp * sm      * sp * tm * tz      / 4.0;

    phi[20+j*27] =   rm      * rp * sm      * sp      * tz * tp / 2.0;
    phi[21+j*27] =   rm      * rp      * sz * sp * tm      * tp / 2.0;
    phi[22+j*27] =   rm * rz      * sm      * sp * tm      * tp / 2.0;
    phi[23+j*27] =   rm      * rp * sm * sz      * tm      * tp / 2.0;
    phi[24+j*27] =        rz * rp * sm      * sp * tm      * tp / 2.0;
    phi[25+j*27] =   rm      * rp * sm      * sp * tm * tz      / 2.0;

    phi[26+j*27] = - rm      * rp * sm      * sp * tm      * tp;
  }

  return phi;
}
//****************************************************************************80

void basis_brick27_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_BRICK27_TEST verifies BASIS_BRICK27.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
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
  int i;
  int j;
  int n;
  int node_num = 27;
  double *p;
  double *phi;
  double phi_sum;
  int seed;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "BASIS_BRICK27_TEST:\n";
  cout << "  Verify basis functions for element BRICK27.\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  n = node_num;

  p = nodes_brick27 ( );

  phi = basis_brick27 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    cout << "  ";
    for ( i = 0; i < node_num; i++ )
    {
      cout << setw(7) << phi[i+j*node_num];
    }
    cout << "\n";
  }

  delete [] p;
  delete [] phi;

  cout << "\n";
  cout << "  The basis function values at ANY point P\n";
  cout << "  should sum to 1:\n";
  cout << "\n";
  cout << "    ------------P-------------     PHI_SUM\n";
  cout << "\n";

  n = test_num;
  seed = 123456789;

  p = r8mat_uniform_01_new ( 3, n, &seed );

  phi = basis_brick27 ( n, p );

  for ( j = 0; j < n; j++ )
  {
    phi_sum = r8vec_sum ( node_num, phi+j*node_num );
    cout << "  " << setw(8) << p[0+j*3]
         << "  " << setw(8) << p[1+j*3]
         << "  " << setw(8) << p[2+j*3]
         << "  " << setw(8) << phi_sum << "\n";
  }

  delete [] p;
  delete [] phi;

  return;
}
//****************************************************************************80

void basis_mn_tet4 ( double t[3*4], int n, double p[], double phi[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_TET4: all bases at N points for a TET4 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a tetrahedron.
//
//    It works directly with these coordinates, and does not refer to a
//    reference element.
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
//  Reference:
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54.
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], the points where the basis functions
//    are to be evaluated.
//
//    Output, double PHI[4*N], the value of the basis functions
//    at the evaluation points.
//
{
  int j;
  double volume;
//
//           | x1 x2 x3 x4 |
//  Volume = | y1 y2 y3 y4 |
//           | z1 z2 z3 z4 |
//           |  1  1  1  1 |
//
  volume =                             
      t[0+0*3] * (                       
        t[1+1*3] * ( t[2+2*3] - t[2+3*3] )   
      - t[1+2*3] * ( t[2+1*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+1*3] - t[2+2*3] ) ) 
    - t[0+1*3] * (                       
        t[1+0*3] * ( t[2+2*3] - t[2+3*3] )   
      - t[1+2*3] * ( t[2+0*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+0*3] - t[2+2*3] ) ) 
    + t[0+2*3] * (                       
        t[1+0*3] * ( t[2+1*3] - t[2+3*3] )   
      - t[1+1*3] * ( t[2+0*3] - t[2+3*3] )   
      + t[1+3*3] * ( t[2+0*3] - t[2+1*3] ) ) 
    - t[0+3*3] * (                       
        t[1+0*3] * ( t[2+1*3] - t[2+2*3] )   
      - t[1+1*3] * ( t[2+0*3] - t[2+2*3] )   
      + t[1+2*3] * ( t[2+0*3] - t[2+1*3] ) );

  if ( volume == 0.0 )
  {
    cerr << "\n";
    cerr << "BASIS_MN_TET4 - Fatal error!\n";
    cerr << "  Element has zero volume.\n";
    exit ( 1 );
  }
//
//             | xp x2 x3 x4 |
//  Phi(1,P) = | yp y2 y3 y4 | / volume
//             | zp z2 z3 z4 |
//             |  1  1  1  1 |
//
  for ( j = 0; j < n; j++ )
  {
    phi[0+j*4] = (                           
        p[0+j*3] * (                         
          t[1+1*3] * ( t[2+2*3] - t[2+3*3] )   
        - t[1+2*3] * ( t[2+1*3] - t[2+3*3] )   
        + t[1+3*3] * ( t[2+1*3] - t[2+2*3] ) ) 
      - t[0+1*3] * (                           
          p[1+j*3] * ( t[2+2*3] - t[2+3*3] )   
        - t[1+2*3] * ( p[2+j*3] - t[2+3*3] )   
        + t[1+3*3] * ( p[2+j*3] - t[2+2*3] ) ) 
      + t[0+2*3] * (                           
          p[1+j*3] * ( t[2+1*3] - t[2+3*3] )   
        - t[1+1*3] * ( p[2+j*3] - t[2+3*3] )   
        + t[1+3*3] * ( p[2+j*3] - t[2+1*3] ) ) 
      - t[0+3*3] * (                           
          p[1+j*3] * ( t[2+1*3] - t[2+2*3] )   
        - t[1+1*3] * ( p[2+j*3] - t[2+2*3] )   
        + t[1+2*3] * ( p[2+j*3] - t[2+1*3] ) ) ) / volume;
//
//             | x1 xp x3 x4 |
//  Phi(2,P) = | y1 yp y3 y4 | / volume
//             | z1 zp z3 z4 |
//             |  1  1  1  1 |
//
    phi[1+j*4] = (                             
        t[0+0*3] * (                             
          p[1+j*3] * ( t[2+2*3] - t[2+3*3] )     
        - t[1+2*3] * ( p[2+j*3] - t[2+3*3] )     
        + t[1+3*3] * ( p[2+j*3] - t[2+2*3] ) )   
      - p[0+j*3]   * (                         
          t[1+0*3] * ( t[2+2*3] - t[2+3*3] )     
        - t[1+2*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - t[2+2*3] ) )   
      + t[0+2*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - p[2+j*3] ) ) 
      - t[0+3*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+2*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+2*3] )     
        + t[1+2*3] * ( t[2+0*3] - p[2+j*3] ) ) ) / volume;
//
//             | x1 x2 xp x4 |
//  Phi(3,P) = | y1 y2 yp y4 | / volume
//             | z1 z2 zp z4 |
//             |  1  1  1  1 |
//
    phi[2+j*4] = (                              
        t[0+0*3] * (                             
          t[1+1*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+1*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+1*3] - p[2+j*3] ) ) 
      - t[0+1*3] * (                             
          t[1+0*3] * ( p[2+j*3] - t[2+3*3] )     
        - p[1+j*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - p[2+j*3] ) ) 
      + p[0+j*3] * (                           
          t[1+0*3] * ( t[2+1*3] - t[2+3*3] )     
        - t[1+1*3] * ( t[2+0*3] - t[2+3*3] )     
        + t[1+3*3] * ( t[2+0*3] - t[2+1*3] ) )   
      - t[0+3*3] * (                             
          t[1+0*3] * ( t[2+1*3] - p[2+j*3] )   
        - t[1+1*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+1*3] ) ) ) / volume;
//
//             | x1 x2 x3 xp |
//  Phi(4,P) = | y1 y2 y3 yp | / volume
//             | z1 z2 z3 zp |
//             |  1  1  1  1 |
//
    phi[3+j*4] = (                             
        t[0+0*3] * (                             
          t[1+1*3] * ( t[2+2*3] - p[2+j*3] )   
        - t[1+2*3] * ( t[2+1*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+1*3] - t[2+2*3] ) )   
      - t[0+1*3] * (                             
          t[1+0*3] * ( t[2+2*3] - p[2+j*3] )   
        - t[1+2*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+2*3] ) )   
      + t[0+2*3] * (                             
          t[1+0*3] * ( t[2+1*3] - p[2+j*3] )   
        - t[1+1*3] * ( t[2+0*3] - p[2+j*3] )   
        + p[1+j*3] * ( t[2+0*3] - t[2+1*3] ) )   
      - p[0+j*3] * (                           
          t[1+0*3] * ( t[2+1*3] - t[2+2*3] )     
        - t[1+1*3] * ( t[2+0*3] - t[2+2*3] )     
        + t[1+2*3] * ( t[2+0*3] - t[2+1*3] ) ) ) / volume;
  }
  return;
}
//****************************************************************************80

void basis_mn_tet4_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_T4_TEST verifies BASIS_MN_TET4.
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
//    None.
//
{
# define NODE_NUM 4

  double *c;
  double c_sum;
  int i;
  int j;
  double *p;
  double phi1[NODE_NUM*1];
  double phi1_sum;
  double phi4[NODE_NUM*NODE_NUM];
  int seed;
  double *t;
  int test;
  int test_num = 5;

  cout << "\n";
  cout << "BASIS_MN_TET4_TEST:\n";
  cout << "  Verify basis functions for element TET4.\n";
  cout << "\n";
  cout << "  Number of nodes = " << NODE_NUM << "\n";

  t = r8mat_uniform_01_new ( 3, 4, &seed );

  cout << "\n";
  cout << "  Tetrahedron Nodes:\n";
  cout << "\n";
  for ( j = 0; j < NODE_NUM; j++ )
  {
    cout << "  "
         << setw(10) << t[0+j*3] << "  "
         << setw(10) << t[1+j*3] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  basis_mn_tet4 ( t, NODE_NUM, t, phi4 );

  for ( j = 0; j < NODE_NUM; j++ )
  {
    for ( i = 0; i < NODE_NUM; i++ )
    {
      cout << "  " << fixed << setprecision(4) << setw(10) << phi4[i+j*NODE_NUM];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The basis function values at ANY point P\n";
  cout << "  should sum to 1:\n";
  cout << "\n";
  cout << "    ------------P-------------    "
       << "-----------------PHI----------------   PHI_SUM\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    c = r8vec_uniform_01_new ( 4, &seed );

    c_sum = r8vec_sum ( 4, c );
    for ( i = 0; i < 4; i++ )
    {
      c[i] = c[i] / c_sum;
    }
    p = r8mat_mv ( 3, 4, t, c );

    basis_mn_tet4 ( t, 1, p, phi1 );

    phi1_sum = r8vec_sum ( NODE_NUM, phi1 );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << p[2]
         << "  " << setw(8) << phi1[0]
         << "  " << setw(8) << phi1[1]
         << "  " << setw(8) << phi1[2]
         << "  " << setw(8) << phi1[3]
         << "  " << setw(8) << phi1_sum << "\n";

    delete [] c;
    delete [] p;
  }

  delete [] t;

  return;
# undef NODE_NUM
}
//****************************************************************************80

void basis_mn_tet10 ( double t[3*4], int n, double p[], double phi[] )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_TET10: all bases at N points for a TET10 element.
//
//  Discussion:
//
//    The routine is given the coordinates of the vertices of a tetrahedron.
//
//    It works directly with these coordinates, and does not refer to a
//    reference element.
//
//    P1 through P4 are vertices.
//
//    P1 <= P5  <= P2
//    P2 <= P6  <= P3
//    P1 <= P7  <= P3
//    P1 <= P8  <= P4
//    P2 <= P9  <= P4
//    P3 <= P10 <= P4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Olgierd Zienkiewicz,
//    The Finite Element Method,
//    Sixth Edition,
//    Butterworth-Heinemann, 2005,
//    ISBN: 0750663200,
//    LC: TA640.2.Z54.
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double P[3*N], the points where the basis functions
//    are to be evaluated.
//
//    Output, double PHI[10*N], the value of the basis functions
//    at the evaluation points.
//
{
  int j;
  double *phi_linear;

  phi_linear = new double[4*n];

  basis_mn_tet4 ( t, n, p, phi_linear );

  for ( j = 0; j < n; j++ )
  {
    phi[0+j*10] = ( 2.0 * phi_linear[0+j*4]  - 1.0 ) * phi_linear[0+j*4];
    phi[1+j*10] = ( 2.0 * phi_linear[1+j*4]  - 1.0 ) * phi_linear[1+j*4];
    phi[2+j*10] = ( 2.0 * phi_linear[2+j*4]  - 1.0 ) * phi_linear[2+j*4];
    phi[3+j*10] = ( 2.0 * phi_linear[3+j*4]  - 1.0 ) * phi_linear[3+j*4];
    phi[4+j*10] =   4.0 * phi_linear[0+j*4]          * phi_linear[1+j*4];
    phi[5+j*10] =   4.0 * phi_linear[1+j*4]          * phi_linear[2+j*4];
    phi[6+j*10] =   4.0 * phi_linear[0+j*4]          * phi_linear[2+j*4];
    phi[7+j*10] =   4.0 * phi_linear[0+j*4]          * phi_linear[3+j*4];
    phi[8+j*10] =   4.0 * phi_linear[1+j*4]          * phi_linear[3+j*4];
    phi[9+j*10] =   4.0 * phi_linear[2+j*4]          * phi_linear[3+j*4];
  }
  delete [] phi_linear;

  return;
}
//****************************************************************************80

void basis_mn_tet10_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BASIS_MN_TET10_TEST verifies BASIS_MN_TET10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2009
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
  double *c;
  double c_sum;
  int i;
  int j;
  double *p;
  double p10[3*10];
  double phi1[10*1];
  double phi1_sum;
  double phi10[10*10];
  int seed;
  double *t;
  int test;
  int test_num = 5;

  seed = 123456789;

  cout << "\n";
  cout << "BASIS_MN_TET10_TEST:\n";
  cout << "  Verify basis functions for element TET10.\n";
  cout << "\n";
  cout << "  Number of nodes = 10.\n";

  t = r8mat_uniform_01_new ( 3, 4, &seed );

  cout << "\n";
  cout << "  Tetrahedron Nodes:\n";
  cout << "\n";
  for ( j = 0; j < 4; j++ )
  {
    cout << "  " << setw(8) << j
         << "  " << setw(14) << t[0+j*3]
         << "  " << setw(14) << t[1+j*3]
         << "  " << setw(14) << t[2+j*3] << "\n";
  }
 
  cout << "\n";
  cout << "  The basis function values at basis nodes\n";
  cout << "  should form the identity matrix.\n";
  cout << "\n";

  for ( i = 0; i < 3; i++ )
  {
    p10[i+0*3] = t[i+0*3];
    p10[i+1*3] = t[i+1*3];
    p10[i+2*3] = t[i+2*3];
    p10[i+3*3] = t[i+3*3];
    p10[i+4*3] = 0.5 * ( t[i+0*3] + t[i+1*3] );
    p10[i+5*3] = 0.5 * ( t[i+1*3] + t[i+2*3] );
    p10[i+6*3] = 0.5 * ( t[i+0*3] + t[i+2*3] );
    p10[i+7*3] = 0.5 * ( t[i+0*3] + t[i+3*3] );
    p10[i+8*3] = 0.5 * ( t[i+1*3] + t[i+3*3] );
    p10[i+9*3] = 0.5 * ( t[i+2*3] + t[i+3*3] );
  }

  basis_mn_tet10 ( t, 10, p10, phi10 );

  for ( i = 0; i < 10; i++ )
  {
    for ( j = 0; j < 10; j++ )
    {
      cout << "  " << fixed << setprecision(4) << setw(7) << phi10[i+j*10];
    }
    cout << "\n";
  }

  cout << "\n";
  cout << "  The basis function values at ANY point P\n";
  cout << "  should sum to 1:\n";
  cout << "\n";
  cout << "    ------------P-------------    ";
  cout << "----------------------------------------------------";
  cout << "PHI-----------------------------------------   PHI_SUM\n";
  cout << "\n";

  for ( test = 1; test <= test_num; test++ )
  {
    c = r8vec_uniform_01_new ( 4, &seed );

    c_sum = r8vec_sum ( 4, c );
    for ( i = 0; i < 4; i++ )
    {
      c[i] = c[i] / c_sum;
    }
    p = r8mat_mv ( 3, 4, t, c );

    basis_mn_tet10 ( t, 1, p, phi1 );
    phi1_sum = r8vec_sum ( 10, phi1 );

    cout << "  " << setw(8) << p[0]
         << "  " << setw(8) << p[1]
         << "  " << setw(8) << p[2]
         << "  " << setw(8) << phi1[0]
         << "  " << setw(8) << phi1[1]
         << "  " << setw(8) << phi1[2]
         << "  " << setw(8) << phi1[3]
         << "  " << setw(8) << phi1[4]
         << "  " << setw(8) << phi1[5]
         << "  " << setw(8) << phi1[6]
         << "  " << setw(8) << phi1[7]
         << "  " << setw(8) << phi1[8]
         << "  " << setw(8) << phi1[9]
         << "  " << setw(8) << phi1_sum << "\n";

    delete [] c;
    delete [] p;
  }

  delete [] t;

  return;
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

double *nodes_brick8 ( )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_BRICK8 returns the natural coordinates of the BRICK8 element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double NODES_BRICK8[3*8], the coordinates.
//
{
  double *p;
  static double p_save[3*8] = {
    -1.0, -1.0, -1.0, 
    +1.0, -1.0, -1.0, 
    +1.0, +1.0, -1.0, 
    -1.0, +1.0, -1.0, 
    -1.0, -1.0, +1.0, 
    +1.0, -1.0, +1.0, 
    +1.0, +1.0, +1.0, 
    -1.0, +1.0, +1.0 };

  p = r8mat_copy_new ( 3, 8, p_save );

  return p;
}
//****************************************************************************80

double *nodes_brick20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_BRICK20 returns the natural coordinates of the BRICK20 element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double NODES_BRICK20[3*20], the coordinates.
//
{
  double *p;
  static double p_save[3*20] = {
    -1.0, -1.0, -1.0, 
    +1.0, -1.0, -1.0, 
    +1.0, +1.0, -1.0, 
    -1.0, +1.0, -1.0, 
    -1.0, -1.0, +1.0, 
    +1.0, -1.0, +1.0, 
    +1.0, +1.0, +1.0, 
    -1.0, +1.0, +1.0, 
     0.0, -1.0, -1.0, 
    +1.0,  0.0, -1.0, 
     0.0, +1.0, -1.0, 
    -1.0,  0.0, -1.0, 
    -1.0, -1.0,  0.0, 
    +1.0, -1.0,  0.0, 
    +1.0, +1.0,  0.0, 
    -1.0, +1.0,  0.0, 
     0.0, -1.0, +1.0, 
    +1.0,  0.0, +1.0, 
     0.0, +1.0, +1.0, 
    -1.0,  0.0, +1.0 };

  p = r8mat_copy_new ( 3, 20, p_save );

  return p;
}
//****************************************************************************80

double *nodes_brick27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    NODES_BRICK27 returns the natural coordinates of the BRICK27 element.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double NODES_BRICK27[3*27], the coordinates.
//
{
  double *p;
  static double p_save[3*27] = {
    -1.0, -1.0, -1.0, 
    +1.0, -1.0, -1.0, 
    +1.0, +1.0, -1.0, 
    -1.0, +1.0, -1.0, 
    -1.0, -1.0, +1.0, 
    +1.0, -1.0, +1.0, 
    +1.0, +1.0, +1.0, 
    -1.0, +1.0, +1.0, 
     0.0, -1.0, -1.0, 
    +1.0,  0.0, -1.0, 
     0.0, +1.0, -1.0, 
    -1.0,  0.0, -1.0, 
    -1.0, -1.0,  0.0, 
    +1.0, -1.0,  0.0, 
    +1.0, +1.0,  0.0, 
    -1.0, +1.0,  0.0, 
     0.0, -1.0, +1.0, 
    +1.0,  0.0, +1.0, 
     0.0, +1.0, +1.0, 
    -1.0,  0.0, +1.0, 
     0.0,  0.0, -1.0, 
     0.0, -1.0,  0.0, 
    +1.0,  0.0,  0.0, 
     0.0, +1.0,  0.0, 
    -1.0,  0.0,  0.0, 
     0.0,  0.0, +1.0, 
     0.0,  0.0,  0.0 };

  p = r8mat_copy_new ( 3, 27, p_save );

  return p;
}
//****************************************************************************80

double *physical_to_reference_tet4 ( double t[], int n, double phy[] )

//****************************************************************************80
//
//  Purpose:
//
//    PHYSICAL_TO_REFERENCE_TET4 maps physical points to reference points.
//
//  Discussion:
//
//    Given the vertices of an order 4 physical tetrahedron and a point 
//    (X,Y,Z) in the physical tetrahedron, the routine computes the value 
//    of the corresponding point (R,S,T) in the reference tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices of the
//    physical tetrahedron.  The vertices are assumed to be the images of
//    (1,0,0), (0,1,0), (0,0,1) and (0,0,0) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double PHY[3*N], the coordinates of physical points
//    to be transformed.
//
//    Output, double PHYSICAL_TO_REFERENCE[3*N], the coordinates of the 
//    corresponding points in the reference tetrahedron.
//
{
  double a[3*3];
  int i;
  int j;
  double *ref;

  for ( j = 0; j < 3; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      a[i+j*3] = t[i+j*3] - t[i+3*3];
    }
  }

  ref = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      ref[i+j*3] = phy[i+j*3] - t[i+3*3];
    }
  }

  r8ge_fss ( 3, a, n, ref );

  return ref;
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

void r8ge_fss ( int n, double a[], int nb, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8GE_FSS factors and solves multiple R8GE systems.
//
//  Discussion:
//
//    The R8GE storage format is used for a "general" M by N matrix.  
//    A physical storage space is made for each logical entry.  The two 
//    dimensional logical array is mapped to a vector, in which storage is 
//    by columns.
//
//    This routine does not save the LU factors of the matrix, and hence cannot
//    be used to efficiently solve multiple linear systems, or even to
//    factor A at one time, and solve a single linear system at a later time.
//
//    This routine uses partial pivoting, but no pivot vector is required.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N].
//    On input, A is the coefficient matrix of the linear system.
//    On output, A is in unit upper triangular form, and
//    represents the U factor of an LU factorization of the
//    original coefficient matrix.
//
//    Input, int NB, the number of right hand sides.
//
//    Input/output, double B[N*NB], on input, the right hand sides.
//    on output, the solutions of the linear systems.
//
{
  int i;
  int ipiv;
  int j;
  int jcol;
  double piv;
  double t;

  for ( jcol = 1; jcol <= n; jcol++ )
  {
//
//  Find the maximum element in column I.
//
    piv = r8_abs ( a[jcol-1+(jcol-1)*n] );
    ipiv = jcol;
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( piv < r8_abs ( a[i-1+(jcol-1)*n] ) )
      {
        piv = r8_abs ( a[i-1+(jcol-1)*n] );
        ipiv = i;
      }
    }

    if ( piv == 0.0 )
    {
      cout << "\n";
      cout << "R8GE_FSS - Fatal error!\n";
      cout << "  Zero pivot on step " << jcol << "\n";
      return;
    }
//
//  Switch rows JCOL and IPIV, and X.
//
    if ( jcol != ipiv )
    {
      for ( j = 1; j <= n; j++ )
      {
        t                 = a[jcol-1+(j-1)*n];
        a[jcol-1+(j-1)*n] = a[ipiv-1+(j-1)*n];
        a[ipiv-1+(j-1)*n] = t;
      }
      for ( j = 0; j < nb; j++ )
      {
        t             = b[jcol-1+j*n];
        b[jcol-1+j*n] = b[ipiv-1+j*n];
        b[ipiv-1+j*n] = t;
      }
    }
//
//  Scale the pivot row.
//
    t = a[jcol-1+(jcol-1)*n];
    a[jcol-1+(jcol-1)*n] = 1.0;
    for ( j = jcol+1; j <= n; j++ )
    {
      a[jcol-1+(j-1)*n] = a[jcol-1+(j-1)*n] / t;
    }
    for ( j = 0; j < nb; j++ )
    {
      b[jcol-1+j*n] = b[jcol-1+j*n] / t;
    }
//
//  Use the pivot row to eliminate lower entries in that column.
//
    for ( i = jcol+1; i <= n; i++ )
    {
      if ( a[i-1+(jcol-1)*n] != 0.0 )
      {
        t = - a[i-1+(jcol-1)*n];
        a[i-1+(jcol-1)*n] = 0.0;
        for ( j = jcol+1; j <= n; j++ )
        {
          a[i-1+(j-1)*n] = a[i-1+(j-1)*n] + t * a[jcol-1+(j-1)*n];
        }
        for ( j = 0; j < nb; j++ )
        {
          b[i-1+j*n] = b[i-1+j*n] + t * b[jcol-1+j*n];
        }
      }
    }
  }
//
//  Back solve.
//
  for ( jcol = n; 2 <= jcol; jcol-- )
  {
    for ( i = 1; i < jcol; i++ )
    {
      for ( j = 0; j < nb; j++ )
      {
        b[i-1+j*n] = b[i-1+j*n] - a[i-1+(jcol-1)*n] * b[jcol-1+j*n];
      }
    }
  }

  return;
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

double r8mat_det_4d ( double a[4*4] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    The two dimensional array is stored as a one dimensional vector,
//    by COLUMNS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80

double *r8mat_mv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MV multiplies a matrix times a vector.
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
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV[M], the product A*X.
//
{
  int i;
  int j;
  double *y;

  y = new double[m];

  for ( i = 0; i < m; i++ )
  {
    y[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}
//****************************************************************************80

int r8mat_solve ( int n, int rhs_num, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
//
//  Discussion: 							    
//
//    A R8MAT is a doubly dimensioned array of double precision values, which
//    may be stored as a vector in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*N]
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
//    Input, int N, the order of the matrix.
//
//    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
//    must be at least 0.
//
//    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
//    to N the coefficient matrix, and in columns N+1 through
//    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
//    area has been destroyed, while the right hand sides have
//    been overwritten with the corresponding solutions.
//
//    Output, int R8MAT_SOLVE, singularity flag.
//    0, the matrix was not singular, the solutions were computed;
//    J, factorization failed on step J, and the solutions could not
//    be computed.
//
{
  double apivot;
  double factor;
  int i;
  int ipivot;
  int j;
  int k;
  double temp;

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a pivot row.
//
    ipivot = j;
    apivot = a[j+j*n];

    for ( i = j; i < n; i++ )
    {
      if ( fabs ( apivot ) < fabs ( a[i+j*n] ) )
      {
        apivot = a[i+j*n];
        ipivot = i;
      }
    }

    if ( apivot == 0.0 )
    {
      return j;
    }
//
//  Interchange.
//
    for ( i = 0; i < n + rhs_num; i++ )
    {
      temp          = a[ipivot+i*n];
      a[ipivot+i*n] = a[j+i*n];
      a[j+i*n]      = temp;
    }
//
//  A(J,J) becomes 1.
//
    a[j+j*n] = 1.0;
    for ( k = j; k < n + rhs_num; k++ )
    {
      a[j+k*n] = a[j+k*n] / apivot;
    }
//
//  A(I,J) becomes 0.
//
    for ( i = 0; i < n; i++ )
    {
      if ( i != j )
      {
        factor = a[i+j*n];
        a[i+j*n] = 0.0;
        for ( k = j; k < n + rhs_num; k++ )
        {
          a[i+k*n] = a[i+k*n] - factor * a[j+k*n];
        }
      }
    }
  }

  return 0;
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a new unit pseudorandom R8MAT.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has 
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
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

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *reference_tet4_sample ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TET4_SAMPLE: sample points in the reference tetrahedron.
//
//  Discussion:
//
//    This sampling method is not uniform.  The algorithm is simple.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double REFERENCE_TO_TET4_SAMPLE[3*N], points in the 
//    reference tetrahedron.
//
{
  double *c;
  double c_sum;
  int i;
  int j;
  double *ref;

  ref = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    c = r8vec_uniform_01_new ( 4, seed );
    c_sum = r8vec_sum ( 4, c );
    for ( i = 0; i < 3; i++ )
    {
      ref[i+j*3] = c[i] / c_sum;
    }
    delete [] c;
  }

  return ref;
}
//****************************************************************************80

double *reference_tet4_uniform ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TET4_UNIFORM: uniform sample points in the reference tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity 
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double REFERENCE_TET4_UNIFORM[3*N];
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    e = r8vec_uniform_01_new ( 4, seed );
    for ( i = 0; i < 4; i++ )
    {
      e[i] = - log ( e[i] );
    }
    e_sum = r8vec_sum ( 4, e );
    for (i = 0; i < 3; i++ )
    {
      x[i+j*3] = e[i] / e_sum;
    }

    delete [] e;
  }

  return x;
}
//****************************************************************************80

double *reference_tet4_uniform2 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TET4_UNIFORM2: uniform sample points in the reference tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Claudio Rocchini, Paolo Cignoni,
//    Generating Random Points in a Tetrahedron,
//    Journal of Graphics Tools,
//    Volume 5, Number 5, 2000, pages 9-12.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, double REFERENCE_TET4_UNIFORM2[3*N], the points.
//
{
  double *c;
  int i;
  int j;
  double t;
  double *x;

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    c = r8vec_uniform_01_new ( 3, seed );

    if ( 1.0 < c[0] + c[1] )
    {
      c[0] = 1.0 - c[0];
      c[1] = 1.0 - c[1];
    }

    if ( 1.0 < c[1] + c[2] )
    {
      t = c[2];
      c[2] = 1.0 - c[0] - c[1];
      c[1] = 1.0 - t;
    }
    else if ( 1.0 < c[0] + c[1] + c[2] )
    {
       t = c[2];
       c[2] = c[0] + c[1] + c[2] - 1.0;
       c[0] = 1.0 - c[1] - t;
    }
    for ( i = 0; i < 3; i++ )
    {
      x[i+j*3] = c[i];
    }
    delete [] c;
  }

  return x;
}
//****************************************************************************80

double *reference_to_physical_tet4 ( double t[], int n, double ref[] )

//****************************************************************************80
//
//  Purpose:
//
//    REFERENCE_TO_PHYSICAL_TET4 maps TET4 reference points to physical points.
//
//  Discussion:
//
//    Given the vertices of an order 4 physical tetrahedron and a point 
//    (R,S,T) in the reference tetrahedron, the routine computes the value 
//    of the corresponding image point (X,Y,Z) in physical space.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[3*4], the coordinates of the vertices.  
//    The vertices are assumed to be the images of (1,0,0), (0,1,0),
//    (0,0,1) and (0,0,0) respectively.
//
//    Input, int N, the number of points to transform.
//
//    Input, double REF[3*N], points in the reference element.
//
//    Output, double REFERENCE_TO_PHYSICAL_TET4[3*N], corresponding points in the
//    physical element.
//
{
  int i;
  int j;
  double *phy;

  phy = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      phy[i+j*3] =
          t[i+0*3] *         ref[0+j*3]
        + t[i+1*3] *                      ref[1+j*3]
        + t[i+2*3] *                                   ref[2+j*3]
        + t[i+3*3] * ( 1.0 - ref[0+j*3] - ref[1+j*3] - ref[2+j*3] );
    }
  }

  return phy;
}
//****************************************************************************80

double *tetrahedron_barycentric ( double tetra[3*4], double p[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_BARYCENTRIC returns the barycentric coordinates of a point.
//
//  Discussion:
//
//    The barycentric coordinates of a point P with respect to
//    a tetrahedron are a set of four values C(1:4), each associated
//    with a vertex of the tetrahedron.  The values must sum to 1.
//    If all the values are between 0 and 1, the point is contained
//    within the tetrahedron.
//
//    The barycentric coordinate of point X related to vertex A can be
//    interpreted as the ratio of the volume of the tetrahedron with 
//    vertex A replaced by vertex X to the volume of the original 
//    tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the vertices of the tetrahedron.
//
//    Input, double P[3], the point to be checked.
//
//    Output, double C[4], the barycentric coordinates of the point with
//    respect to the tetrahedron.
//
{
# define N 3
# define RHS_NUM 1

  double a[N*(N+RHS_NUM)];
  double *c;
  int info;
//
//  Set up the linear system
//
//    ( X2-X1  X3-X1  X4-X1 ) C1    X - X1
//    ( Y2-Y1  Y3-Y1  Y4-Y1 ) C2  = Y - Y1
//    ( Z2-Z1  Z3-Z1  Z4-Z1 ) C3    Z - Z1
//
//  which is satisfied by the barycentric coordinates.
//

  a[0+0*N] = tetra[0+1*3] - tetra[0+0*3];
  a[1+0*N] = tetra[1+1*3] - tetra[1+0*3];
  a[2+0*N] = tetra[2+1*3] - tetra[2+0*3];

  a[0+1*N] = tetra[0+2*3] - tetra[0+0*3];
  a[1+1*N] = tetra[1+2*3] - tetra[1+0*3];
  a[2+1*N] = tetra[2+2*3] - tetra[2+0*3];

  a[0+2*N] = tetra[0+3*3] - tetra[0+0*3];
  a[1+2*N] = tetra[1+3*3] - tetra[1+0*3];
  a[2+2*N] = tetra[2+3*3] - tetra[2+0*3];

  a[0+3*N] = p[0]         - tetra[0+0*3];
  a[1+3*N] = p[1]         - tetra[1+0*3];
  a[2+3*N] = p[2]         - tetra[2+0*3];
//
//  Solve the linear system.
//
  info = r8mat_solve ( N, RHS_NUM, a );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TETRAHEDRON_BARYCENTRIC - Fatal error!\n";
    cout << "  The linear system is singular.\n";
    cout << "  The input data does not form a proper tetrahedron.\n";
    exit ( 1 );
  }

  c = new double[4];

  c[1] = a[0+3*N];
  c[2] = a[1+3*N];
  c[3] = a[2+3*N];

  c[0] = 1.0 - c[1] - c[2] - c[3];

  return c;
# undef N
# undef RHS_NUM
}
//****************************************************************************80

double tetrahedron_volume ( double tetra[3*4] )

//****************************************************************************80
//
//  Purpose:
//
//    TETRAHEDRON_VOLUME computes the volume of a tetrahedron in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double TETRA[3*4], the coordinates of the vertices.
//
//    Output, double TETRAHEDRON_VOLUME, the volume of the tetrahedron.
//
{
  double a[4*4];
  int i;
  int j;
  double volume;

  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < 4; j++ )
    { 
      a[i+j*4] = tetra[i+j*3];
    }
  }

  i = 3;
  for ( j = 0; j < 4; j++ )
  {
    a[i+j*4] = 1.0;
  }

  volume = fabs ( r8mat_det_4d ( a ) ) / 6.0;

  return volume;
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
