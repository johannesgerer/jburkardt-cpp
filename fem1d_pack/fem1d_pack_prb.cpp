# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem1d_pack.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_PACK_PRB.
//
//  Discussion:
//
//    FEM1D_PACK_PRB tests the FEM1D_PACK library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "FEM1D_PACK_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FEM1D_PACK library.\n";

  test01 ( );
//test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM1D_PACK_PRB\n";
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
//    TEST01 verifies LOCAL_BASIS_1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2011
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

  int i;
  int j;
  int node_num = NODE_NUM;
  double node_x[NODE_NUM] = { 1.0, 2.0, 4.0, 4.5 };
  double *phi;
  double phi_matrix[NODE_NUM*NODE_NUM];
  double s;
  int seed;
  double x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  LOCAL_BASIS_1D evaluates the local basis functions\n";
  cout << "  for a 1D element.\n";
  cout << "\n";
  cout << "  Test that the basis functions, evaluated at the nodes,\n";
  cout << "  form the identity matrix.\n";
  cout << "\n";
  cout << "  Number of nodes = " << node_num << "\n";

  cout << "\n";
  cout << "  Node coordinates:\n";
  cout << "\n";
  for ( j = 0; j < node_num; j++ )
  {
    cout << "  " << setw(8) << j
         << "  " << setw(7) << node_x[j] << "\n";
  }

  for ( j = 0; j < node_num; j++ )
  {
    x = node_x[j];
    phi = local_basis_1d ( node_num, node_x, x );
    for ( i = 0; i < node_num; i++ )
    {
      phi_matrix[i+j*node_num] = phi[i];
    }
    delete [] phi;
  }

  r8mat_print ( node_num, node_num, phi_matrix, "  A(I,J) = PHI(I) at node (J):" );

  seed = 123456789;

  cout << "\n";
  cout << "  The PHI functions should sum to 1 at random X values:\n";
  cout << "\n";
  cout << "       X        Sum ( PHI(:)(X) )\n";
  cout << "\n";

  for ( j = 1; j <= 5; j++ )
  {
    x = r8_uniform ( 1.0, 4.5, &seed );
    phi = local_basis_1d ( node_num, node_x, x );
    s = r8vec_sum ( node_num, phi );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << s << "\n";
    delete [] phi;
  }

  return;
# undef NODE_NUM
}
