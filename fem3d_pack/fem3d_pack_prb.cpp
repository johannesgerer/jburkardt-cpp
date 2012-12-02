# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem3d_pack.hpp"

int main ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM3D_PACK_PRB calls the various FEM3D_PACK tests.
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
{
  int i;
  timestamp ( );

  cout << "\n";
  cout << "FEM3D_PACK_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the FEM3D_PACK library.\n";

  basis_mn_tet4_test ( );
  basis_mn_tet10_test ( );
  basis_brick8_test ( );
  basis_brick20_test ( );
  basis_brick27_test ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM3D_PACK_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PHYSICAL_TO_REFERENCE_TET4 and REFERENCE_TO_PHYSICAL_TET4.
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
{
  int i;
  int j;
  int n = 10;
  double *phy;
  double *ref;
  double *ref2;
  int seed;
  double t[3*4] = {
    1.0,  2.0,  3.0,
    4.0,  1.0,  2.0,
    2.0,  4.0,  4.0,
    3.0,  2.0,  5.0 };

  seed = 123456789;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  For an order 4 tetrahedron,\n";
  cout << "  PHYSICAL_TO_REFERENCE_TET4 maps a physical point to\n";
  cout << "    a reference point.\n";
  cout << "  REFERENCE_TO_PHYSICAL_TET4 maps a reference point to\n";
  cout << "    a physical point.\n";
  cout << "\n";
  cout << "    ( R       S       T ) ==> ( X       Y       Z ) ";
  cout << "==> ( R2      S2      T2 )\n";
  cout << "\n";

  ref = reference_tet4_uniform ( n, &seed );

  phy = reference_to_physical_tet4 ( t, n, ref );
  ref2 = physical_to_reference_tet4 ( t, n, phy );

  for ( j = 0; j < n; j++ )
  {
    cout << "  " << setw(8) << ref[0+j*3]
         << "  " << setw(8) << ref[1+j*3]
         << "  " << setw(8) << ref[2+j*3]
         << "  " << setw(8) << phy[0+j*3]
         << "  " << setw(8) << phy[1+j*3]
         << "  " << setw(8) << phy[2+j*3]
         << "  " << setw(8) << ref2[0+j*3]
         << "  " << setw(8) << ref2[1+j*3]
         << "  " << setw(8) << ref2[2+j*3] << "\n";
  }

  delete [] phy;
  delete [] ref;
  delete [] ref2;

  return;
}
