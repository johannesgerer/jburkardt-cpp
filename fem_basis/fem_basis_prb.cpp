# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "fem_basis.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FEM_BASIS_PRB
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "FEM_BASIS_PRB:\n";
  cout << "  C++ version.\n";
  cout << "  Test the FEM_BASIS library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FEM_BASIS_PRB:\n";
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
//    TEST01 tests FEM_BASIS_1D
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  double lij;
  double x1;
  double x2;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  FEM_BASIS_1D evaluates an arbitrary\n";
  cout << "  basis function over an interval.\n";

  i1 = 2;
  j1 = 1;
  d = i1 + j1;
  x1 = r8_fraction ( i1, d );
  cout << "\n";
  cout << "   I   J        X      L(I,J)(X)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1
       << "  " << setw(2) << j1
       << "  " << setw(10) << x1
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( i2 = 0; i2 <= d; i2++ )
  {
    j2 = d - i2;
    x2 = r8_fraction ( i2, d );
    lij = fem_basis_1d ( i1, j1, x2 );
    cout << "  " << setw(2) << i2
         << "  " << setw(2) << j2
         << "  " << setw(10) << x2
         << "  " << setw(14) << lij << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests FEM_BASIS_2D
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  int k1;
  int k2;
  double lijk;
  double x1;
  double x2;
  double y1;
  double y2;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  FEM_BASIS_2D evaluates an arbitrary triangular\n";
  cout << "  basis function.\n";

  i1 = 1;
  j1 = 0;
  k1 = 2;
  d = i1 + j1 + k1;
  x1 = r8_fraction ( i1, d );
  y1 = r8_fraction ( j1, d );
  cout << "\n";
  cout << "   I   J   K        X           Y      L(I,J,K)(X,Y)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1
       << "  " << setw(2) << j1 
       << "  " << setw(2) << k1
       << "  " << setw(10) << x1
       << "  " << setw(10) << y1
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( j2 = 0; j2 <= d; j2++ )
  {
    for ( i2 = 0; i2 <= d - j2; i2++ )
    {
      k2 = d - i2 - j2;
      x2 = r8_fraction ( i2, d );
      y2 = r8_fraction ( j2, d );
      lijk = fem_basis_2d ( i1, j1, k1, x2, y2 );
      cout << "  " << setw(2) << i2
           << "  " << setw(2) << j2 
           << "  " << setw(2) << k2
           << "  " << setw(10) << x2
           << "  " << setw(10) << y2
           << "  " << setw(14) << lijk << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests FEM_BASIS_3D
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i1;
  int i2;
  int j1;
  int j2;
  int k1;
  int k2;
  int l1;
  int l2;
  double lijkl;
  double x1;
  double x2;
  double y1;
  double y2;
  double z1;
  double z2;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  FEM_BASIS_3D evaluates an arbitrary tetrahedral\n";
  cout << "  basis function.\n";

  i1 = 1;
  j1 = 0;
  k1 = 2;
  l1 = 1;
  d = i1 + j1 + k1 + l1;
  x1 = r8_fraction ( i1, d );
  y1 = r8_fraction ( j1, d );
  z1 = r8_fraction ( k1, d );
  cout << "\n";
  cout << "   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1 
       << "  " << setw(2) << j1 
       << "  " << setw(2) << k1
       << "  " << setw(2) << l1
       << "  " << setw(10) << x1
       << "  " << setw(10) << y1
       << "  " << setw(10) << z1
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( k2 = 0; k2 <= d; k2++ )
  {
    for ( j2 = 0; j2 <= d - k2; j2++ )
    {
      for ( i2 = 0; i2 <= d - j2 - k2; i2++ )
      {
        l2 = d - i2 - j2 - k2;
        x2 = r8_fraction ( i2, d );
        y2 = r8_fraction ( j2, d );
        z2 = r8_fraction ( k2, d );
        lijkl = fem_basis_3d ( i1, j1, k1, l1, x2, y2, z2 );
        cout << "  " << setw(2) << i2 
             << "  " << setw(2) << j2 
             << "  " << setw(2) << k2
             << "  " << setw(2) << l2
             << "  " << setw(10) << x2
             << "  " << setw(10) << y2
             << "  " << setw(10) << z2
             << "  " << setw(14) << lijkl << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests FEM_BASIS_MD, repeating TEST01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int i1[2];
  int i2[2];
  double l;
  int m = 1;
  int p1;
  double x1[1];
  double x2[1];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  FEM_BASIS_MD evaluates an arbitrary\n";
  cout << "  basis function over an M-dimensional simplex.\n";

  i1[0] = 2;
  i1[1] = 1;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  cout << "\n";
  cout << "   I   J        X      L(I,J)(X)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1[0]
       << "  " << setw(2) << i1[1]
       << "  " << setw(10) << x1[0]
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( p1 = 0; p1 <= d; p1++ )
  {
    i2[0] = p1;
    i2[1] = d - i2[0];
    for ( i = 0; i < m; i++ )
    {
      x2[i] = r8_fraction ( i2[i], d );
    }
    l = fem_basis_md ( m, i1, x2 );
    cout << "  " << setw(2) << i2[0]
         << "  " << setw(2) << i2[1]
         << "  " << setw(10) << x2[0]
         << "  " << setw(14) << l << "\n";
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests FEM_BASIS_MD, repeating TEST02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int i1[3];
  int i2[3];
  double l;
  int m = 2;
  int p1;
  int p2;
  double x1[2];
  double x2[2];

  cout << "\n";
  cout << "TEST05\n";
  cout << "  FEM_BASIS_MD evaluates an arbitrary\n";
  cout << "  basis function over an M-dimensional simplex.\n";

  i1[0] = 1;
  i1[1] = 0;
  i1[2] = 2;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  cout << "\n";
  cout << "   I   J   K        X           Y      L(I,J,K)(X,Y)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1[0]
       << "  " << setw(2) << i1[1]
       << "  " << setw(2) << i1[2]
       << "  " << setw(10) << x1[0]
       << "  " << setw(10) << x1[1]
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( p2 = 0; p2 <= d; p2++ )
  {
    i2[1] = p2;
    for ( p1 = 0; p1 <= d - p2; p1++ )
    {
      i2[0] = p1;
      i2[2] = d - i2[0] - i2[1];
      for ( i = 0; i < m; i++ )
      {
        x2[i] = r8_fraction ( i2[i], d );
      }
      l = fem_basis_md ( m, i1, x2 );
      cout << "  " << setw(2) << i2[0]
           << "  " << setw(2) << i2[1]
           << "  " << setw(2) << i2[2]
           << "  " << setw(10) << x2[0]
           << "  " << setw(10) << x2[1]
           << "  " << setw(14) << l << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests FEM_BASIS_MD, repeating TEST03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int i;
  int i1[4];
  int i2[4];
  double l;
  int m = 3;
  int p1;
  int p2;
  int p3;
  double x1[3];
  double x2[3];

  cout << "\n";
  cout << "TEST06\n";
  cout << "  FEM_BASIS_MD evaluates an arbitrary\n";
  cout << "  basis function over an M-dimensional simplex.\n";

  i1[0] = 1;
  i1[1] = 0;
  i1[2] = 2;
  i1[3] = 1;
  d = i4vec_sum ( m + 1, i1 );
  for ( i = 0; i < m; i++ )
  {
    x1[i] = r8_fraction ( i1[i], d );
  }
  cout << "\n";
  cout << "   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1[0]
       << "  " << setw(2) << i1[1]
       << "  " << setw(2) << i1[2]
       << "  " << setw(2) << i1[3]
       << "  " << setw(10) << x1[0]
       << "  " << setw(10) << x1[1]
       << "  " << setw(10) << x1[2]
       << "  " << setw(14) << 1.0 << "\n";
  cout << "\n";
  for ( p3 = 0; p3 <= d; p3++ )
  {
    i2[2] = p3;
    for ( p2 = 0; p2 <= d - p3; p2++ )
    {
      i2[1] = p2;
      for ( p1 = 0; p1 <= d - p3 - p2; p1++ )
      {
        i2[0] = p1;
        i2[3] = d - i2[0] - i2[1] - i2[2];
        for ( i = 0; i < m; i++ )
        {
          x2[i] = r8_fraction ( i2[i], d );
        }
        l = fem_basis_md ( m, i1, x2 );
        cout << "  " << setw(2) << i2[0]
             << "  " << setw(2) << i2[1]
             << "  " << setw(2) << i2[2]
             << "  " << setw(2) << i2[3]
             << "  " << setw(10) << x2[0]
             << "  " << setw(10) << x2[1]
             << "  " << setw(10) << x2[2]
             << "  " << setw(14) << l << "\n";
      }
    }
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests FEM_BASIS_PRISM_TRIANGLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int di;
  int dj;
  int i1[3] = { 2, 0, 0 };
  int i2[3];
  int i_0;
  int i_1;
  int j1[2] = { 1, 1 };
  int j2[2];
  int j_0;
  double xyz1[3];
  double xyz2[3];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary\n";
  cout << "  basis function over a right triangular prism.\n";
  cout << "\n";
  cout << "  Here, we generate basis functions which can be\n";
  cout << "  up to degree 2 in X and Y, and up to degree 2 in Z.\n";
  cout << "\n";
  cout << "  Choose a node N1, define the basis function associated\n";
  cout << "  with that node, and then evaluate it at all other nodes.\n";

  di = i4vec_sum ( 3, i1 );
  xyz1[0] = r8_fraction ( i1[0], di );
  xyz1[1] = r8_fraction ( i1[1], di );

  dj = i4vec_sum ( 2, j1 );
  xyz1[2] = r8_fraction ( j1[0], dj );

  cout << "\n";
  cout << "  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1[0]
       << "  " << setw(2) << i1[1]
       << "  " << setw(2) << i1[2]
       << "  " << setw(2) << j1[0]
       << "  " << setw(2) << j1[1]
       << "  " << setw(10) << xyz1[0]
       << "  " << setw(10) << xyz1[1]
       << "  " << setw(10) << xyz1[2]
       << "  " << setw(10) << 1.0 << "\n";

  cout << "\n";

  for ( i_0 = 0; i_0 <= di; i_0++ )
  {
    i2[0] = i_0;
    xyz2[0] = r8_fraction ( i2[0], di );
    for ( i_1 = 0; i_1 <= di - i2[0]; i_1++ )
    {
      i2[1] = i_1;
      xyz2[1] = r8_fraction ( i2[1], di );
      i2[2] = di - i2[0] - i2[1];
      for( j_0 = 0; j_0 <= dj; j_0++ )
      {
        j2[0] = j_0;
        j2[1] = dj - j2[0];
        xyz2[2] = r8_fraction ( j2[0], dj );

        b = fem_basis_prism_triangle ( i1, j1, xyz2 );

        cout << "  " << setw(2) << i2[0]
             << "  " << setw(2) << i2[1]
             << "  " << setw(2) << i2[2]
             << "  " << setw(2) << j2[0]
             << "  " << setw(2) << j2[1]
             << "  " << setw(10) << xyz2[0]
             << "  " << setw(10) << xyz2[1]
             << "  " << setw(10) << xyz2[2]
             << "  " << setw(10) << b << "\n";
     }
    }
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests FEM_BASIS_PRISM_TRIANGLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2010
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int di;
  int dj;
  int i1[3] = { 2, 0, 1 };
  int i2[3];
  int i_0;
  int i_1;
  int j1[2] = { 1, 0 };
  int j2[2];
  int j_0;
  double xyz1[3];
  double xyz2[3];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary\n";
  cout << "  basis function over a right triangular prism.\n";
  cout << "\n";
  cout << "  Here, we generate basis functions which can be\n";
  cout << "  up to degree 3 in X and Y, and up to degree 1 in Z.\n";
  cout << "\n";
  cout << "  Choose a node N1, define the basis function associated\n";
  cout << "  with that node, and then evaluate it at all other nodes.\n";

  di = i4vec_sum ( 3, i1 );
  xyz1[0] = r8_fraction ( i1[0], di );
  xyz1[1] = r8_fraction ( i1[1], di );

  dj = i4vec_sum ( 2, j1 );
  xyz1[2] = r8_fraction ( j1[0], dj );

  cout << "\n";
  cout << "  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)\n";
  cout << "\n";
  cout << "  " << setw(2) << i1[0]
       << "  " << setw(2) << i1[1]
       << "  " << setw(2) << i1[2]
       << "  " << setw(2) << j1[0]
       << "  " << setw(2) << j1[1]
       << "  " << setw(10) << xyz1[0]
       << "  " << setw(10) << xyz1[1]
       << "  " << setw(10) << xyz1[2]
       << "  " << setw(10) << 1.0 << "\n";
  cout << "\n";

  for ( i_0 = 0; i_0 <= di; i_0++ )
  {
    i2[0] = i_0;
    xyz2[0] = r8_fraction ( i2[0], di );
    for ( i_1 = 0; i_1 <= di - i2[0]; i_1++ )
    {
      i2[1] = i_1;
      xyz2[1] = r8_fraction ( i2[1], di );
      i2[2] = di - i2[0] - i2[1];
      for( j_0 = 0; j_0 <= dj; j_0++ )
      {
        j2[0] = j_0;
        j2[1] = dj - j2[0];
        xyz2[2] = r8_fraction ( j2[0], dj );

        b = fem_basis_prism_triangle ( i1, j1, xyz2 );

        cout << "  " << setw(2) << i2[0]
             << "  " << setw(2) << i2[1]
             << "  " << setw(2) << i2[2]
             << "  " << setw(2) << j2[0]
             << "  " << setw(2) << j2[1]
             << "  " << setw(10) << xyz2[0]
             << "  " << setw(10) << xyz2[1]
             << "  " << setw(10) << xyz2[2]
             << "  " << setw(10) << b << "\n";
      }
    }
  }
  return;
}
