# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "hypersphere_properties.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HYPERSPHERE_PROPERTIES_PRB.
//
//  Discussion:
//
//    HYPERSPHERE_PROPERTIES_PRB tests the HYPERSPHERE_PROPERTIES library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HYPERSPHERE_PROPERTIES_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the HYPERSPHERE_PROPERTIES library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HYPERSPHERE_PROPERTIES_PRB:\n";
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
//    TEST01 tests the coordinate conversion routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  double err;
  int m;
  int n;
  double *r;
  int seed;
  int test;
  double *theta;
  double *x;
  double *x2;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the coordinate conversion routines:\n";
  cout << "  CARTESIAN_TO_HYPERSPHERE: X       -> R,Theta\n";
  cout << "  HYPERSPHERE_TO_CARTESIAN: R,Theta -> X.\n";
  cout << "\n";
  cout << "  Pick a random X, and compute X2 by converting X\n";
  cout << "  to hypersphere and back.  Consider norm of difference.\n";
  cout << "\n";
  cout << "  M    || X - X2 ||\n";

  seed = 123456789;

  n = 1;
  r = new double[n];

  for ( m = 1; m <= 5; m++ )
  {
    cout << "\n";

    theta = new double[(m-1)*n];

    for ( test = 1; test <= 5; test++ )
    {
      x = r8mat_uniform_01_new ( m, n, seed );
      c = r8vec_uniform_01_new ( m, seed );
      cartesian_to_hypersphere ( m, n, c, x, r, theta );
      x2 = hypersphere_to_cartesian ( m, n, c, r, theta );
      err = r8mat_norm_fro_affine ( m, n, x, x2 );
      cout << "  " << setw(2) << m
           << "  " << setw(14) << err << "\n";
      delete [] c;
      delete [] x;
      delete [] x2;
    }
    delete [] theta;
  }

  delete [] r;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HYPERSPHERE_01_SURFACE_UNIFORM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int seed;
  int test;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HYPERSPHERE_01_SURFACE_UNIFORM samples uniformly from the\n";
  cout << "  surface of the unit hypersphere\n";

  seed = 123456789;

  n = 1;
  for ( m = 1; m <= 5; m++ )
  {
    for ( test = 1; test <= 3; test++ )
    {
      x = hypersphere_01_surface_uniform ( m, n, seed );
      r8vec_transpose_print ( m, x, "  Random hypersphere point:" );
      delete [] x;
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
//    TEST03 tests HYPERSPHERE_01_AREA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  double area2;
  int m;
  int n_data;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  HYPERSPHERE_01_AREA evaluates the area of the unit\n";
  cout << "  hypersphere in M dimensions.\n";
  cout << "\n";
  cout << "       M      Exact       Computed\n";
  cout << "              Area        Area\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypersphere_01_area_values ( n_data, m, area );

    if ( n_data == 0 )
    {
      break;
    }

    area2 = hypersphere_01_area ( m );

    cout << "  " << setw(6) << m
         << "  " << setw(10) << area
         << "  " << setw(10) << area2 << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests HYPERSPHERE_01_VOLUME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n_data;
  double volume;
  double volume2;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  HYPERSPHERE_01_VOLUME evaluates the area of the unit\n";
  cout << "  hypersphere in M dimensions.\n";
  cout << "  HYPERSPHERE_01_VOLUME_VALUES returns some test values.\n";
  cout << "\n";
  cout << "       M      Exact       Computed\n";
  cout << "              Volume      Volume\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    hypersphere_01_volume_values ( n_data, m, volume );

    if ( n_data == 0 )
    {
      break;
    }

    volume2 = hypersphere_01_volume ( m );

    cout << "  " << setw(6) << m
         << "  " << setw(10) << volume
         << "  " << setw(10) << volume2 << "\n";
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HYPERSPHERE_AREA, HYPERSPHERE_VOLUME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 May 2013
//
//  Author:
//
//    John Burkardt
//
{
  double area;
  int m;
  double r;
  double volume;

  r = 1.5;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  For a hypersphere in M dimensions:\n";
  cout << "  HYPERSPHERE_AREA computes the area\n";
  cout << "  HYPERSPHERE_VOLUME computes the volume.\n";
  cout << "\n";
  cout << "  Notice that both quantities eventually decrease\n";
  cout << "\n";
  cout << "  We use a radius of R = " << r << "\n";
  cout << "\n";
  cout << "    M        Area          Volume    Area / Volume \n";
  cout << "\n";

  for ( m = 1; m <= 20; m++ )
  {
    area = hypersphere_area ( m, r );
    volume = hypersphere_volume ( m, r );
    cout << "  " << setw(3) << m
         << "  " << setw(14) << area
         << "  " << setw(14) << volume
         << "  " << setw(14) << area / volume << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests the stereographic mapping.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 December 2013
//
//  Author:
//
//    John Burkardt
//
{
  double err;
  int m;
  int n;
  int seed;
  int test;
  double *x1;
  double *x2;
  double *x3;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Test the stereographic mapping:\n";
  cout << "  HYPERSPHERE_STEREOGRAPH maps hypersphere points to the plane.\n";
  cout << "  HYPERSPHERE_STEREOGRAPH_INVERSE inverts the mapping.\n";
  cout << "\n";
  cout << "  Pick a random X1 on the hypersphere.\n";
  cout << "  Map it to a point X2 on the plane.\n";
  cout << "  Map it back to a point X3 on the hypersphere.\n";
  cout << "  Consider norm of difference.\n";
  cout << "\n";
  cout << "  M    || X1 - X3 ||\n";

  seed = 123456789;

  n = 1;
  for ( m = 2; m <= 5; m++ )
  {
    cout << "\n";
    for ( test = 1; test <= 5; test++ )
    {
      x1 = hypersphere_01_surface_uniform ( m, n, seed );
      x2 = hypersphere_stereograph ( m, n, x1 );
      x3 = hypersphere_stereograph_inverse ( m, n, x2 );
      err = r8mat_norm_fro_affine ( m, n, x1, x3 );
      cout << "  " << setw(2) << m
           << "  " << setw(14) << err << "\n";
      delete [] x1;
      delete [] x2;
      delete [] x3;
    }
  }
  return;
}

