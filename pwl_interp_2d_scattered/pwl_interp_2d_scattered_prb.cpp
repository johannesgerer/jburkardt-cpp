# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "pwl_interp_2d_scattered.hpp"
# include "test_interp_2d.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( int prob );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PWL_INTERP_2D_SCATTERED_PRB.
//
//  Discussion:
//
//    PWL_INTERP_2D_SCATTERED_PRB tests the PWL_INTERP_2D_SCATTERED library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prob;
  int prob_num;

  timestamp ( );
  cout << "\n";
  cout << "PWL_INTERP_2D_SCATTERED_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the PWL_INTERP_2D_SCATTERED library.\n";
  cout << "  The R8LIB library is needed.\n";
  cout << "  This test also needs the TEST_INTERP_2D library.\n";

  test01 ( );
  test02 ( );
//
//  Numerical tests.
//
  prob_num = f00_num ( );

  for ( prob = 1; prob <= prob_num; prob++ )
  {
    test03 ( prob );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "PWL_INTERP_2D_SCATTERED_PRB:\n";
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
//    TEST01 tests R8TRIS2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int element_neighbor[3*2*9];
  int element_num;
  int element_order = 3;
  int node_num = 9;
  double node_xy[2*9] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5, 
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5, 
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle[3*2*9];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  R8TRIS2 computes the Delaunay triangulation of\n";
  cout << "  a set of nodes in 2D.\n";
//
//  Set up the Delaunay triangulation.
//
  r8tris2 ( node_num, node_xy, element_num, triangle, element_neighbor );

  triangulation_order3_print ( node_num, element_num, node_xy, 
    triangle, element_neighbor );

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PWL_INTERP_2D_SCATTERED_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int element_neighbor[3*2*9];
  int element_num;
  int element_order = 3;
  int i;
  int j;
  int k;
  int ni = 25;
  int node_num = 9;
  double node_xy[2*9] = {
       0.0, 0.0, 
       0.0, 1.0, 
       0.2, 0.5, 
       0.3, 0.6, 
       0.4, 0.5, 
       0.6, 0.4, 
       0.6, 0.5,
       1.0, 0.0, 
       1.0, 1.0 };
  int triangle[3*2*9];
  double x;
  double xyi[2*25];
  double y;
  double zd[9];
  double ze;
  double *zi;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PWL_INTERP_2D_SCATTERED_VALUE evaluates a\n";
  cout << "  piecewise linear interpolant to scattered data.\n";
//
//  Set up the Delaunay triangulation.
//
  r8tris2 ( node_num, node_xy, element_num, triangle, element_neighbor );

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( 0 < element_neighbor[i+j*3] )
      {
        element_neighbor[i+j*3] = element_neighbor[i+j*3] - 1;
      }
    }
  }
  triangulation_order3_print ( node_num, element_num, node_xy, 
    triangle, element_neighbor );
//
//  Define the Z data.
//
  for ( i = 0; i < node_num; i++ )
  {
    x = node_xy[0+i*2];
    y = node_xy[1+i*2];
    zd[i] = x + 2.0 * y;
  }
//
//  Define the interpolation points.
//
  k = 0;
  for ( i = 0; i <= 4; i++ )
  {
    for ( j = 0; j <= 4; j++ )
    {
      xyi[0+k*2] = ( i - 1 ) / 4.0;
      xyi[1+k*2] = ( j - 1 ) / 4.0;
      k = k + 1;
    }
  }
//
//  Evaluate the interpolant.
//
  zi = pwl_interp_2d_scattered_value ( node_num, node_xy, zd, element_num, 
    triangle, element_neighbor, ni, xyi );

  cout << "\n";
  cout << "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n";
  cout << "\n";
  for ( k = 0; k < ni; k++ )
  {
    ze = xyi[0+k*2] + 2.0 * xyi[1+k*2];
    cout << "  " << setw(4) << k
         << "  " << setw(10) << xyi[0+k*2]
         << "  " << setw(10) << xyi[1+k*2]
         << "  " << setw(10) << zi[k]
         << "  " << setw(10) << ze << "\n";
  }

  delete [] zi;

  return;
}
//****************************************************************************80

void test03 ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PWL_INTERP_2D_SCATTERED_VALUE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int *element_neighbor;
  int element_num;
  int g;
  int i;
  int j;
  int k;
  int nd;
  int ni = 25;
  double rms;
  int *triangle;
  double x;
  double *xd;
  double xi[25];
  double *xyd;
  double xyi[2*25];
  double y;
  double *yd;
  double yi[25];
  double *zd;
  double *ze;
  double *zi;

  g = 2;
  nd = g00_size ( g );

  cout << "\n";
  cout << "TEST03\n";
  cout << "  PWL_INTERP_2D_SCATTERED_VALUE evaluates a\n";
  cout << "  piecewise linear interpolant to scattered data.\n";
  cout << "  Here, we use grid number " << g << "\n";
  cout << "  with " << nd << " scattered points in the unit square\n";
  cout << "  on problem " << prob << "\n";
//
//  Get the data points and evaluate the function there.
//
  xd = new double[nd];
  yd = new double[nd];

  g00_xy ( g, nd, xd, yd );

  zd = new double[nd];
  f00_f0 ( prob, nd, xd, yd, zd );

  xyd = new double[2*nd];

  for ( i = 0; i < nd; i++ )
  {
    xyd[0+i*2] = xd[i];
    xyd[1+i*2] = yd[i];
  }
//
//  Set up the Delaunay triangulation.
//
  element_neighbor = new int[3*2*nd];
  triangle = new int[3*2*nd];

  r8tris2 ( nd, xyd, element_num, triangle, element_neighbor );

  for ( j = 0; j < element_num; j++ )
  {
    for ( i = 0; i < 3; i++ )
    {
      if ( 0 < element_neighbor[i+j*3] )
      {
        element_neighbor[i+j*3] = element_neighbor[i+j*3] - 1;
      }
    }
  }
//
//  Define the interpolation points.
//
  k = 0;
  for ( i = 1; i <= 5; i++ )
  {
    for ( j = 1; j <= 5; j++ )
    {
      xyi[0+k*2] = ( 2 * i - 1 ) / 10.0;
      xyi[1+k*2] = ( 2 * j - 1 ) / 10.0;
      k = k + 1;
    }
  }

  for ( k = 0; k < ni; k++ )
  {
    xi[k] = xyi[0+k*2];
    yi[k] = xyi[1+k*2];
  }
  ze = new double[ni];
  f00_f0 ( prob, ni, xi, yi, ze );
//
//  Evaluate the interpolant.
//
  zi = pwl_interp_2d_scattered_value ( nd, xyd, zd, element_num, 
    triangle, element_neighbor, ni, xyi );

  rms = 0.0;
  for ( k = 0; k < ni; k++ )
  {
    rms = rms + pow ( zi[k] - ze[k], 2 );
  }
  rms = sqrt ( rms / ( double ) ( ni ) );

  cout << "\n";
  cout << "  RMS error is " << rms << "\n";

  cout << "\n";
  cout << "     K      Xi(K)       Yi(K)       Zi(K)       Z(X,Y)\n";
  cout << "\n";

  for ( k = 0; k < ni; k++ )
  {
    cout << "  " << setw(4) << k
         << "  " << setw(10) << xyi[0+k*2]
         << "  " << setw(10) << xyi[1+k*2]
         << "  " << setw(10) << zi[k]
         << "  " << setw(10) << ze[k] << "\n";
  }

  delete [] element_neighbor;
  delete [] triangle;
  delete [] xd;
  delete [] xyd;
  delete [] yd;
  delete [] zd;
  delete [] ze;

  return;
}
