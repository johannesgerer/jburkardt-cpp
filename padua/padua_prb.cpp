# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "padua.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PADUA_PRB.
//
//  Discussion:
//
//    PADUA_PRB tests the PADUA library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "PADUA_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the PADUA library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PADUA_PRB\n";
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
//    TEST01 tests PADUA_ORDER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  int l;
  int n;

  cout << " \n";
  cout << "TEST01\n";
  cout << "  PADUA_ORDER converts the level L into the order N\n";
  cout << "  of any Padua rule.\n";
  cout << " \n";
  cout << "     L         N\n";
  cout << " \n";

  for ( l = 0; l <= 10; l++ )
  {
    n = padua_order ( l );
    cout << "  " << setw(4) << l
         << "  " << setw(8) << n << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PADUA_POINTS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  int l;
  string label;
  int n;
  double *xy;

  cout << " \n";
  cout << "TEST02\n";
  cout << "  PADUA_POINTS returns the points of a Padua rule.\n";

  for ( l = 0; l <= 10; l++ )
  {
    n = padua_order ( l );
    xy = padua_points ( l );
    label = " Level " + i4_to_string ( l ) + " Padua points:";
    r8mat_transpose_print ( 2, n, xy, label );
    delete [] xy;
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests PADUA_PLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  string filename;
  int l;
 
  cout << " \n";
  cout << "TEST03\n";
  cout << "  PADUA_PLOT plots the Padua points.\n";

  filename = "padua_00";

  for ( l = 0; l <= 10; l++ )
  {;
    padua_plot ( l, filename );
    filename_inc ( &filename );
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests PADUA_POINTS and PADUA_POINTS_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  int l;
  int n;
  double *xy1;
  double *x2;
  double *y2;

  cout << " \n";
  cout << "TEST04\n";
  cout << "  PADUA_POINTS computes the points of a Padua rule.\n";
  cout << "  PADUA_POINTS_SET looks them up in a table.\n";

  for ( l = 3; l <= 4; l++ )
  {
    n = padua_order ( l );
    xy1 = padua_points ( l );
    x2 = new double[n];
    y2 = new double[n];
    padua_points_set ( l, x2, y2 );
    cout << "\n";
    cout << "  Level " << l << " Padua points.\n";
    cout << "\n";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(4) << j << "  "
           << setw(14) << xy1[0+j*2] << "  "
           << setw(14) << xy1[1+j*2] << "\n";
      cout << "    " << "  "
           << setw(14) << x2[j] << "  "
           << setw(14) << y2[j] << "\n";
    }
    delete [] xy1;
    delete [] x2;
    delete [] y2;
  }

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests PADUA_WEIGHTS and PADUA_WEIGHTS_SET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2014
//
//  Author:
//
//    John Burkardt
//
{
  double diff;
  int j;
  int l;
  int n;
  double *w1;
  double *w2;

  cout << " \n";
  cout << "TEST05\n";
  cout << "  PADUA_WEIGHTS computes quadrature points of a Padua rule.\n";
  cout << "  PADUA_WEIGHTS_SET looks them up in a table.\n";

  for ( l = 3; l <= 4; l++ )
  {
    n = padua_order ( l );
    w1 = padua_weights ( l );
    w2 = padua_weights_set ( l );
    cout << "\n";
    cout << "  Level " << l << " Padua quadrature weights.\n";
    cout << "\n";
    diff = 0.0;
    for ( j = 0; j < n; j++ )
    {
      cout << setw(4) << j << "  "
           << setw(14) << w1[j] << "  "
           << setw(14) << w2[j] << "\n";
      diff = r8_max ( diff, fabs ( w1[j] - w2[j] ) );
    }
    cout << "\n";
    cout << "  Maximum difference = " << diff << "\n";
    delete [] w1;
    delete [] w2;
  }

  return;
}
