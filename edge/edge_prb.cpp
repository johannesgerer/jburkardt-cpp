# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "edge.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test035 ( );
void test036 ( );
void test037 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    EDGE_PRB tests the EDGE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "EDGE_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the EDGE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test035 ( );
  test036 ( );
  test037 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "EDGE_PRB\n";
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
//    TEST01 plots functions with jump discontinuities.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double *f;
  string header;
  int i;
  int n;
  int seed;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Plot 1D test functions.\n";

  test_num = 7;

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx1";
      f = fx1_vec ( n, x );
      title = "1D Test Function #1";
    }
    else if ( test == 2 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx2";
      f = fx2_vec ( n, x );
      title = "1D Test Function #2";
    }
    else if ( test == 3 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx3";
      f = fx3_vec ( n, x );
      title = "1D Test Function #3";
    }
    else if ( test == 4 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx4";
      f = fx4_vec ( n, x );
      title = "1D Test Function #4";
    }
    else if ( test == 5 )
    {
      n = 101;
      x_min = -1.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx5";
      f = fx5_vec ( n, x );
      title = "1D Test Function #5";
    }
    else if ( test == 6 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx6";
      f = fx6_vec ( n, x );
      title = "1D Test Function #6";
    }
    else if ( test == 7 )
    {
      n = 101;
      x_min = 0.0;
      x_max = +1.0;
      x = r8vec_linspace_new ( n, x_min, x_max );
      header = "fx7";
      f = fx7_vec ( n, x );
      title = "1D Test Function #7";
    }

    data_filename = header + "_data.txt";
    data_unit.open ( data_filename.c_str ( ) );
    for ( i = 0; i < n; i++ )
    {
      data_unit << x[i] 
                << "  " << f[i] << "\n";
    }
    data_unit.close ( );
    cout << "  Created data file '" << data_filename << "'\n";

    command_filename = header + "_commands.txt";
    command_unit.open ( command_filename.c_str ( ) );
    command_unit << "# " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "# Usage:\n";
    command_unit << "#  gnuplot < " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "set term png\n";
    command_unit << "set output '" << header << ".png'\n";
    command_unit << "set xlabel '<--- X --->'\n";
    command_unit << "set ylabel '<--- Y --->'\n";
    command_unit << "set title '" << title << "'\n";
    command_unit << "set grid\n";
    command_unit << "set style data lines\n";
    command_unit << "plot '" << data_filename 
                 << "' using 1:2 with points lt 3 pt 4 linecolor rgb 'blue'\n";
    command_unit << "quit\n";
    command_unit.close ( );
    cout << "  Created command file '" << command_filename << "'\n";

    delete [] f;
  }

  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 plots a function with a jump discontinuity along a circle.
//
//  Discussion:
//
//    This is example 4.1 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxy;
  string header;
  int i;
  int j;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Plot 2D test function #1 with jump along circle.\n";

  header = "fxy1";
  title = "2D test function #1 with discontinuity along circle";

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy1 ( x[i], y[j] );
      data_unit << x[i]
                << "  " << y[j]
                << "  " << fxy << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file '" << data_filename << "'\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << header << ".png'\n";
  command_unit << "set view 120, 77\n";
  command_unit << "set hidden3d\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << data_filename << "' with lines\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 plots a function with a jump discontinuity along a circle.
//
//  Discussion:
//
//    This is example 4.2 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxy;
  string header;
  int i;
  int j;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Plot 2D test function #2, the Shepp Logan phantom.\n";

  header = "fxy2";
  title = "2D test function #2, the Shepp Logan phantom";

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy2 ( x[i], y[j] );
      data_unit << x[i]
                << "  " << y[j]
                << "  " << fxy << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file '" << data_filename << "'\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << header << ".png'\n";
  command_unit << "set view 30, 75\n";
  command_unit << "set hidden3d\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << data_filename << "' with lines\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 plots a function with a jump discontinuity along a circle.
//
//  Discussion:
//
//    This is example 3.2 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the location of discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxy;
  string header;
  int i;
  int j;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST035:\n";
  cout << "  Plot 2D test function #3, the modified 2D Harten function.\n";

  header = "fxy3";
  title = "2D test function #3, the modified 2D Harten function";

  n = 101;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.0;
  y_max = +1.0;
  y = r8vec_linspace_new ( n, y_min, y_max );

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy3 ( x[i], y[j] );
      data_unit << x[i]
                << "  " << y[j]
                << "  " << fxy << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file '" << data_filename << "'\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << header << ".png'\n";
  command_unit << "set view 30, 75\n";
  command_unit << "set hidden3d\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << data_filename << "' with lines\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 plots a function with a derivative discontinuity.
//
//  Discussion:
//
//    This is example 3.1 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the location of discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxy;
  string header;
  int i;
  int j;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST036:\n";
  cout << "  Plot 2D test function #4, the discontinuous medium wave, P(x,t).\n";

  header = "fxy4";
  title = "2D test function #4, the discontinuous medium wave, P(x,t)";

  n = 101;
  x_min = -1.0;
  x_max = 0.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = 0.0;
  y_max = 0.1;
  y = r8vec_linspace_new ( n, y_min, y_max );

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy4 ( x[i], y[j] );
      data_unit << x[i]
                << "  " << y[j]
                << "  " << fxy << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file '" << data_filename << "'\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << header << ".png'\n";
  command_unit << "set view 30, 45\n";
  command_unit << "set hidden3d\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << data_filename << "' with lines\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 plots a function with a derivative discontinuity.
//
//  Discussion:
//
//    This is example 3.1 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the location of discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxy;
  string header;
  int i;
  int j;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double *y;
  double y_max;
  double y_min;

  cout << "\n";
  cout << "TEST037:\n";
  cout << "  Plot 2D test function #5, the discontinuous medium wave, U(x,t).\n";

  header = "fxy5";
  title = "2D test function #5, the discontinuous medium wave, U(x,t)";

  n = 101;
  x_min = -1.0;
  x_max = 0.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = 0.0;
  y_max = 0.1;
  y = r8vec_linspace_new ( n, y_min, y_max );

  data_filename = header + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      fxy = fxy5 ( x[i], y[j] );
      data_unit << x[i]
                << "  " << y[j]
                << "  " << fxy << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file '" << data_filename << "'\n";

  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << header << ".png'\n";
  command_unit << "set view 30, 45\n";
  command_unit << "set hidden3d\n";
  command_unit << "set timestamp\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "splot '" << data_filename << "' with lines\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file '" << command_filename << "'\n";

  delete [] x;
  delete [] y;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 plots slices of a 3D function.
//
//  Discussion:
//
//    Although the slice plots look uninteresting, there is a lot of detail
//    hidden in the data in variations that are not obvious at first.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Larry Shepp,
//    Computerized tomography and nuclear magnetic resonance,
//    Journal of Computer Assisted Tomography,
//    Volume 4, Number 1, February 1980, pages 94-107.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  double fxyz;
  string header;
  int i;
  int j;
  int k;
  int n;
  int test;
  int test_num;
  string title;
  double *x;
  double x_max;
  double x_min;
  double x_val;
  double *y;
  double y_max;
  double y_min;
  double y_val;
  double *z;
  double z_max;
  double z_min;
  double z_val;

  cout << "\n";
  cout << "TEST04:\n";
  cout << "  Plot 3D test function #1, the Shepp Logan 3D phantom.\n";

  test_num = 3;

  n = 101;
  x_min = -1.5;
  x_max = +1.5;
  x = r8vec_linspace_new ( n, x_min, x_max );
  y_min = -1.5;
  y_max = +1.5;
  y = r8vec_linspace_new ( n, y_min, y_max );
  z_min = -1.5;
  z_max = +1.5;
  z = r8vec_linspace_new ( n, z_min, z_max );

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      x_val = 0.0;
      title = "Slice X = 0.0";
      header = "fxyz1_x";
    }
    else if ( test == 2 )
    {
      y_val = 0.0;
      title = "Slice Y = 0.0";
      header = "fxyz1_y";
    }
    else if ( test == 3 )
    {
      z_val = - 0.1;
      title = "Slice Z = - 0.1";
      header = "fxyz1_z";
    }

    data_filename = header + "_data.txt";
    data_unit.open ( data_filename.c_str ( ) );
    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        if ( test == 1 )
        {
          fxyz = fxyz1 ( x_val, y[j], z[i] );
          data_unit <<  y[j]
                    << "  " << z[i]
                    << "  " << fxyz << "\n";
        }
        else if ( test == 2 )
        {
          fxyz = fxyz1 ( x[j], y_val, z[i] );
          data_unit << x[j]
                    << "  " << z[i]
                    << "  " << fxyz << "\n";
        }
        else if ( test == 3 )
        {
          fxyz = fxyz1 ( x[j], y[i], z_val );
          data_unit << x[j]
                    << "  " << y[i] 
                    << "  " << fxyz << "\n";
        }
      }
      data_unit << "\n";
    }
    data_unit.close ( );
    cout << "  Created data file '" << data_filename << "'\n";

    command_filename = header + "_commands.txt";
    command_unit.open ( command_filename.c_str ( ) );
    command_unit << "# " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "# Usage:\n";
    command_unit << "#  gnuplot < " << command_filename << "\n";
    command_unit << "#\n";
    command_unit << "set term png\n";
    command_unit << "set output '" << header << ".png'\n";
    command_unit << "set view 20, 75\n";
    command_unit << "set hidden3d\n";
    command_unit << "set timestamp\n";
    if ( test == 1 )
    {
      command_unit << "set xlabel '<--- Y --->'\n";
      command_unit << "set ylabel '<--- Z --->'\n";
      command_unit << "set zlabel '<--- X --->'\n";
    }
    else if ( test == 2 )
    {
      command_unit << "set xlabel '<--- X --->'\n";
      command_unit << "set ylabel '<--- Z --->'\n";
      command_unit << "set zlabel '<--- Y --->'\n";
    }
    else if ( test == 3 )
    {
      command_unit << "set xlabel '<--- X --->'\n";
      command_unit << "set ylabel '<--- Y --->'\n";
      command_unit << "set zlabel '<--- Z --->'\n";
    }
    command_unit << "set title '" << title << "'\n";
    command_unit << "set grid\n";
    command_unit << "set style data lines\n";
    command_unit << "splot '" << data_filename << "' with lines\n";
    command_unit << "quit\n";
    command_unit.close ( );
    cout << "  Created command file '" << command_filename << "'\n";
  }

  delete [] x;
  delete [] y;
  delete [] z;

  return;
}
