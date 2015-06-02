# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "stochastic_diffusion.hpp"

int main ( );
void bnt_contour ( );
void elman_contour ( );
void ntw_contour ( );
void xk_contour ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for STOCHASTIC_DIFFUSION_PRB.
//
//  Discussion:
//
//    STOCHASTIC_DIFFUSION_PRB tests the STOCHASTIC_DIFFUSION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "STOCHASTIC_DIFFUSION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the STOCHASTIC_DIFFUSION library.\n";

  bnt_contour ( );
  elman_contour ( );
  ntw_contour ( );
  xk_contour ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "STOCHASTIC_DIFFUSION_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void bnt_contour ( )

//****************************************************************************80
//
//  Purpose:
//
//    BNT_CONTOUR displays contour plots of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is compute by DIFFUSIVITY_2D_BNT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ivo Babuska, Fabio Nobile, Raul Tempone,
//    A stochastic collocation method for elliptic partial differential equations
//    with random input data,
//    SIAM Journal on Numerical Analysis,
//    Volume 45, Number 3, 2007, pages 1005-1034.
//
{
  string command_filename = "bnt_commands.txt";
  ofstream command_unit;
  string data_filename = "bnt_data.txt";
  ofstream data_unit;
  double *dc;
  double dc0;
  int i;
  int j;
  int m = 4;
  int n;
  int nx = 41;
  int ny = 31;
  double *omega;
  int seed;
  double *xmat;
  double *xvec;
  double *ymat;
  double *yvec;

  cout << "\n";
  cout << "BNT_CONTOUR\n";
  cout << "  Display contour or surface plots of the stochastic\n";
  cout << "  diffusivity function defined by DIFFUSIVITY_2D_BNT.\n";
  cout << "\n";
  cout << "  The first plot uses uniform random values for OMEGA.\n";
  cout << "  The second uses Gaussian (normal) random values.\n";
//
//  Set the spatial grid.
//
  xvec = r8vec_linspace_new ( nx, -1.5, 0.0 );
  yvec = r8vec_linspace_new ( ny, -0.4, 0.8 );

  xmat = new double[nx*ny];
  ymat = new double[nx*ny];
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
//
//  Sample OMEGA.
//
  seed = 123456789;
  omega = r8vec_uniform_01_new ( m, seed );
//
//  Compute the diffusivity field.
//
  dc0 = 10.0;
  n = nx * ny;
  dc = diffusivity_2d_bnt ( dc0, omega, n, xmat, ymat );
//
//  Create a data file.
//
  data_unit.open ( data_filename.c_str() );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      data_unit 
        << "  " << setw(14) << xmat[i+j*nx]
        << "  " << setw(14) << ymat[i+j*nx]
        << "  " << setw(14) << dc[i+j*nx] << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'.\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str() );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'bnt_contour.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set zlabel '<---DC(X,Y)--->'\n";
  command_unit << "set title 'BNT Stochastic diffusivity function'\n";
  command_unit << "set contour\n";
  command_unit << "set timestamp\n";
  command_unit << "set cntrparam levels 10\n";
  command_unit << "#set view map\n";
  command_unit << "set view 75, 75\n";
  command_unit << "unset key\n";
  command_unit << "splot '" << data_filename << "'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] dc;
  delete [] omega;
  delete [] xmat;
  delete [] xvec;
  delete [] ymat;
  delete [] yvec;

  return;
}
//****************************************************************************80

void elman_contour ( )

//****************************************************************************80
//
//  Purpose:
//
//    ELMAN_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is compute by DIFFUSIVITY_2D_ELMAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Howard Elman, Darran Furnaval,
//    Solving the stochastic steady-state diffusion problem using multigrid,
//    IMA Journal on Numerical Analysis,
//    Volume 27, Number 4, 2007, pages 675-688.
//
{
  double a;
  double cl;
  string command_filename = "elman_commands.txt";
  ofstream command_unit;
  string data_filename = "elman_data.txt";
  ofstream data_unit;
  double *dc;
  double dc0;
  int i;
  int j;
  int m;
  int m_1d = 5;
  int nx = 51;
  int ny = 51;
  double *omega;
  int seed;
  double *xmat;
  double *xvec;
  double *ymat;
  double *yvec;

  cout << "\n";
  cout << "ELMAN_CONTOUR\n";
  cout << "  Display contour or surface plots of the stochastic\n";
  cout << "  diffusivity function defined by DIFFUSIVITY_2D_ELMAN.\n";
//
//  Set the spatial grid.
//
  a = 1.0;
  xvec = r8vec_linspace_new ( nx, -a, a );
  yvec = r8vec_linspace_new ( ny, -a, a );

  xmat = new double[nx*ny];
  ymat = new double[nx*ny];
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
//
//  Sample OMEGA.
//
  seed = 123456789;
  omega = r8vec_normal_01_new ( m_1d * m_1d, seed );
//
//  Compute the diffusivity field.
//
  cl = 0.1;
  dc0 = 10.0;
  dc = diffusivity_2d_elman ( a, cl, dc0, m_1d, omega, nx, nx, xmat, ymat );
//
//  Create a data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      data_unit
        << "  " << setw(14) << xmat[i+j*nx]
        << "  " << setw(14) << ymat[i+j*nx]
        << "  " << setw(14) << dc[i+j*nx] << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str() );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'elman_contour.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set zlabel '<---DC(X,Y)--->'\n";
  command_unit << "set title 'Elman Stochastic diffusivity function'\n";
  command_unit << "set contour\n";
  command_unit << "set timestamp\n";
  command_unit << "set cntrparam levels 10\n";
  command_unit << "#set view map\n";
  command_unit << "set view 75, 75\n";
  command_unit << "unset key\n";
  command_unit << "splot '" << data_filename << "'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] dc;
  delete [] omega;
  delete [] xmat;
  delete [] xvec;
  delete [] ymat;
  delete [] yvec;

  return;
}
//****************************************************************************80

void ntw_contour ( )

//****************************************************************************80
//
//  Purpose:
//
//    NTW_CONTOUR displays a contour plot of a 2D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is compute by DIFFUSIVITY_2D_NTW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Fabio Nobile, Raul Tempone, Clayton Webster,
//    A Sparse Grid Stochastic Collocation Method for Partial Differential
//    Equations with Random Input Data,
//    SIAM Journal on Numerical Analysis,
//    Volume 46, Number 5, 2008, pages 2309-2345.
//
{
  double cl;
  string command_filename = "ntw_commands.txt";
  ofstream command_unit;
  double d;
  string data_filename = "ntw_data.txt";
  ofstream data_unit;
  double *dc;
  double dc0;
  int i;
  int j;
  int m = 21;
  int nx = 101;
  int ny = 101;
  double *omega;
  int seed;
  double *xmat;
  double *xvec;
  double *ymat;
  double *yvec;

  cout << "\n";
  cout << "NTW_CONTOUR\n";
  cout << "  Display contour or surface plots of the stochastic\n";
  cout << "  diffusivity function defined by DIFFUSIVITY_2D_NTW.\n";
//
//  Set the spatial grid.
//
  d = 1.0;
  xvec = r8vec_linspace_new ( nx, 0.0, d );
  yvec = r8vec_linspace_new ( ny, 0.0, d );

  xmat = new double[nx*ny];
  ymat = new double[nx*ny];
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
//
//  Sample OMEGA.
//  We rescale to  [-sqrt(3),sqrt(3)].
//
  seed = 123456789;
  omega = r8vec_uniform_01_new ( m, seed );
  for ( i = 0; i < m; i++ )
  {
    omega[i] = ( 1.0 - omega[i] ) * ( - sqrt ( 3.0 ) ) 
             +         omega[i]   *     sqrt ( 3.0 );
  }
//
//  Evaluate the diffusivity field.
//
  cl = 0.1;
  dc0 = 0.5;
  dc = diffusivity_2d_ntw ( cl, dc0, m, omega, nx * ny, xmat, ymat );
//
//  Create a data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      data_unit 
        << "  " << setw(14) << xmat[i+j*nx]
        << "  " << setw(14) << ymat[i+j*nx]
        << "  " << setw(14) << dc[i+j*nx] << "\n";
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
//
//  Create the command file.
//
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'ntw_contour.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set zlabel '<---DC(X,Y)--->'\n";
  command_unit << "set title 'NTW Stochastic diffusivity function'\n";
  command_unit << "set contour\n";
  command_unit << "set timestamp\n";
  command_unit << "set cntrparam levels 15\n";
  command_unit << "#set view map\n";
  command_unit << "set view 65, 65\n";
  command_unit << "set key\n";
  command_unit << "splot '" << data_filename << "'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'.\n";

  delete [] dc;
  delete [] omega;
  delete [] xmat;
  delete [] xvec;
  delete [] ymat;
  delete [] yvec;

  return;
}
//****************************************************************************80

void xk_contour ( )

//****************************************************************************80
//
//  Purpose:
//
//    XK_CONTOUR displays contour plots of a 1D stochastic diffusivity function.
//
//  Discussion:
//
//    The diffusivity function is compute by DIFFUSIVITY_1D_XK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongbin Xiu, George Karniadakis,
//    Modeling uncertainty in steady state diffusion problems via
//    generalized polynomial chaos,
//    Computer Methods in Applied Mechanics and Engineering,
//    Volume 191, 2002, pages 4927-4948.
//
{
  double *dc;
  double dc_max;
  double dc0;
  string command_filename = "xk_commands.txt";
  ofstream command_unit;
  string data_filename = "xk_data.txt";
  ofstream data_unit;
  int j;
  int m;
  int n;
  int seed;
  double *omega;
  double *x;
  double x_min;
  double x_max;

  cout << "\n";
  cout << "XK_CONTOUR\n";
  cout << "  Plot the stochastic diffusivity function\n";
  cout << "  defined by DIFFUSIVITY_1D_XK.\n";
//
//  Set up the spatial grid.
//
  n = 51;
  x_min = -1.0;
  x_max = +1.0;
  x = r8vec_linspace_new ( n, x_min, x_max );
//
//  Sample the OMEGA values.
//
  m = 5;
  seed = 123456789;
  omega = r8vec_normal_01_new ( m, seed );
//
//  Compute the diffusivity field.
//
  dc0 = 10.0;
  dc = diffusivity_1d_xk ( dc0, m, omega, n, x );
//
//  Create data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < n; j++ )
  {
    data_unit
      << "  " << setw(14) << x[j]
      << "  " << setw(14) << dc[j] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
//
//  Create the command file.
//
  dc_max = r8vec_max ( n, dc );

  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'xk_contour.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---DC(X)--->'\n";
  command_unit << "set yrange [0.0:" << dc_max << "]\n";
  command_unit << "set title 'XK Stochastic diffusivity function'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'red'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] dc;
  delete [] omega;
  delete [] x;

  return;
}
