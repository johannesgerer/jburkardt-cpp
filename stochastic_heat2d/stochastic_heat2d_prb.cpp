# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "stochastic_heat2d.hpp"

int main ( );
void test01 ( );
void test02 ( );
void boundary ( int nx, int ny, double x[], double y[], int n, double a[], 
  double rhs[] );
double test01_f ( double x, double y );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for STOCHASTIC_HEAT2D_PRB.
//
//  Discussion:
//
//    STOCHASTIC_HEAT2D_PRB tests the STOCHASTIC_HEAT2D library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "STOCHASTIC_HEAT2D_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the STOCHASTIC_HEAT2D library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "STOCHASTIC_HEAT2D_PRB:\n";
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
//    TEST01 plots a sample solution of a 2D stochastic diffusivity equation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  string command_filename = "solution_commands.txt";
  ofstream command_unit;
  string data_filename = "solution_data.txt";
  ofstream data_unit;
  int i;
  int j;
  int nx;
  int ny;
  double *omega;
  int seed;
  double *umat;
  double u_mean;
  double *xmat;
  double xmax;
  double xmin;
  double *xvec;
  double *ymat;
  double ymax;
  double ymin;
  double *yvec;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Consider the steady heat equation in the unit square,\n";
  cout << "  with 0 Dirichlet boundary conditions, \n";
  cout << "  and a heat source term F that is a Gaussian centered at (0.60,0.80).\n";
  cout << "\n";
  cout << "  Model the diffusivity coefficient as spatially varying,\n";
  cout << "  with a stochastic dependence on parameters OMEGA(1:4),\n";
  cout << "  as described in Babuska, Nobile, Tempone (BNT).\n";
  cout << "\n";
  cout << "  Compute and display the solution U for a given choice\n";
  cout << "  of the parameters OMEGA.\n";
//
//  Create the X and Y coordinate vectors.
//
  nx = 21;
  xmin = 0.0;
  xmax = 1.0;
  xvec = r8vec_linspace_new ( nx, xmin, xmax );

  ny = 21;
  ymin = 0.0;
  ymax = 1.0;
  yvec = r8vec_linspace_new ( ny, ymin, ymax );
//
//  Create the X and Y coordinate matrices.
//
  xmat = new double[nx * ny];
  ymat = new double[nx * ny];
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
//
//  Sample OMEGA:
//
  seed = 123456789;
  omega = r8vec_normal_01_new ( 4, seed );
  for ( i = 0; i < 4; i++ )
  {
    omega[i] = 2.0 * omega[i];
  }

  r8vec_print ( 4, omega, "  Sampled OMEGA values:" );
//
//  Solve the finite difference approximation to the steady 2D heat equation
//  for this set of OMEGA values.
//
  umat = stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f );
//
//  Create a data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      data_unit << "  " << xmat[i+j*nx]
                << "  " << ymat[i+j*nx]
                << "  " << umat[i+j*nx] << "\n";
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
  command_unit << "set output 'solution.png'\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Y--->'\n";
  command_unit << "set zlabel '<---U(X,Y)--->'\n";
  command_unit << "set title 'Sample Solution'\n";
  command_unit << "set contour\n";
  command_unit << "set timestamp\n";
  command_unit << "set cntrparam levels 10\n";
  command_unit << "set view 75, 75\n";
  command_unit << "unset key\n";
  command_unit << "splot '" << data_filename << "'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'\n";
//
//  Report the average value of U.
//
  u_mean = r8mat_mean ( nx, ny, umat );

  cout << "\n";
  cout << "  Mean value of U is " << u_mean << "\n";
//
//  Free memory.
//
  delete [] omega;
  delete [] umat;
  delete [] xmat;
  delete [] xvec;
  delete [] ymat;
  delete [] yvec;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 looks at mean temperature as a function of OMEGA(1) and OMEGA(2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  string command_filename = "umean_commands.txt";
  ofstream command_unit;
  string data_filename = "umean_data.txt";
  ofstream data_unit;
  int i;
  int j;
  int nx;
  int ny;
  double omega[4];
  double *omega1_mat;
  double omega1_max;
  double omega1_min;
  int omega1_num;
  double *omega1_vec;
  double *omega2_mat;
  double omega2_max;
  double omega2_min;
  int omega2_num;
  double *omega2_vec;
  double *umat;
  double *u_mean_mat;
  double u_mean_max;
  double *xmat;
  double xmax;
  double xmin;
  double *xvec;
  double *ymat;
  double ymax;
  double ymin;
  double *yvec;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Fix OMEGA(3) = 4, OMEGA(4) = 0, and\n";
  cout << "  examine dependence of average temperature on OMEGA(1) and OMEGA(2)\n";
  cout << "  over the range [-10,+10].\n";
//
//  Create the X and Y coordinate vectors.
//
  nx = 21;
  xmin = 0.0;
  xmax = 1.0;
  xvec = r8vec_linspace_new ( nx, xmin, xmax );

  ny = 21;
  ymin = 0.0;
  ymax = 1.0;
  yvec = r8vec_linspace_new ( ny, ymin, ymax );
//
//  Create the X and Y coordinate matrices.
//
  xmat = new double[nx * ny];
  ymat = new double[nx * ny];
  r8vec_mesh_2d ( nx, ny, xvec, yvec, xmat, ymat );
//
//  Create OMEGA1 and OMEGA2 vectors.
//
  omega1_num = 21;
  omega1_min = -10.0;
  omega1_max = +10.0;
  omega1_vec = r8vec_linspace_new ( omega1_num, omega1_min, omega1_max );

  omega2_num = 21;
  omega2_min = -10.0;
  omega2_max = +10.0;
  omega2_vec = r8vec_linspace_new ( omega2_num, omega2_min, omega2_max );
//
//  Create the OMEGA1 and OMEGA2 coordinate matrices.
//
  omega1_mat = new double[omega1_num * omega2_num];
  omega2_mat = new double[omega1_num * omega2_num];
  r8vec_mesh_2d ( omega1_num, omega2_num, omega1_vec, omega2_vec, omega1_mat, omega2_mat );
//
//  Set OMEGA(3) and OMEGA(4).
//
  omega[2] = 4.0;
  omega[3] = 0.0;

  cout << "\n";
  cout << "  Omega(3) fixed at " << omega[2] << "\n";
  cout << "  Omega(4) fixed at " << omega[3] << "\n";
//
//  Solve the finite difference approximation to the steady 2D heat equation,
//  and save the mean value of the solution, which is a slightly biased
//  estimate of the heat integral over the unit square.
//
  u_mean_mat = new double[omega1_num * omega2_num];

  for ( j = 0; j < omega2_num; j++ )
  {
    omega[1] = omega2_vec[j];
    for ( i = 0; i < omega1_num; i++ )
    {
      omega[0] = omega1_vec[i];
      umat = stochastic_heat2d ( omega, nx, ny, xvec, yvec, test01_f );
      u_mean_mat[i+j*omega1_num] = r8mat_mean ( nx, ny, umat );
      delete [] umat;
    }
  }
//-
//  Create a data file.
//
  data_unit.open ( data_filename.c_str ( ) );
  for ( j = 0; j < ny; j++ )
  {
    for ( i = 0; i < nx; i++ )
    {
      data_unit << "  " << omega1_mat[i+j*omega1_num]
                << "  " << omega2_mat[i+j*omega1_num]
                << "  " << u_mean_mat[i+j*omega1_num] << "\n";
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
  command_unit << "set output 'umean.png'\n";
  command_unit << "set xlabel '<---OMEGA1--->'\n";
  command_unit << "set ylabel '<---OMEGA2--->'\n";
  command_unit << "set zlabel '<---U_MEAN(OMEGA1,OMEGA2)--->'\n";
  command_unit << "set title 'Solution Mean as Function of Omega1, Omega2'\n";
  command_unit << "set contour\n";
  command_unit << "set timestamp\n";
  command_unit << "set cntrparam levels 10\n";
  command_unit << "set view 75, 75\n";
  command_unit << "unset key\n";
  command_unit << "splot '" << data_filename << "'\n";

  command_unit.close ( );

  cout << "  Created graphics command file '" << command_filename << "'\n";
//
//  Print the maximum value of the mean.
//
  u_mean_max = r8mat_max ( omega1_num, omega2_num, u_mean_mat );

  cout << "\n";
  cout << "  U_Mean_Max = " << u_mean_max << "\n";
//
//  Free memory.
//
  delete [] omega1_mat;
  delete [] omega1_vec;
  delete [] omega2_mat;
  delete [] omega2_vec;
  delete [] u_mean_mat;
  delete [] xmat;
  delete [] xvec;
  delete [] ymat;
  delete [] yvec;

  return;
}
//****************************************************************************80

void boundary ( int nx, int ny, double x[], double y[], int n, double a[], 
  double rhs[] )

//****************************************************************************80
//
//  Purpose:
//
//    BOUNDARY sets up the matrix and right hand side at boundary nodes.
//
//  Discussion:
//
//    For this simple problem, the boundary conditions specify that the solution
//    is 100 on the left side, and insulated on the right, top and bottom.
//
//    Nodes are assigned a single index K, which increases as:
//
//    (NY-1)*NX+1  (NY-1)*NX+2  ...  NY * NX
//           ....         ....  ...    .....
//           NX+1         NX+2  ...   2 * NX
//              1            2  ...       NX
//
//    The index K of a node on the lower boundary satisfies:
//      1 <= K <= NX
//    The index K of a node on the upper boundary satisfies:
//      (NY-1)*NX+1 <= K <= NY * NX
//    The index K of a node on the left boundary satisfies:
//      mod ( K, NX ) = 1
//    The index K of a node on the right boundary satisfies:
//      mod ( K, NX ) = 0
//
//    If we number rows from bottom I = 1 to top I = NY
//    and columns from left J = 1 to right J = NX, then the relationship
//    between the single index K and the row and column indices I and J is:
//      K = ( I - 1 ) * NX + J
//    and
//      J = 1 + mod ( K - 1, NX )
//      I = 1 + ( K - J ) / NX
//      
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NX, NY, the number of grid points in X and Y.
//
//    Input, double X[NX], Y[NY], the coordinates of grid lines.
//
//    Input, int N, the number of nodes.
//
//    Input/output, double A[N*N].  On input, the system matrix, with the 
//    entries for the interior nodes filled in.  On output, the entries for
//    the boundary nodes have been set as well.
//
//    Input, double RHS[N], on input, the system right hand side, 
//    with the entries for the interior nodes filled in.  On output, the entries for
//    the boundary nodes have been set as well.
//
{
  int i;
  int j;
  int kc;
//
//  Left boundary.
//
  j = 0;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
//
//  Right boundary.
//
  j = nx - 1;
  for ( i = 1; i < ny - 1; i++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
//
//  Lower boundary.
//
  i = 0;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }
//
//  Upper boundary.
//
  i = ny - 1;
  for ( j = 0; j < nx; j++ )
  {
    kc = i * nx + j;
    a[kc+kc*n] = a[kc+kc*n] + 1.0;
    rhs[kc] = 0.0;
  }

  return;
}
//****************************************************************************80

double test01_f ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01_F evaluates the heat source term.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 September 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the evaluation point.
//
//    Output, double TEST01_F, the value of the heat source term at (X,Y).
//
{
  double arg;
  double v;
  double value;

  v = 0.05;
  arg = ( pow ( x - 0.60, 2 ) + pow ( y - 0.80, 2 ) ) / pow ( v, 2 );
  value = 2000.0 * exp ( - arg );

  return value;
}
