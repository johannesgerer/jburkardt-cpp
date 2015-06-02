# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iostream>

using namespace std;

# include "ffmsh_io.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FFMSH_IO_PRB.
//
//  Discussion:
//
//    FFMSH_IO_PRB tests the FFMSH_IO library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "FFMSH_IO_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the FFMSH_IO library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "FFMSH_IO_PRB\n";
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
//    TEST01 gets the example data and prints it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e_l;
  int e_num;
  int *e_v;
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Get the example 2D data and print it.\n";
//
//  Get example sizes.
//
  ffmsh_2d_size_example ( v_num, e_num, t_num );
//
//  Print example sizes.
//
  ffmsh_2d_size_print ( "  Example Sizes:", v_num, e_num, t_num );
//
//  Allocate memory.
//
  v_xy = new double[2*v_num];
  v_l = new int[v_num];
  e_v = new int[2*e_num];
  e_l = new int[e_num];
  t_v = new int[3*t_num];
  t_l = new int[t_num];
//
//  Get example data.
//
  ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
//
//  Print example data.
//
  ffmsh_2d_data_print ( "  Example data:", v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );
//
//  Free memory.
//
  delete [] e_l;
  delete [] e_v;
  delete [] t_l;
  delete [] t_v;
  delete [] v_l;
  delete [] v_xy;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 writes the example data to a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e_l;
  int e_num;
  int *e_v;
  string filename = "output.msh";
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  Get the example 2D data and print it.\n";
//
//  Get example sizes.
//
  ffmsh_2d_size_example ( v_num, e_num, t_num );
//
//  Allocate memory.
//
  v_xy = new double[2*v_num];
  v_l = new int[v_num];
  e_v = new int[2*e_num];
  e_l = new int[e_num];
  t_v = new int[3*t_num];
  t_l = new int[t_num];
//
//  Get example data.
//
  ffmsh_2d_data_example ( v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
//
//  Write the data to a file.
//
  ffmsh_2d_write ( filename, v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );

  cout << "\n";
  cout << "  The data was written to '" << filename << "'\n";
//
//  Free memory.
//
  delete [] e_l;
  delete [] e_v;
  delete [] t_l;
  delete [] t_v;
  delete [] v_l;
  delete [] v_xy;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 gets the example data from a file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int *e_l;
  int e_num;
  int *e_v;
  string filename = "input.msh";
  int *t_l;
  int t_num;
  int *t_v;
  int *v_l;
  int v_num;
  double *v_xy;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  Read the example 2D data from a file.\n";
//
//  Read sizes.
//
  ffmsh_2d_size_read ( filename, v_num, e_num, t_num );
//
//  Print sizes.
//
  ffmsh_2d_size_print ( "  Example Sizes:", v_num, e_num, t_num );
//
//  Allocate memory.
//
  v_xy = new double[2*v_num];
  v_l = new int[v_num];
  e_v = new int[2*e_num];
  e_l = new int[e_num];
  t_v = new int[3*t_num];
  t_l = new int[t_num];
//
//  Read data.
//
  ffmsh_2d_data_read ( filename, v_num, e_num, t_num, v_xy, v_l, e_v, e_l, 
    t_v, t_l );
//
//  Print data.
//
  ffmsh_2d_data_print ( "  Example data:", v_num, e_num, t_num, v_xy, 
    v_l, e_v, e_l, t_v, t_l );
//
//  Free memory.
//
  delete [] e_l;
  delete [] e_v;
  delete [] t_l;
  delete [] t_v;
  delete [] v_l;
  delete [] v_xy;

  return;
}
