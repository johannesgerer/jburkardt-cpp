# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <string>

using namespace std;

# include "burgers_solution.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BURGERS_SOLUTION_PRB.
//
//  Discussion:
//
//    BURGERS_SOLUTION_PRB tests the BURGERS_SOLUTION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "BURGERS_SOLUTION_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BURGERS_SOLUTION library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BURGERS_SOLUTION_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//.
//  Purpose:
//
//    TEST01 tests sets up a small test case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "burgers_test01.txt";
  double nu;
  double pi = 3.141592653589793;
  double thi;
  double tlo;
  double *vu;
  double *vt;
  int vtn = 11;
  double *vx;
  int vxn = 11;
  double xhi;
  double xlo;

  nu = 0.01 / pi;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Compute an analytic solution to the Burgers equation.\n";
  cout << "\n";
  cout << "  Viscosity NU = " << nu << "\n";
  cout << "  NX = " << vxn << "\n";
  cout << "  NT = " << vtn << "\n";

  xlo = -1.0;
  xhi = +1.0;
  vx = r8vec_even_new ( vxn, xlo, xhi );
  r8vec_print ( vxn, vx, "  X grid points:" );

  tlo = 0.0;
  thi = 3.0 / pi;
  vt = r8vec_even_new ( vtn, tlo, thi );
  r8vec_print ( vtn, vt, "  T grid points:" );

  vu = burgers_solution ( nu, vxn, vx, vtn, vt );

  r8mat_print ( vxn, vtn, vu, "  U(X,T) at grid points:" );

  r8mat_write ( filename, vxn, vtn, vu );

  cout << "\n";
  cout << "  Data written to file \"" << filename << "\".\n";

  delete [] vt;
  delete [] vu;
  delete [] vx;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests sets up a finer test case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "burgers_test02.txt";
  double nu;
  double pi = 3.141592653589793;
  double thi;
  double tlo;
  double *vu;
  double *vt;
  int vtn = 41;
  double *vx;
  int vxn = 41;
  double xhi;
  double xlo;

  nu = 0.01 / pi;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Compute an analytic solution to the Burgers equation.\n";
  cout << "\n";
  cout << "  Viscosity NU = " << nu << "\n";
  cout << "  NX = " << vxn << "\n";
  cout << "  NT = " << vtn << "\n";

  xlo = -1.0;
  xhi = +1.0;
  vx = r8vec_even_new ( vxn, xlo, xhi );
  r8vec_print ( vxn, vx, "  X grid points:" );

  tlo = 0.0;
  thi = 3.0 / pi;
  vt = r8vec_even_new ( vtn, tlo, thi );
  r8vec_print ( vtn, vt, "  T grid points:" );

  vu = burgers_solution ( nu, vxn, vx, vtn, vt );

  r8mat_write ( filename, vxn, vtn, vu );

  cout << "\n";
  cout << "  Data written to file \"" << filename << "\".\n";

  delete [] vt;
  delete [] vu;
  delete [] vx;

  return;
  return;
}
