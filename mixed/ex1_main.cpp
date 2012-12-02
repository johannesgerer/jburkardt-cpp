# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

extern "C"
{
  float double_add_func ( double *r1, double *r2 );
  void  double_add_sub ( double *r1, double *r2, double *r3 );
  int  int_add_func ( int *i1, int *i2 );
  void int_add_sub ( int *i1, int *i2, int *i3 );
  float real_add_func ( float *r1, float *r2 );
  void  real_add_sub ( float *r1, float *r2, float *r3 );
}

int main ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    This is a C++ main program, which will call a C subroutine.
//
//  Modified:
//
//    04 January 2006
//
//  Author:
//
//    John Burkardt
//
{
  double d1;
  double d2;
  double d3;
  int i1;
  int i2;
  int i3;
  float r1;
  float r2;
  float r3;

  cout << "\n";
  cout << "MAIN\n";
  cout << "\n";
  cout << "  Demonstrate how a C++ main program\n";
  cout << "  may call a C routine, and live!\n";
  cout << "\n";
  cout << "  The C++ main program will now call the C routine.\n";

  cout << "\n";
  cout << "MAIN:\n";
  cout << "  Set integers I1 and I2, then call the C routine \n";
  cout << "  INT_ADD_SUB to compute I3 = I1 + I2.\n";
  i1 = 42;
  i2 = 22;
  i3 = 0;
  cout << "\n";
  cout << "  I1 = " << i1 << "\n";
  cout << "  I2 = " << i2 << "\n";

  int_add_sub ( &i1, &i2, &i3 );

  cout << "\n";
  cout << "  INT_ADD_SUB returned I3 = " << i3 << "\n";

  cout << "\n";
  cout << "  Now repeat, but using a C function.\n";

  i3 = int_add_func ( &i1, &i2 );

  cout << "\n";
  cout << "  INT_ADD_FUNC returned I3 = " << i3 << "\n";

  cout << "\n";
  cout << "MAIN:\n";
  cout << "  Set floats R1 and R2, then call the C routine \n";
  cout << "  REAL_ADD_SUB to compute R3 = R1 + R2.\n";
  r1 = 42.0;
  r2 = 22.0;
  r3 = 0.0;

  real_add_sub ( &r1, &r2, &r3 );

  cout << "\n";
  cout << "  REAL_ADD_SUB returned R3 = " << r3 << "\n";

  cout << "\n";
  cout << "  Now repeat, but using a C function.\n";

  r3 = real_add_func ( &r1, &r2 );

  cout << "\n";
  cout << "  REAL_ADD_FUNC returned R3 = " << r3 << "\n";

  cout << "\n";
  cout << "MAIN:\n";
  cout << "  Set doubles D1 and D2, then call the C routine \n";
  cout << "  DOUBLE_ADD_SUB to compute D3 = D1 + D2.\n";
  d1 = 42.0;
  d2 = 22.0;
  d3 = 0.0;

  double_add_sub ( &d1, &d2, &d3 );

  cout << "\n";
  cout << "  DOUBLE_ADD_SUB returned D3 = " << d3 << "\n";

  cout << "\n";
  cout << "  Now repeat, but using a C function.\n";

  d3 = double_add_func ( &d1, &d2 );

  cout << "\n";
  cout << "  REAL_ADD_FUNC returned R3 = " << d3 << "\n";

  return 0;
}
