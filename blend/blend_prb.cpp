# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "blend.hpp"

int main ( );
void cubic_rs ( double r, double s, int i, double *xi );
void identity_r ( double r, int i, double *xi );
void identity_rs ( double r, double s, int i, double *xi );
void identity_rst ( double r, double s, double t, int i, double *xi );
void quad_rst ( double r, double s, double t, int i, double *xi );
void stretch_r ( double r, int i, double *xi );
void stretch_rs ( double r, double s, int i, double *xi );
void stretch_rst ( double r, double s, double t, int i, double *xi );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BLEND_PRB.
//
//  Discussion:
//
//    BLEND_PRB tests routines from BLEND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "BLEND_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the BLEND library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "BLEND_PRB:\n";
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
//    TEST01 checks out BLEND_R_0DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double x[1];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Identity test on BLEND_R_0DN.\n";

  n = 1;

  r = 0.0;
  blend_r_0dn ( r, x, n, identity_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  r = 1.0;
  blend_r_0dn ( r, x, n, identity_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  r = 0.5;
  blend_r_0dn ( r, x, n, identity_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 checks out BLEND_RS_0DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double x[2];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Identity test on BLEND_RS_0DN.\n";

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.5;
  s = 0.5;
  blend_rs_0dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 checks out BLEND_RS_1DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double x[2];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Identity test on BLEND_RS_1DN.\n";

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.5;
  s = 0.5;
  blend_rs_1dn ( r, s, x, n, identity_rs );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 checks out BLEND_RST_0DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Identity test on BLEND_RST_0DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_0dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 checks out BLEND_RST_1DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Identity test on BLEND_RST_1DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_1dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 checks out BLEND_RST_2DN on the identity map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n" ;
  cout << "TEST06\n";
  cout << "  Identity test on BLEND_RST_2DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_2dn ( r, s, t, x, n, identity_rst );
  cout                    << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 checks out BLEND_R_0DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double x[1];

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Stretch test on BLEND_R_0DN.\n";

  n = 1;

  r = 0.0;
  blend_r_0dn ( r, x, n, stretch_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  r = 1.0;
  blend_r_0dn ( r, x, n, stretch_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  r = 0.5;
  blend_r_0dn ( r, x, n, stretch_r );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << x[0] << "\n";

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 checks out BLEND_RS_0DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double x[2];

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Stretch test on BLEND_RS_0DN.\n";

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 0.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 1.0;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.5;
  s = 0.5;
  blend_rs_0dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 checks out BLEND_RS_1DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double x[2];

  cout << "\n";
  cout << "TEST09\n";
  cout << "  Stretch test on BLEND_RS_1DN.\n";

  n = 2;

  r = 0.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 0.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 1.0;
  s = 1.0;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  r = 0.5;
  s = 0.5;
  blend_rs_1dn ( r, s, x, n, stretch_rs );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "\n";

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 checks out BLEND_RST_0DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n";
  cout << "TEST10\n";
  cout << "  Stretch test on BLEND_RST_0DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_0dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 checks out BLEND_RST_1DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  Stretch test on BLEND_RST_1DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_1dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 checks out BLEND_RST_2DN on the stretch map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double r;
  double s;
  double t;
  double x[3];

  cout << "\n";
  cout << "TEST12\n";
  cout << "  Stretch test on BLEND_RST_2DN.\n";

  n = 3;

  r = 0.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 0.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 1.0;
  t = 0.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.0;
  s = 0.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 1.0;
  s = 1.0;
  t = 1.0;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  r = 0.5;
  s = 0.5;
  t = 0.5;
  blend_rst_2dn ( r, s, t, x, n, stretch_rst );
  cout                     << "  "
       << setw(8) << r    << "  "
       << setw(8) << s    << "  "
       << setw(8) << t    << "  "
       << setw(8) << x[0] << "  "
       << setw(8) << x[1] << "  "
       << setw(8) << x[2] << "\n";

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 checks out BLEND_I_0D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int m;
  double x[5];

  m = 5;
  x[0] = 100.0;
  x[m-1] = 100.0 + ( m - 1 ) * 5;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  BLEND_I_0D1 interpolates data in a vector.\n";
  cout << "\n";
  cout << "  X[0] = " << x[0] << "\n";
  cout << "  X(" << m-1 << ")= " << x[m-1] << "\n";
  cout << "\n";
  cout << "  Interpolated values:\n";
  cout << "\n";

  blend_i_0d1 ( x, m );

  for ( i = 0; i < m; i++ )
  {
    cout                     << "  "
         << setw(6)  << i    << "  "
         << setw(8) << x[i] << "\n";
  }
  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 checks out BLEND_IJ_0D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m1 = 5;
  int m2 = 4;
  double r;
  double s;
  double temp;
  double x[20];

  cout << "\n" ;
  cout << "TEST14\n";
  cout << "  BLEND_IJ_0D1 interpolates data in a table,\n";
  cout << "  from corner data.\n";
  cout << "\n";
  cout << "  The table is " << m1 << " rows by " << m2 << " columns.\n";
//
//  Load data in the corners only.
//
  i = 0;
  j = 0;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = m1 - 1;
  j = 0;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = 0;
  j = m2 - 1;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  i = m1 - 1;
  j = m2 - 1;
  r = ( double ) i / ( double ) ( m1 - 1 );
  s = ( double ) j / ( double ) ( m2 - 1 );
  cubic_rs ( r, s, 1, &temp );
  x[i*m2+j] = temp;

  blend_ij_0d1 ( x, m1, m2 );

  cout << "\n";
  cout << "  Values interpolated by BLEND_IJ_0D1:\n";
  cout << "\n";

  for ( i = 0; i < m1; i++ )
  {
    cout                         << "  "
         << setw(8) << x[i*m2]   << "  "
         << setw(8) << x[i*m2+1] << "  "
         << setw(8) << x[i*m2+2] << "  "
         << setw(8) << x[i*m2+3] << "\n";
  }

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 checks out BLEND_IJ_1D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int m1 = 5;
  int m2 = 4;
  double r;
  double s;
  double temp;
  double x[20];

  cout << "\n";
  cout << "TEST15\n";
  cout << "  BLEND_IJ_1D1 interpolates data in a table,\n";
  cout << "  from edge data.\n";
  cout << "\n";
  cout << "  The table is " << m1 << " rows by " << m2 << " columns.\n";
//
//  Load data in the edges only.
//
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );

    j = 0;
    s = ( double ) j / ( double ) ( m2 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;

    j = m2 - 1;
    s = ( double ) j / ( double ) ( m2 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;
  }

  for ( j = 0; j < m2; j++ )
  {
    s = ( double ) j / ( double ) ( m2 - 1 );

    i = 0;
    r = ( double ) i / ( double ) ( m1 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;

    i = m1 - 1;
    r = ( double ) i / ( double ) ( m1 - 1 );
    cubic_rs ( r, s, 1, &temp );
    x[i*m2+j] = temp;
  }

  blend_ij_1d1 ( x, m1, m2 );

  cout << "\n";
  cout << "  Values interpolated by BLEND_IJ_1D1:\n";
  cout << "\n";

  for ( i = 0; i < m1; i++ )
  {
    cout                         << "  "
         << setw(8) << x[i*m2]   << "  "
         << setw(8) << x[i*m2+1] << "  "
         << setw(8) << x[i*m2+2] << "  "
         << setw(8) << x[i*m2+3] << "\n";
  }

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 checks out BLEND_IJK_0D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  int num_extreme;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  cout << "\n";
  cout << "TEST16\n";
  cout << "  BLEND_IJK_0D1 interpolates data in a 3D table,\n";
  cout << "  from corner data.\n";
  cout << "\n";
  cout << "  The table is " << m1 << " rows by " << m2 <<
    " columns by " << m3 << " layers.\n";
//
//  Load data in the faces only.
//
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        num_extreme = 0;
        if ( i == 0 || i == m1 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( j == 0 || j == m2 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( k == 0 || k == m3 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( num_extreme >= 3 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Data given to BLEND_IJK_0D1:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  blend_ijk_0d1 ( x, m1, m2, m3 );

  cout << "\n";
  cout << "  Values interpolated by BLEND_IJK_0D1:\n" ;
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }

  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Exact data:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 checks out BLEND_IJK_1D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  int num_extreme;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  cout << "\n";
  cout << "TEST17\n";
  cout << "  BLEND_IJK_1D1 interpolates data in a 3D table,\n";
  cout << "  from edge data.\n";
  cout << "\n";
  cout << "  The table is " << m1 << " rows by " << m2
    << " columns by " << m3 << " layers.\n";
//
//  Load data in the faces only.
//
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        num_extreme = 0;
        if ( i == 0 || i == m1 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( j == 0 || j == m2 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( k == 0 || k == m3 - 1 )
        {
          num_extreme = num_extreme + 1;
        }
        if ( num_extreme >= 2 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Data given to BLEND_IJK_1D1:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  blend_ijk_1d1 ( x, m1, m2, m3 );

  cout << "\n";
  cout << "  Values interpolated by BLEND_IJK_1D1:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }

  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Exact data:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n" ;

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 checks out BLEND_IJK_2D1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int k;
  int m1 = 4;
  int m2 = 3;
  int m3 = 3;
  double r;
  double s;
  double t;
  double temp;
  double x[36];

  cout << "\n";
  cout << "TEST18\n";
  cout << "  BLEND_IJK_2D1 interpolates data in a 3D table,\n";
  cout << "  from face data.\n";
  cout << "\n";
  cout << "  The table is " << m1 << " rows by " << m2
       << " columns by " << m3 << " layers.\n";
//
//  Load data in the faces only.
//
  for ( i = 0; i < m1; i++ )
  {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        if ( i == 0 || i == m1 - 1 ||
             j == 0 || j == m2 - 1 ||
             k == 0 || k == m3 - 1 )
        {
          quad_rst ( r, s, t, 1, &temp );
        }
        else
        {
          temp = 0.0;
        }
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Data given to BLEND_IJK_2D1:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                             << "  "
        << setw(8) << x[(i*m3+0)*m2+k] << "  "
        << setw(8) << x[(i*m3+1)*m2+k] << "  "
        << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  blend_ijk_2d1 ( x, m1, m2, m3 );

  cout << "\n";
  cout << "  Values interpolated by BLEND_IJK_2D1:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }

  for ( i = 0; i < m1; i++ )
   {
    r = ( double ) i / ( double ) ( m1 - 1 );
    for ( j = 0; j < m2; j++ )
    {
      s = ( double ) j / ( double ) ( m2 - 1 );
      for ( k = 0; k < m3; k++ )
      {
        t = ( double ) k / ( double ) ( m3 - 1 );
        quad_rst ( r, s, t, 1, &temp );
        x[(i*m3+j)*m2+k] = temp;
      }
    }
  }

  cout << "\n";
  cout << "  Exact data:\n";
  cout << "\n";

  for ( k = 0; k < m3; k++ )
  {
    cout << "\n";
    cout << "  Layer K = " << k << "\n";
    cout << "\n";

    for ( i = 0; i < m1; i++ )
    {
      cout                                << "  "
           << setw(8) << x[(i*m3+0)*m2+k] << "  "
           << setw(8) << x[(i*m3+1)*m2+k] << "  "
           << setw(8) << x[(i*m3+2)*m2+k] << "\n";
    }
  }
  return;
}
//****************************************************************************80

void cubic_rs ( double r, double s, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    CUBIC_RS evaluates a function of R and S used for some tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  *xi = 20.0 * ( r * r * s * s * s );

  return;
}
//****************************************************************************80

void quad_rst ( double r, double s, double t, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_RST evaluates a function of R, S and T used for some tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
{
  *xi = 18.0 * ( r * r + s + t );

  return;
}
//****************************************************************************80

void identity_r ( double r, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    IDENTITY_R returns a data component given (R).
//
//  Discussion:
//
//    This is a dummy routine, which simply returns (R).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the coordinate of a point that lies on the
//    boundary of the cube.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R), which is on an endpoint of the unit line segment.
//
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else
  {
    cout << "\n";
    cout << "IDENTITY_R - Fatal error!\n";
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
//****************************************************************************80

void identity_rs ( double r, double s, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    IDENTITY_RS returns a data component given (R,S).
//
//  Discussion:
//
//    This is a dummy routine, which simply returns (R,S).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the coordinates of a point that lies on the
//    boundary of the square.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R,S), which is on a corner, or edge, of the unit square.
//
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else if ( i == 1 )
  {
    *xi = s;
  }
  else
  {
    cout << "\n";
    cout << "IDENTITY_RS - Fatal error!\n" ;
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
//****************************************************************************80

void identity_rst ( double r, double s, double t, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    IDENTITY_RST returns a data component given (R,S,T).
//
//  Discussion:
//
//    This is a dummy routine, which simply returns (R,S,T).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, T, the coordinates of a point that lies on the
//    boundary of the cube.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R,S), which is on a corner, edge or face of the unit cube.
//
{
  if ( i == 0 )
  {
    *xi = r;
  }
  else if ( i == 1 )
  {
    *xi = s;
  }
  else if ( i == 2 )
  {
    *xi = t;
  }
  else
  {
    cout << "\n";
    cout << "IDENTITY_RST - Fatal error!\n";
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
//****************************************************************************80

void stretch_r ( double r, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    STRETCH_R returns a data component given (R).
//
//  Discussion:
//
//    This routine shifts by 1 and stretches by 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, the coordinate of a point that lies on the
//    boundary of the cube.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R), which is on an endpoint of the unit line segment.
//
{
  if ( i == 0 )
  {
    *xi = 2.0 * r + 1.0;
  }
  else
  {
    cout << "\n";
    cout << "STRETCH_R - Fatal error\n";
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
//****************************************************************************80

void stretch_rs ( double r, double s, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    STRETCH_RS returns a data component given (R,S).
//
//  Discussion:
//
//    This routine shifts by (1,2) and stretches by (3,4).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, the coordinates of a point that lies on the
//    boundary of the square.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R,S), which is on a corner, or edge, of the unit square.
//
{
  if ( i == 0 )
  {
    *xi = 3.0 * r + 1.0;
  }
  else if ( i == 1 )
  {
    *xi = 4.0 * s + 2.0;
  }
  else
  {
    cout << "\n";
    cout << "STRETCH_RS - Fatal error!\n";
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
//****************************************************************************80

void stretch_rst ( double r, double s, double t, int i, double *xi )

//****************************************************************************80
//
//  Purpose:
//
//    STRETCH_RST returns a data component given (R,S,T).
//
//  Discussion:
//
//    This routine shifts by (1,2,3) and stretches by (4,5,6)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 October 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double R, S, T, the coordinates of a point that lies on the
//    boundary of the cube.
//
//    Input, int I, the component of the data which is to be returned.
//
//    Output, double *XI, the I-th component of the data vector X, evaluated
//    at the point (R,S), which is on a corner, edge or face of the unit cube.
//
{
  if ( i == 0 )
  {
    *xi = 4.0 * r + 1.0;
  }
  else if ( i == 1 )
  {
    *xi = 5.0 * s + 2.0;
  }
  else if ( i == 2 )
  {
    *xi = 6.0 * t + 3.0;
  }
  else
  {
    cout << "\n";
    cout << "STRETCH_RST - Fatal error\n";
    cout << "  Illegal component index I = " << i << "\n";
    *xi = 0.0;
  }

  return;
}
