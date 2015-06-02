# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_svg.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE_SVG_PRB.
//
//  Discussion:
//
//    TRIANGLE_SVG_PRB tests the TRIANGLE_SVG library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE_SVG_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TRIANGLE_SVG library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE_SVG_PRB\n";
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
//    TEST01 calls TRIANGLE_SVG to plot a triangle and some points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 April 2014
//
//  Author:
//
//    John Burkardt
//
{ 
  double p[2*13] = {
    0.333333333333333, 0.333333333333333,
    0.479308067841923, 0.260345966079038,
    0.260345966079038, 0.479308067841923,
    0.260345966079038, 0.260345966079038,
    0.869739794195568, 0.065130102902216,
    0.065130102902216, 0.869739794195568,
    0.065130102902216, 0.065130102902216,
    0.638444188569809, 0.312865496004875,
    0.638444188569809, 0.048690315425316,
    0.312865496004875, 0.638444188569809,
    0.312865496004875, 0.048690315425316,
    0.048690315425316, 0.638444188569809,
    0.048690315425316, 0.312865496004875 };
  int p_num = 13;
  string plot_filename = "test01.svg";
  double t[2*3] = {
    0.0, 0.0, 
    1.0, 0.0, 
    0.0, 1.0 };

  triangle_svg ( plot_filename, t, p_num, p );

  return;
}
