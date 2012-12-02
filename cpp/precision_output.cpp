# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>

using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PRECISION_OUTPUT.
//
//  Discussion:
//
//    In "iomanip", there is a manipulator function called "setprecision"
//    which allows you to print floating point numbers to a particular precision.
//
//    Unfortunately, once you've called setprecision, ALL floating point
//    numbers will thereafter be printed to that precision.  So it's not like
//    using a format in C, where you can specify the precision to be used for
//    a single output statement, without affecting any others.
//
//    Thus, it's helpful to use the "precision" function which can capture
//    the current precision setting, and then be used to restore it.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 March 2012
//
//  Author:
//
//    John Burkardt
//
{
  int prec;
  double x;
  double y;

  cout << "\n";
  cout << "PRECISION_OUTPUT:\n";
  cout << "  C++ version.\n";
  cout << "  Demonstrate the precision() function for controlling the default\n";
  cout << "  precision of floating point output.\n";

  prec = cout.precision ( );

  cout << "\n";
  cout << "  Save prec = cout.precision() = " << prec << "\n";

  x = sqrt ( 2.0 );
  y = sqrt ( 3.0 );

  cout << "\n";
  cout << "  Sqrt ( 2.0 ) = " << x << "\n";
  cout << "  Sqrt ( 3.0 ) = " << y << "\n";
  cout << "\n";
  cout << "  Use setprecision(16), but only on the Sqrt ( 2.0 ) statement.\n";
  cout << "  However, all floating point output will now have the new precision.\n";
  cout << "\n";
  cout << "  Sqrt ( 2.0 ) = " << setprecision(16) << x << "\n";
  cout << "  Sqrt ( 3.0 ) = " << y << "\n";
  cout << "\n";
  cout << "  Use cout.precision(10) to request 10 digits of precision.\n";
  cout.precision ( 10 );
  cout << "\n";
  cout << "  Sqrt ( 2.0 ) = " << x << "\n";
  cout << "  Sqrt ( 3.0 ) = " << y << "\n";
  cout << "\n";
  cout << "  Use cout.precision(prec) to restore default precision.\n";
  cout.precision ( prec );
  cout << "\n";
  cout << "  Sqrt ( 2.0 ) = " << x << "\n";
  cout << "  Sqrt ( 3.0 ) = " << y << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "PRECISION_OUTPUT\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
