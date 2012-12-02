# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cpp_simple.H"

//******************************************************************************

int main ( void )

//******************************************************************************
//
//  Purpose:
//
//    MAIN is the main program for the CPP_SIMPLE example.
//
//  Modified:
//
//    03 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Local, float A, B, the endpoints of the interval of integration.
//
//    Local, float F ( float T ), the name of the function to be integrated.
//
//    Local, int INT_NUM, the number of intervals to be used.
//
//    Local, float QUAD, the approximate value of the integral.
//
{
  float a;
  float b;
  int int_num;
  float quad;
  int test;

  cout << "\n";
  cout << "CPP_SIMPLE\n";
  cout << "  A simple C++ program to demonstrate\n";
  cout << "  the use of makefiles.\n";

  cout << "\n";
  cout << "  Estimate the integral from 0 to 1000, of\n";
  cout << "  F(T) = (4+T/365+1/2 sin(pi*T/91)) * (2+exp(-sin(2*pi*T)))\n";
  cout << "  a function which models daily power consumption.\n";
  cout << "\n";
  cout << "  quad = midpoint ( a, b, f, int_num )\n";
  cout << "  estimates the integral using the midpoint rule.\n";
  cout << "\n";
  cout << "  f ( t )\n";
  cout << "  evaluates the integrand.\n";
  cout << "\n";
  cout << "  Intervals   Estimate\n";
  cout << "\n";

  a = 0.0;
  b = 1000.0;
  int_num = 100;

  for ( test = 1; test <= 3; test++ )
  {
    quad = midpoint ( a, b, &f, int_num );

    cout << "  " << setw(8) << int_num
         << "  " << setw(14) << quad << "\n";

    int_num = int_num * 100;
  }
 
  cout << "\n";
  cout << "CPP_SIMPLE:\n";
  cout << "  Normal end of execution.\n";
 
  return 0;
}
