# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

//******************************************************************************

float f ( float t )

//******************************************************************************
//
//  Purpose:
//
//     F evaluates the power consumption function F(T).
//
//  Modified:
//
//    03 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rubin Landau,
//    A First Course in Scientific Computing,
//    Princeton, 2005,
//    ISBN: 0-691-12183-4
//    LC: Q183.9.L36
//
//  Parameters:
//
//    Input, float T, the argument of the function.
//
//    Output, float F, the value of the function.
//
{
  float pi = 3.14159265;
  float value;

  value = ( 4.0 + t / 365.0 + 0.5 * sin ( pi * t / 91.0 ) ) * 
          ( 2.0 + exp ( - sin ( 2.0 * pi * t ) ) );

  return value;
}
