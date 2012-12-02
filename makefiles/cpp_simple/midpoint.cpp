# include <cstdlib>
# include <iostream>
# include <cmath>

using namespace std;

//******************************************************************************

float midpoint ( float a, float b, float f ( float t ), int int_num )

//******************************************************************************
//
//  Purpose:
//
//    MIDPOINT approximates an integral using the composite midpoint rule.
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
//    Input, float A, B, the endpoints of the interval of integration.
//
//    Input, float F ( float T ), the name of the function to be integrated.
//
//    Input, int INT_NUM, the number of intervals to be used.
//
//    Output, float QUAD, the approximate value of the integral.
//
{
  int i;
  float int_width;
  float quad;
  float t;

  quad = 0.0;

  for ( i = 1; i <= int_num; i++ )
  { 
    t = ( ( float ) ( 2 * int_num - 2 * i + 1 ) * a   
        + ( float ) (               2 * i - 1 ) * b ) 
        / ( float ) ( 2 * int_num );

    quad = quad + f ( t );
  }

  int_width = ( b - a ) / float ( int_num );

  quad = quad * int_width;

  return quad;
}
