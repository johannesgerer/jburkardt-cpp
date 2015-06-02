# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "edge.hpp"

//****************************************************************************80

double fx1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX1 is the first 1D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX1, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx1_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX2 is the second 1D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX2, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx2_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx3 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX3 is the third 1D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX3, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx3_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx4 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX4 is the fourth 1D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX4, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx4_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx5 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX5 is 1D example #5, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX5, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx5_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx6 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX6 is 1D example #6, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX6, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx6_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fx7 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    FX7 is 1D example #7, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double FX7, the function value.
//
{
  double value;
  double *value_vec;
  double x_vec[1];

  x_vec[0] = x;
  value_vec = fx7_vec ( 1, x_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxy1 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FXY1 is the first 2D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  value_vec = fxy1_vec ( 1, x_vec, y_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxy2 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FXY2 is the second 2D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  value_vec = fxy2_vec ( 1, x_vec, y_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxy3 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FXY3 is the third 2D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  value_vec = fxy3_vec ( 1, x_vec, y_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxy4 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FXY4 is the fourth 2D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  value_vec = fxy4_vec ( 1, x_vec, y_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxy5 ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    FXY5 is the fifth 2D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  value_vec = fxy5_vec ( 1, x_vec, y_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double fxyz1 ( double x, double y, double z )

//****************************************************************************80
//
//  Purpose:
//
//    FXYZ1 is the first 3D example, scalar version.
//
//  Discussion:
//
//    This function allows the user a more convenient interface when
//    only a single input argument is supplied.  See FXYZ1_VEC for details.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
{
  double value;
  double *value_vec;
  double x_vec[1];
  double y_vec[1];
  double z_vec[1];

  x_vec[0] = x;
  y_vec[0] = y;
  z_vec[0] = z;
  value_vec = fxyz1_vec ( 1, x_vec, y_vec, z_vec );
  value = value_vec[0];
  delete [] value_vec;

  return value;
}
//****************************************************************************80

double *fx1_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX1_VEC is the first 1D example, vector version.
//
//  Discussion:
//
//    This is example 3.1 in the reference.
//
//    The function should be plotted over [-1.0,+1.0].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX1_VEC[N], the function values.
//
//  Local parameters:
//
//    Local, real STEEP, controls the steepness of the slope.
//    The default value is a moderate 5.  For a sharp rise, use 25 instead.  
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;
  const double steep = 5.0;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] < 0.0 )
    {
      f[i] = cos ( 3.0 * r8_pi * x[i] );
    }
    else if ( 0.0 <= x[i] )
    {
      f[i] = - 1.0 + 2.0 / ( 1.0 + 3.0 * exp ( - steep * ( 2.0 * x[i] - 1.0 ) ) );
    }
  }

  return f;
}
//****************************************************************************80

double *fx2_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX2_VEC is the second 1D example, vector version.
//
//  Discussion:
//
//    The function should be plotted over [-1,+1].
//
//    The "internal" coordinate range will be [-2.0,6.0*pi].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX2_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;
  double x2;

  f = new double[n];
//
//  Map from the convenient range [-1,+1] to the physical range [-2,6pi].
//
  for ( i = 0; i < n; i++ )
  {
    x2 = ( ( 1.0 - x[i] ) * ( - 2.0 )  
         + ( 1.0 + x[i] ) * 6.0 * r8_pi ) 
         /   2.0;

    if ( x2 < 0.0 )
    {
      f[i] = exp ( x2 );
    }
    else if ( 0.0 <= x2 && x2 < 3.0 * r8_pi / 2.0 )
    {
      f[i] = - exp ( - x2 );
    }
    else if ( 3.0 * r8_pi / 2.0 <= x2 )
    {
      f[i] = -1.5 * sin ( x2 );
    }
  }

  return f;
}
//****************************************************************************80

double *fx3_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX3_VEC is the third 1D example, vector version.
//
//  Discussion:
//
//    The function should be plotted over [-1.0,+1.0].
//
//    Internally, this range is mapped to [-3.0,+3.0].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX3_VEC[N], the function values.
//
{
  double *f;
  int i;
  double x2;

  f = new double[n];
//
//  Map from the convenient range [-1,+1] to the physical range [-3,+3].
//
  for ( i = 0; i < n; i++ )
  {
    x2 = ( ( 1.0 - x[i] ) * ( -3.0 )   
         + ( 1.0 + x[i] ) * ( +3.0 ) ) 
         /   2.0;

    if ( -2.0 <= x2 && x2 <= -1.0 )
    {
      f[i] = 1.0;
    }
    else if ( -0.5 <= x2 && x2 <= 0.5 )
    {
      f[i] = 0.5 + 4.0 * pow ( x2 + 0.5, 2 );
    }
    else if ( 1.25 <= x2 && 3.0 * x2 <= 7.0 )
    {
      f[i] = 3.0 * ( 2.0 - x2 );
    }
    else
    {
      f[i] = 0.0;
    }
  }

  return f;
}
//****************************************************************************80

double *fx4_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX4_VEC is the fourth 1D example, vector version.
//
//  Discussion:
//
//    The function should be plotted over [0.0,+1.0].
//
//    The function is continuous, but the derivative has a discontinuity at 0.5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX4_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;
  double x2;

  f = new double[n];
//
//  Convert from -1 <= x <= 1 to 0 <= x <= 1:
//
  for ( i = 0; i < n; i++ )
  {
    x2 = ( x[i] + 1.0 ) / 2.0;

    if ( x2 <= 0.5 )
    {
      f[i] = - ( x2 - 0.5 ) + sin ( 4.0 * r8_pi * x2 ) / 6.0;
    }
    else if ( 0.5 < x2 )
    {
      f[i] =   ( x2 - 0.5 ) + sin ( 4.0 * r8_pi * x2 ) / 6.0;
    }
  }

  return f;
}
//****************************************************************************80

double *fx5_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX5_VEC is 1D example #5, vector version.
//
//  Discussion:
//
//    The function should be plotted over [-1.0,+1.0].
//
//    The function actually has no discontinuities, but does have a
//    steep rise.  The local parameter S controls the steepness of the rise.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX5_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double steep = 20.0;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = tanh ( steep * x[i] );
  }

  return f;
}
//****************************************************************************80

double *fx6_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX6_VEC is 1D example #6, vector version.
//
//  Discussion:
//
//    This is example 2.1 in the reference.
//
//    The function should be plotted over [0.0,+1.0].
//
//    The function has a discontinuous first derivative at 1/2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the location of discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX6_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = sin ( 2.0 * r8_pi * x[i] ) / 6.0;

    if ( x[i] < 0.5 )
    {
      f[i] = f[i] - ( x[i] - 0.5 );
    }
    else
    {
      f[i] = f[i] + ( x[i] - 0.5 );
    }
  }

  return f;
}
//****************************************************************************80

double *fx7_vec ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX7_VEC is 1D example #7, vector version.
//
//  Discussion:
//
//    This is example 2.1 in the reference.
//
//    The function should be plotted over [0.0,+1.0].
//
//    The function has a discontinuous second derivative at 1/2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the location of discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], the arguments.
//
//    Output, double FX6_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = sin ( 2.0 * r8_pi * x[i] ) / 6.0;

    if ( x[i] < 0.5 )
    {
      f[i] = f[i] - 0.5 * pow ( x[i] - 0.5, 2 );
    }
    else
    {
      f[i] = f[i] + 0.5 * pow ( x[i] - 0.5, 2 );
    }
  }

  return f;
}
//****************************************************************************80

double *fxy1_vec ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXY1_VEC is the first 2D example, vector version.
//
//  Discussion:
//
//    This is example 4.1 in the reference.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the arguments.
//
//    Output, double FXY1_VEC[N], the function values.
//
{
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = x[i] * y[i] + cos ( 2.0 * r8_pi * x[i] * x[i] ) 
      - sin ( 2.0 * r8_pi * x[i] * x[i] );

    if ( 0.25 < x[i] * x[i] + y[i] * y[i] )
    {
      f[i] = f[i] + 10.0 * x[i] - 5.0;
    }
  }

  return f;
}
//****************************************************************************80

double *fxy2_vec ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXY2_VEC is the second 2D example, vector version.
//
//  Discussion:
//
//    This is example 4.2 in the reference.
//
//    It is known as the Shepp-Logan phantom.
//
//    It should be plotted on [-1,+1] x [-1,+1].
//
//    Note that the Archibald reference incorrectly describes the divisor
//    of x in the second ellipse as 0.06624, when it should be 0.6624.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Polynomial fitting for edge detection in irregularly sampled signals 
//    and images,
//    SIAM Journal on Numerical Analysis,
//    Volume 43, Number 1, 2006, pages 259-279.
//
//    Larry Shepp, Ben Logan,
//    The Fourier reconstruction of a head section,
//    IEEE Transactions on Nuclear Science,
//    Volume  NS-21, June 1974, pages 21-43.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the arguments.
//
//    Output, double FXY2_VEC[N], the function values.
//
//  Local parameters:
//
//    Local, integer CHOICE:
//    1, use Archibald's (and Shepp and Logan's) level values;
//    2, use Matlab's level values;
//    3, use Matlab's enhanced contrast level values.
//
{
  double *c;
  double c1[4] = {
    2.0, -0.98, -0.02, +0.01 };
  double c2[4] = { 
    1.0, -0.98, -0.02, +0.01 };
  double c3[4] = {
    1.0, -0.8, -0.2, +0.1 }; 
  int choice;
  double eta1;
  double eta2;
  double *f;
  int i;
  const double r8_pi = 3.141592653589793;
  double xi1;
  double xi2;

  f = new double[n];

  choice = 3;

  if ( choice == 1 )
  {
    c = r8vec_copy_new ( 4, c1 );
  }
  else if ( choice == 2 )
  {
    c = r8vec_copy_new ( 4, c2 );
  }
  else
  {
    c = r8vec_copy_new ( 4, c3 );
  }

  for ( i = 0; i < n; i++ )
  {
    f[i] = 0.0;

    xi1  =   ( x[i] - 0.22 ) * cos ( 0.4 * r8_pi ) 
             + y[i]          * sin ( 0.4 * r8_pi );
    eta1 = - ( x[i] - 0.22 ) * sin ( 0.4 * r8_pi ) 
             + y[i]          * cos ( 0.4 * r8_pi );

    xi2  =   ( x[i] + 0.22 ) * cos ( 0.6 * r8_pi ) 
             + y[i]          * sin ( 0.6 * r8_pi );
    eta2 = - ( x[i] + 0.22 ) * sin ( 0.6 * r8_pi ) 
             + y[i]          * cos ( 0.6 * r8_pi );

    if ( pow ( x[i] / 0.69, 2 ) + pow ( y[i] / 0.92, 2 ) <= 1.0 )
    {
      f[i] = f[i] + c[0];
    }

    if ( pow ( x[i] / 0.6624, 2 ) 
       + pow ( ( y[i] + 0.0184 ) / 0.874, 2 ) <= 1.0 )
    {
      f[i] = f[i] + c[1];
    }

    if ( pow ( xi1 / 0.31, 2 ) + pow ( eta1 / 0.11, 2 ) <= 1.0 || 
         pow ( xi2 / 0.41, 2 ) + pow ( eta2 / 0.16, 2 ) <= 1.0 )
    {
      f[i] = f[i] + c[2];
    }

    if ( ( pow ( ( x[i] - 0.35 )  / 0.3, 2 )
         + pow (   y[i]               / 0.6, 2 ) <= 1.0 ) || 
         ( pow (   x[i]               / 0.21, 2 ) 
         + pow ( ( y[i] - 0.35  ) / 0.25, 2 ) <= 1.0 ) || 
         ( pow (   x[i]               / 0.046, 2 ) 
         + pow ( ( y[i] - 0.1   ) / 0.046, 2 ) <= 1.0 ) || 
         ( pow (   x[i]               / 0.046, 2 ) 
         + pow ( ( y[i] + 0.1   ) / 0.046, 2 ) <= 1.0 ) || 
         ( pow ( ( x[i] + 0.08  ) / 0.046, 2 ) 
         + pow ( ( y[i] + 0.605 ) / 0.023, 2 ) <= 1.0 ) || 
         ( pow (   x[i]               / 0.023, 2 ) 
         + pow ( ( y[i] + 0.605 ) / 0.023, 2 ) <= 1.0 ) || 
         ( pow ( ( x[i] - 0.06  ) / 0.023, 2 ) 
         + pow ( ( y[i] + 0.605 ) / 0.023, 2 ) <= 1.0 ) )
    {
      f[i] = f[i] + c[3];
    }
  }

  delete [] c;

  return f;
}
//****************************************************************************80

double *fxy3_vec ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXY3_VEC is the third 2D example, vector version.
//
//  Discussion:
//
//    This is example 3.2 in the reference.
//
//    It is known as the modified two-dimensional Harten function.
//
//    It should be plotted on [-1,+1] x [-1,+1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the locations and discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the arguments.
//
//    Output, double FXY3_VEC[N], the function values.
//
//  Local parameters:
//
//    Local, integer CHOICE:
//    1, use Archibald's (and Shepp and Logan's) level values;
//    2, use Matlab's level values;
//    3, use Matlab's enhanced contrast level values.
//
{
  double *f;
  int i;
  double r;
  const double r8_pi = 3.141592653589793;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    r = ( 4.0 * x[i] * x[i] + 4.0 * y[i] * y[i] - 1.0 ) / 6.0;

    if ( 3.0 * r <= -1.0 )
    {
      f[i] = - r * sin ( 0.5 * r8_pi * r * r );
    }
    else if ( 3.0 * r < 1.0 )
    {
      f[i] = fabs ( sin ( 2.0 * r8_pi * r ) );
    }
    else
    {
      f[i] = 2.0 * r - 1.0 - sin ( 3.0 * r8_pi * r ) / 6.0;
    }
  }

  return f;
}
//****************************************************************************80

double *fxy4_vec ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXY4_VEC is the fourth 2D example, vector version.
//
//  Discussion:
//
//    This is example 3.1 in the reference.
//
//    It is known as the discontinuous medium wave function.
//
//    Here, we are computing the first component of the solution, P(X,Y).
//
//    It should be plotted on (x,y) in [-1,0]x[0,0.1].
//
//    The second variable y actually represents time.
//
//    Note that in the reference, the formula reads:
//     f(i) = 2.0D+00 * rhor * cr / ( rhol * cl + rhor * cr ) 
//          * sin ( r8_pi * omega * ( y(i) - ( x(i) + 0.5D+00 ) / cr ) )
//    but I believe this should be:
//     f(i) = 2.0D+00 * rhor * cr / ( rhol * cl + rhor * cr ) 
//          * sin ( r8_pi * omega * ( y(i) - ( x(i) + 0.5D+00 ) / cl ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the locations and discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the arguments.
//
//    Output, double FXY3_VEC[N], the function values.
//
{
  const double cl = 0.87879;
  const double cr = 1.0;
  double *f;
  int i;
  const double omega = 12.0;
  const double r8_pi = 3.141592653589793;
  const double rhol = 0.55556;
  const double rhor = 1.0;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= -0.5 )
    {
      f[i] = sin ( r8_pi * omega * ( y[i] - ( x[i] + 0.5 ) / cl ) ) 
           - ( rhol * cl - rhor * cr ) / ( rhol * cl + rhor * cr )      
           * sin ( r8_pi * omega * ( y[i] + ( x[i] + 0.5 ) / cl ) );
    }
    else
    {
      f[i] = 2.0 * rhor * cr / ( rhol * cl + rhor * cr ) 
           * sin ( r8_pi * omega * ( y[i] - ( x[i] + 0.5 ) / cl ) );
    }
  }

  return f;
}
//****************************************************************************80

double *fxy5_vec ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXY5_VEC is the fifth 2D example, vector version.
//
//  Discussion:
//
//    This is example 3.1 in the reference.
//
//    It is known as the discontinuous medium wave function.
//
//    Here, we are computing the second component of the solution, U(X,Y).
//
//    It should be plotted on (x,y) in [-1,0]x[0,0.1].
//
//    The second variable y actually represents time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Rick Archibald, Anne Gelb, Jungho Yoon,
//    Determining the locations and discontinuities in the derivatives
//    of functions,
//    Applied Numerical Mathematics,
//    Volume 58, 2008, pages 577-592.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the arguments.
//
//    Output, double FXY3_VEC[N], the function values.
//
{
  const double cl = 0.87879;
  const double cr = 1.0;
  double *f;
  int i;
  const double omega = 12.0;
  const double r8_pi = 3.141592653589793;
  const double rhol = 0.55556;
  const double rhor = 1.0;

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( x[i] <= -0.5 )
    {
      f[i] = sin ( r8_pi * omega * ( y[i] - ( x[i] + 0.5 ) / cl ) ) 
           + ( rhol * cl - rhor * cr ) / ( rhol * cl + rhor * cr ) 
           / ( rhol * cl ) 
           * sin ( r8_pi * omega * ( y[i] + ( x[i] + 0.5 ) / cl ) );
    }
    else
    {
      f[i] = 2.0 / ( rhol * cl + rhor * cr ) 
           * sin ( r8_pi * omega * ( y[i] - ( x[i] + 0.5 ) / cl ) );
    }
  }

  return f;
}
//****************************************************************************80

double *fxyz1_vec ( int n, double x[], double y[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXYZ1_VEC is the first 3D example, vector version.
//
//  Discussion:
//
//    This example is known as the 3D Shepp-Logan phantom.
//
//    It should be plotted on [-1,+1] x [-1,+1.5] x [-1.5,+1.5].
//
//    Seventeen objects are modeled by ellipses of various gray levels,
//    including:
//
//     1: Outer skull
//     2: Inner skull
//     3: Left eye
//     4: Right eye
//     5: Nose
//     6: Mouth
//     7: Left ear
//     8: Right ear
//     9: Left small tumor
//    10: Center small tumor
//    11: Right small tumor
//    12: Old f
//    13: Old g
//    14: Old e
//    15: Right ventricle
//    16: Left ventricle
//    17: Blood clot
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 February 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Larry Shepp,
//    Computerized tomography and nuclear magnetic resonance,
//    Journal of Computer Assisted Tomography,
//    Volume 4, Number 1, February 1980, pages 94-107.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], Z[N], the arguments.
//
//    Output, double FXYZ1_VEC[N], the function values.
//
{
  double a1[17] = { 
          0.7233,  0.7008,  0.1270,  0.1270,  0.1270, 
          0.4575,  0.0635,  0.0635,  0.0460,  0.0230, 
          0.0230,  0.0460,  0.2100,  0.1100,  0.1600, 
          0.1600,  0.0300 };
  double a2[17] = { 
          0.9644,  0.9246,  0.1270,  0.1270,  0.3400, 
          0.6099,  0.3175,  0.3175,  0.0230,  0.0230, 
          0.0460,  0.0460,  0.2581,  0.2500,  0.3100, 
          0.4100,  0.2000 };
  double a3[17] = {
          1.2700,  1.2241,  0.1270,  0.1270,  0.1700, 
          0.5080,  0.3175,  0.3175,  0.0230,  0.0460, 
          0.0230,  0.0460,  0.2581,  0.2300,  0.2540, 
          0.3810,  0.2000 };
  double c;
  double *f;
  double g[17] = {
          2.0000, -0.9800, -1.0000, -1.0000,  1.5000, 
         -1.0000,  1.0000,  1.0000,  0.0100,  0.0100, 
          0.0100,  0.0100,  0.0100,  0.0100, -0.0200, 
         -0.0200,  0.0300 };
  int e;
  int i;
  double v11[17] = { 
          1.0000,  1.0000,  1.0000,  1.0000,  1.0000, 
          1.0000,  0.9903, -0.9903,  1.0000,  1.0000, 
          1.0000,  1.0000,  1.0000,  1.0000,  0.9511, 
         -0.9511,  0.9192 };
  double v12[17] = { 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000, -0.1085, -0.1085,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000, -0.3090, 
         -0.3090, -0.3381 };
  double v13[17] = { 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000, -0.0865, -0.0865,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.2020 };
  double v21[17] = { 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.1089, -0.1089,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.3090, 
         -0.3090,  0.3452 };
  double v22[17] = {
          1.0000,  1.0000,  1.0000,  1.0000,  0.5446, 
          1.0000,  0.9941,  0.9941,  1.0000,  1.0000, 
          1.0000,  1.0000,  1.0000,  1.0000,  0.9511, 
          0.9511,  0.9385 };
  double v23[17] = { 
          0.0000,  0.0000,  0.0000,  0.0000, -0.8387, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.0000 };
  double v31[17] = { 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.0860, -0.0860,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000,  0.1896 };
  double v32[17] = {
          0.0000,  0.0000,  0.0000,  0.0000,  0.8387, 
          0.0000, -0.0094, -0.0094,  0.0000,  0.0000, 
          0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 
          0.0000, -0.0697 };
  double v33[17] = {
          1.0000,  1.0000,  1.0000,  1.0000,  0.5446, 
          1.0000,  0.9963,  0.9963,  1.0000,  1.0000, 
          1.0000,  1.0000,  1.0000,  1.0000,  1.0000, 
          1.0000, -0.9794 };
  double x0[17] = {
          0.0000,  0.0000,  0.2583, -0.2583,  0.0000, 
          0.0000,  0.7076, -0.7076, -0.0800,  0.0000, 
          0.0600,  0.0000,  0.0000,  0.0000,  0.2200, 
         -0.2200,  0.5600 };
  double y0[17] = { 
          0.0000, -0.0184,  0.7534,  0.7534,  1.1398, 
          0.0000, -0.1378, -0.1378, -0.6050, -0.6050, 
         -0.6050,  0.1000, -0.1000,  0.3500,  0.0000, 
          0.0000, -0.4000 };
  double z0[17] = {
          0.0000, -0.0185,  0.0000,  0.0000, -0.1957, 
         -0.7620, -0.1905, -0.1905,  0.3810,  0.3810, 
          0.3810,  0.3810,  0.1270,  0.3810,  0.3810, 
          0.3810,  0.3810 };

  f = new double[n];

  for ( i = 0; i < n; i++ )
  {
    f[i] = 0.0;

    for ( e = 0; e < 17; e++ ) 
    {
      c = pow ( ( ( x[i] - x0[e] ) * v11[e] 
                + ( y[i] - y0[e] ) * v12[e] 
                + ( z[i] - z0[e] ) * v13[e] ) / a1[e], 2 )
        + pow ( ( ( x[i] - x0[e] ) * v21[e] 
                + ( y[i] - y0[e] ) * v22[e] 
                + ( z[i] - z0[e] ) * v23[e] ) / a2[e], 2 )
        + pow ( ( ( x[i] - x0[e] ) * v31[e] 
                + ( y[i] - y0[e] ) * v32[e] 
                + ( z[i] - z0[e] ) * v33[e] ) / a3[e], 2 );

      if ( c <= 1.0 )
      {
        f[i] = f[i] + g[e];
      }

    }

  }

  return f;
}
//****************************************************************************80

double *r8vec_copy_new ( int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY_NEW copies an R8VEC to a new R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double R8VEC_COPY_NEW[N], the copy of A1.
//
{
  double *a2;
  int i;

  a2 = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return a2;
}
//****************************************************************************80

double *r8vec_linspace_new ( int n, double a_first, double a_last )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
//
//    In other words, the interval is divided into N-1 even subintervals,
//    and the endpoints of intervals are used as the points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first 
             + ( double ) (         i ) * a_last ) 
             / ( double ) ( n - 1     );
    }
  }
  return a;
}
//****************************************************************************80

double *r8vec_uniform_ab_new ( int n, double a, double b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_AB_NEW returns a scaled pseudorandom R8VEC.
//
//  Discussion:
//
//    Each dimension ranges from A to B.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A, B, the lower and upper limits of the pseudorandom values.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_AB_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_AB_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = a + ( b - a ) * ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
