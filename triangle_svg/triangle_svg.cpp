# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "triangle_svg.hpp"

//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  }
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

int r8_to_i4 ( double xmin, double xmax, double x, int ixmin, int ixmax )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_I4 maps real X in [XMIN, XMAX] to integer IX in [IXMIN, IXMAX].
//
//  Discussion:
//
//    IX := IXMIN + ( IXMAX - IXMIN ) * ( X - XMIN ) / ( XMAX - XMIN )
//    IX := min ( IX, max ( IXMIN, IXMAX ) )
//    IX := max ( IX, min ( IXMIN, IXMAX ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double XMIN, XMAX, the real range.  XMAX and XMIN must not be
//    equal.  It is not necessary that XMIN be less than XMAX.
//
//    Input, double X, the real number to be converted.
//
//    Input, int IXMIN, IXMAX, the allowed range of the output
//    variable.  IXMAX corresponds to XMAX, and IXMIN to XMIN.
//    It is not necessary that IXMIN be less than IXMAX.
//
//    Output, int R8_TO_I4, the value in the range [IXMIN,IXMAX] that
//    corresponds to X.
//
{
  int ix;
  double temp;

  if ( xmax == xmin )
  {
    cerr << "\n";
    cerr << "R8_TO_I4 - Fatal error!\n";
    cerr << "  XMAX = XMIN, making a zero divisor.\n";
    cerr << "  XMAX = " << xmax << "\n";
    cerr << "  XMIN = " << xmin << "\n";
    exit ( 1 );
  }

  temp =
      ( ( xmax - x        ) * ( double ) ixmin
      + (        x - xmin ) * ( double ) ixmax )
      / ( xmax     - xmin );

  if ( 0.0 <= temp )
  {
    temp = temp + 0.5;
  }
  else
  {
    temp = temp - 0.5;
  }

  ix = ( int ) temp;

  return ix;
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
//****************************************************************************80

void triangle_svg ( string plot_filename, double t[], int p_num, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SVG plots a triangle and points in SVG format.
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
//  Parameters:
//
//    Input, string PLOT_FILENAME, the name of the output file.
//
//    Input, double T[2*3], points forming a triangle.
//
//    Input, int P_NUM, the number of points.
//
//    Input, double P[2*P_NUM], the points.
//
{
  int i;
  int i4;
  int i4_max;
  int i4_min;
  int ii;
  int j;
  int j4;
  int j4_max;
  int j4_min;
  int node;
  ofstream output;
  int r;
  double x;
  double x_max;
  double x_min;
  double x_scale;
  double y;
  double y_max;
  double y_min;
  double y_scale;
//
//  Determine SCALE, the maximum data range.
//
  x_max = p[0+0*2];
  x_min = p[0+0*2];
  for ( j = 0; j < p_num; j++ )
  {
    x_max = r8_max ( x_max, p[0+j*2] );
    x_min = r8_min ( x_min, p[0+j*2] );
  }
  for ( j = 0; j < 3; j++ )
  {
    x_max = r8_max ( x_max, t[0+j*2] );
    x_min = r8_min ( x_min, t[0+j*2] );
  }
  x_scale = x_max - x_min;
  x_max = x_max + 0.05 * x_scale;
  x_min = x_min - 0.05 * x_scale;
  x_scale = x_max - x_min;

  y_max = p[1+0*2];
  y_min = p[1+0*2];
  for ( j = 0; j < p_num; j++ )
  {
    y_max = r8_max ( y_max, p[1+j*2] );
    y_min = r8_min ( y_min, p[1+j*2] );
  }
  for ( j = 0; j < 3; j++ )
  {
    y_max = r8_max ( y_max, t[1+j*2] );
    y_min = r8_min ( y_min, t[1+j*2] );
  }
  y_scale = y_max - y_min;
  y_max = y_max + 0.05 * y_scale;
  y_min = y_min - 0.05 * y_scale;
  y_scale = y_max - y_min;

  i4_min = 1;
  j4_min = 1;
  if ( x_scale < y_scale )
  {
    i4_max = ( int ) ( 0.5 + 500.0 * x_scale / y_scale );
    j4_max = 500;
  }
  else
  {
    i4_max = 500;
    j4_max = ( int ) ( 0.5 + 500.0 * y_scale / x_scale );
  }
//
//  Open the file.
//
  output.open ( plot_filename.c_str ( ) );
//
//  Write that junk.
//
  output << "<?xml version = \"1.0\" standalone=\"no\"?>\n";
  output << "\n";
  output << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n";
  output << "  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
  output << "\n";
  output << "<svg\n";
  output << "  width=\"" << i4_max << "\"\n";
  output << "  height=\"" << j4_max << "\"\n";
  output << "  viewbox=\"" << i4_min 
         << "," << j4_min 
         << "," << i4_max 
         << "," << j4_max << "\"\n";
  output << "  xmlns=\"http://www.w3.org/2000/svg\"\n";
  output << "  version=\"1.1\">\n";
  output << "  <desc>\n";
  output << "    Triangulation created by triangle_svg.c\n";
  output << "  </desc>\n";
//
//  Draw the triangle.
//
  output << "  <polygon\n";
  output << "    fill=\"pink\"\n";
  output << "    stroke=\"black\"\n";
  output << "    stroke-width=\"2\"\n";
  output << "    points=\"\n";

  for ( j = 0; j < 3; j++ )
  {
    i4 = r8_to_i4 ( x_min, x_max, t[0+j*2], i4_min, i4_max );
    j4 = r8_to_i4 ( y_max, y_min, t[1+j*2], j4_min, j4_max );
    output << "      "<< i4 << "," << j4 << "\n";
  }
  output << "  \" />\n";
//
//  Draw points.
//
  for ( j = 0; j < p_num; j++ )
  {
    i4 = r8_to_i4 ( x_min, x_max, p[0+j*2], i4_min, i4_max );
    j4 = r8_to_i4 ( y_max, y_min, p[1+j*2], j4_min, j4_max );
    r = 5;

    output << "  <circle\n";
    output << "    cx=\"" << i4 << "\"\n";
    output << "    cy=\"" << j4 << "\"\n";
    output << "    r=\"" << r << "\"\n";
    output << "    fill=\"blue\"\n";
    output << "    stroke=\"black\"\n";
    output << "    stroke-width=\"2\"\n";
    output << "  />\n";
  }
//
//  End of plot.
//
  output << "</svg>\n";

  output.close ( );

  cout << "\n";
  cout << "  Graphics data written to file \"" << plot_filename << "\"\n";

  return;
}
