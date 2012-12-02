# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( );
bool ppma_write ( string file_out_name, int xsize, int ysize, int *r, 
  int *g, int *b );
bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *r,
  int *g, int *b );
bool ppma_write_header ( ofstream &file_out, string file_out_name, int xsize, 
  int ysize, int rgb_max );
void timestamp ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MANDELBROT.
//
//  Discussion:
//
//    MANDELBROT computes an image of the Mandelbrot set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 July 2010
//
//  Author:
//
//    John Burkardt
//
//  Local Parameters:
//
//    Local, integer COUNT_MAX, the maximum number of iterations taken
//    for a particular pixel.
//
{
  int *b;
  int c;
  int c_max;
  int *count;
  int count_max = 400;
  string filename = "mandelbrot.ppm";
  int *g;
  int i;
  bool ierror;
  int j;
  int k;
  int n = 501;
  int *r;
  double x;
  double x_max =   1.25;
  double x_min = - 2.25;
  double x1;
  double x2;
  double y;
  double y_max =   1.75;
  double y_min = - 1.75;
  double y1;
  double y2;

  timestamp ( );
  cout << "\n";
  cout << "MANDELBROT\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Create an ASCII PPM image of the Mandelbrot set.\n";
  cout << "\n";
  cout << "  For each point C = X + i*Y\n";
  cout << "  with X range [" << x_min << "," << x_max << "]\n";
  cout << "  and  Y range [" << y_min << "," << y_max << "]\n";
  cout << "  carry out " << count_max << " iterations of the map\n";
  cout << "  Z(n+1) = Z(n)^2 + C.\n";
  cout << "  If the iterates stay bounded (norm less than 2)\n";
  cout << "  then C is taken to be a member of the set.\n";
  cout << "\n";
  cout << "  An ASCII PPM image of the set is created using\n";
  cout << "    N = " << n << " pixels in the X direction and\n";
  cout << "    N = " << n << " pixels in the Y direction.\n";
//
//  Carry out the iteration for each pixel, determining COUNT.
//
  count = new int[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      x = ( ( double ) (     j     ) * x_max   
          + ( double ) ( n - j - 1 ) * x_min ) 
          / ( double ) ( n     - 1 );

      y = ( ( double ) (     i     ) * y_max   
          + ( double ) ( n - i - 1 ) * y_min ) 
          / ( double ) ( n     - 1 );

      count[i+j*n] = 0;

      x1 = x;
      y1 = y;

      for ( k = 1; k <= count_max; k++ )
      {
        x2 = x1 * x1 - y1 * y1 + x;
        y2 = 2 * x1 * y1 + y;

        if ( x2 < -2.0 || 2.0 < x2 || y2 < -2.0 || 2.0 < y2 )
        {
          count[i+j*n] = k;
          break;
        }
        x1 = x2;
        y1 = y2;
      }
    }
  }
//
//  Determine the coloring of each pixel.
//
  c_max = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( c_max < count[i+j*n] )
      {
        c_max = count[i+j*n];
      }
    }
  }
//
//  Set the image data.
//
  r = new int[n*n];
  g = new int[n*n];
  b = new int[n*n];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      if ( count[i+j*n] % 2 == 1 )
      {
        r[i+j*n] = 255;
        g[i+j*n] = 255;
        b[i+j*n] = 255;
      }
      else
      {
        c = ( int ) ( 255.0 * sqrt ( sqrt ( sqrt (
          ( ( double ) ( count[i+j*n] ) / ( double ) ( c_max ) ) ) ) ) );
        r[i+j*n] = 3 * c / 5;
        g[i+j*n] = 3 * c / 5;
        b[i+j*n] = c;
      }
    }
  }
//
//  Write an image file.
//
  ierror = ppma_write ( filename, n, n, r, g, b );

  cout << "\n";
  cout << "  ASCII PPM image data stored in \"" << filename << "\".\n";

  delete [] b;
  delete [] count;
  delete [] g;
  delete [] r;
//
//  Terminate.
//
  cout << "\n";
  cout << "MANDELBROT\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

bool ppma_write ( string file_out_name, int xsize, int ysize, int *r, 
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.
// 
//  Example:
//
//    P3
//    # feep.ppm
//    4 4
//    15
//     0  0  0    0  0  0    0  0  0   15  0 15
//     0  0  0    0 15  7    0  0  0    0  0  0
//     0  0  0    0  0  0    0 15  7    0  0  0
//    15  0 15    0  0  0    0  0  0    0  0  0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
// 
//    28 February 2003
// 
//  Author:
// 
//    John Burkardt
//
//  Parameters:
//
//    Input, string FILE_OUT_NAME, the name of the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE, is
//    true, if an error was detected, or
//    false, if the file was written.
//
{
  int *b_index;
  bool error;
  ofstream file_out;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
//
//  Open the output file.
//
  file_out.open ( file_out_name.c_str ( ) );

  if ( !file_out )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  Cannot open the output file \"" << file_out_name << "\".\n";
    return true;
  }
//
//  Compute the maximum.
//
  rgb_max = 0;
  r_index = r;
  g_index = g;
  b_index = b;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      if ( rgb_max < *r_index )
      {
        rgb_max = *r_index;
      }
      r_index = r_index + 1;

      if ( rgb_max < *g_index )
      {
        rgb_max = *g_index;
      }
      g_index = g_index + 1;

      if ( rgb_max < *b_index )
      {
        rgb_max = *b_index;
      }
      b_index = b_index + 1;
    }
  }
//
//  Write the header.
//
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, rgb_max );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_HEADER failed.\n";
    return true;
  }
//
//  Write the data.
//
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

  if ( error )
  {
    cout << "\n";
    cout << "PPMA_WRITE - Fatal error!\n";
    cout << "  PPMA_WRITE_DATA failed.\n";
    return true;
  }
//
//  Close the file.
//
  file_out.close ( );

  return false;
}
//****************************************************************************80

bool ppma_write_data ( ofstream &file_out, int xsize, int ysize, int *r,
  int *g, int *b )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.
//
//    Output, bool PPMA_WRITE_DATA, is
//    true, if an error was detected, or
//    false, if the data was written.
//
{
  int *b_index;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_num;

  r_index = r;
  g_index = g;
  b_index = b;
  rgb_num = 0;

  for ( j = 0; j < ysize; j++ )
  {
    for ( i = 0; i < xsize; i++ )
    {
      file_out << *r_index << " " << *g_index << " " << *b_index;
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        file_out << "\n";
      }
      else
      {
        file_out << " ";
      }
    }
  }
  return false;
}
//****************************************************************************80

bool ppma_write_header ( ofstream &file_out, string file_out_name, int xsize, 
  int ysize, int rgb_max )

//****************************************************************************80
//
//  Purpose:
//
//    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
//    portable pixel map data.
//
//    Input, string FILE_OUT_NAME, the name of the file.
//
//    Input, int XSIZE, YSIZE, the number of rows and columns of data.
//
//    Input, int RGB_MAX, the maximum RGB value.
//
//    Output, bool PPMA_WRITE_HEADER, is
//    true, if an error was detected, or
//    false, if the header was written.
//
{
  file_out << "P3\n";
  file_out << "# " << file_out_name << " created by PPMA_WRITE.C.\n";
  file_out << xsize << "  " << ysize << "\n";
  file_out << rgb_max << "\n";

  return false;
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
