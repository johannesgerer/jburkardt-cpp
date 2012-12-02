# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

int main ( void );
int ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
  int *g, int *b );
int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r,
  int *g, int *b );
int ppma_write_header ( FILE *file_out, char *file_out_name, int xsize,
  int ysize, int rgb_max );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

    MAIN is the main program for MANDELBROT.

  Discussion:

    MANDELBROT computes an image of the Mandelbrot set.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    22 July 2010

  Author:

    John Burkardt

  Local Parameters:

    Local, integer COUNT_MAX, the maximum number of iterations taken
    for a particular pixel.
*/
{
  int *b;
  int c;
  int c_max;
  int *count;
  int count_max = 400;
  char *filename = "mandelbrot.ppm";
  int *g;
  int i;
  int ierror;
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
  printf ( "\n" );
  printf ( "MANDELBROT\n" );
  printf ( "  C version\n" );
  printf ( "\n" );
  printf ( "  Create an ASCII PPM image of the Mandelbrot set.\n" );
  printf ( "\n" );
  printf ( "  For each point C = X + i*Y\n" );
  printf ( "  with X range [%f,%f]\n", x_min, x_max );
  printf ( "  and  Y range [%f,%f]\n", y_min, y_max );
  printf ( "  carry out %d iterations of the map\n", count_max );
  printf ( "  Z(n+1) = Z(n)^2 + C.\n" );
  printf ( "  If the iterates stay bounded (norm less than 2)\n" );
  printf ( "  then C is taken to be a member of the set.\n" );
  printf ( "\n" );
  printf ( "  An ASCII PPM image of the set is created using\n" );
  printf ( "    N = %d pixels in the X direction and\n", n );
  printf ( "    N = %d pixels in the Y direction.\n", n );
/*
  Carry out the iteration for each pixel, determining COUNT.
*/
  count = ( int * ) malloc ( n * n * sizeof ( int ) );

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
/*
  Determine the coloring of each pixel.
*/
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
/*
  Set the image data.
*/
  r = ( int * ) malloc ( n * n * sizeof ( int ) );
  g = ( int * ) malloc ( n * n * sizeof ( int ) );
  b = ( int * ) malloc ( n * n * sizeof ( int ) );

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
/*
  Write an image file.
*/
  ierror = ppma_write ( filename, n, n, r, g, b );

  printf ( "\n" );
  printf ( "  ASCII PPM image data stored in \"%s\".\n", filename );

  free ( b );
  free ( count );
  free ( g );
  free ( r );
/*
  Terminate.
*/
  printf ( "\n" );
  printf ( "MANDELBROT\n" );
  printf ( "  Normal end of execution.\n" );
  printf ( "\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

int ppma_write ( char *file_out_name, int xsize, int ysize, int *r,
  int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE writes the header and data for an ASCII portable pixel map file.

  Example:

    P3
    # feep.ppm
    4 4
    15
     0  0  0    0  0  0    0  0  0   15  0 15
     0  0  0    0 15  7    0  0  0    0  0  0
     0  0  0    0  0  0    0 15  7    0  0  0
    15  0 15    0  0  0    0  0  0    0  0  0

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, char *FILE_OUT_NAME, the name of the file to contain the ASCII
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.

    Output, int PPMA_WRITE, is
    true, if an error was detected, or
    false, if the file was written.
*/
{
  int *b_index;
  int error;
  FILE *file_out;
  int *g_index;
  int i;
  int j;
  int *r_index;
  int rgb_max;
/*
  Open the output file.
*/
  file_out = fopen ( file_out_name, "wt" );

  if ( !file_out )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  Cannot open the output file \"%s\".\n", file_out_name );
    return 1;
  }
/*
  Compute the maximum.
*/
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
/*
  Write the header.
*/
  error = ppma_write_header ( file_out, file_out_name, xsize, ysize, rgb_max );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_HEADER failed.\n" );
    return 1;
  }
/*
  Write the data.
*/
  error = ppma_write_data ( file_out, xsize, ysize, r, g, b );

  if ( error )
  {
    printf ( "\n" );
    printf ( "PPMA_WRITE - Fatal error!\n" );
    printf ( "  PPMA_WRITE_DATA failed.\n" );
    return 1;
  }
/*
  Close the file.
*/
  fclose ( file_out );

  return 0;
}
/******************************************************************************/

int ppma_write_data ( FILE *file_out, int xsize, int ysize, int *r,
  int *g, int *b )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE_DATA writes the data for an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, ofstream &FILE_OUT, a pointer to the file to contain the ASCII
    portable pixel map data.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int *R, *G, *B, the arrays of XSIZE by YSIZE data values.

    Output, int PPMA_WRITE_DATA, is
    true, if an error was detected, or
    false, if the data was written.
*/
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
      fprintf ( file_out, "%d  %d  %d", *r_index, *g_index, *b_index );
      rgb_num = rgb_num + 3;
      r_index = r_index + 1;
      g_index = g_index + 1;
      b_index = b_index + 1;

      if ( rgb_num % 12 == 0 || i == xsize - 1 || rgb_num == 3 * xsize * ysize )
      {
        fprintf ( file_out, "\n" );
      }
      else
      {
        fprintf ( file_out, " " );
      }
    }
  }
  return 0;
}
/******************************************************************************/

int ppma_write_header ( FILE *file_out, char *file_out_name, int xsize,
  int ysize, int rgb_max )

/******************************************************************************/
/*
  Purpose:

    PPMA_WRITE_HEADER writes the header of an ASCII portable pixel map file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    28 February 2003

  Author:

    John Burkardt

  Parameters:

    Input, FILE *FILE_OUT, a pointer to the file to contain the ASCII
    portable pixel map data.

    Input, char *FILE_OUT_NAME, the name of the file.

    Input, int XSIZE, YSIZE, the number of rows and columns of data.

    Input, int RGB_MAX, the maximum RGB value.

    Output, int PPMA_WRITE_HEADER, is
    true, if an error was detected, or
    false, if the header was written.
*/
{
  fprintf ( file_out, "P3\n" );
  fprintf ( file_out, "# %s created by ppma_write.c.\n", file_out_name );
  fprintf ( file_out, "%d  %d\n", xsize, ysize );
  fprintf ( file_out, "%d\n", rgb_max );

  return 0;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
