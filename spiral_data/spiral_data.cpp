# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iostream>

using namespace std;

# include "spiral_data.hpp"

//****************************************************************************80

void grid_2d ( int x_num, double x_lo, double x_hi, int y_num, double y_lo, 
  double y_hi, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_2D returns a regular 2D grid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int X_NUM, the number of X values to use.
//
//    Input, double X_LO, X_HI, the range of X values.
//
//    Input, int Y_NUM, the number of Y values to use.
//
//    Input, double Y_LO, Y_HI, the range of Y values.
//
//    Output, double X[X_NUM*Y_NUM], Y[X_NUM*Y_NUM], 
//    the coordinates of the grid.
//
{
  int i;
  int j;
  double xi;
  double yj;

  if ( x_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        x[i+j*x_num] = ( x_lo + x_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( i = 0; i < x_num; i++ )
    {
      xi = ( ( double ) ( x_num - i - 1 ) * x_lo   
           + ( double ) (         i     ) * x_hi ) 
           / ( double ) ( x_num     - 1 );
      for ( j = 0; j < y_num; j++ )
      {
        x[i+j*x_num] = xi;
      }
    }
  }

  if ( y_num == 1 )
  {
    for ( j = 0; j < y_num; j++ )
    {
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = ( y_lo + y_hi ) / 2.0;
      }
    }
  }
  else
  {
    for ( j = 0; j < y_num; j++ )
    {
      yj = ( ( double ) ( y_num - j - 1 ) * y_lo   
           + ( double ) (         j     ) * y_hi ) 
           / ( double ) ( y_num     - 1 );
      for ( i = 0; i < x_num; i++ )
      {
        y[i+j*x_num] = yj;
      }
    }
  }

  return;
}
//****************************************************************************80

double r8vec_amax ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMAX returns the maximum absolute value in an R8VEC.
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
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, double AMAX, the value of the entry
//    of largest magnitude.
//
{
  double amax;
  int i;

  amax = 0.0;
  for ( i = 0; i < n; i++ )
  {
    if ( amax < fabs ( a[i] ) )
    {
      amax = fabs ( a[i] );
    }
  }

  return amax;
}
//****************************************************************************80

double r8vec_amin ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_AMIN returns the minimum absolute value in an R8VEC.
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
//    18 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], the array.
//
//    Output, double R8VEC_AMIN, the value of the entry
//    of smallest magnitude.
//
{
  double amin;
  int i;
  const double r8_huge = 1.79769313486231571E+308;

  for ( i = 0; i < n; i++ )
  {
    if ( fabs ( a[i] ) < amin )
    {
      amin = fabs ( a[i] );
    }
  }

  return amin;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
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
//    22 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.  This
//    is set to 0.0 if N <= 0.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
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
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double value;

  value = r8vec[0];

  for ( i = 1; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
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
  const int i4_huge = 2147483647;
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

void resid_spiral ( int n, double x[], double y[], double c, double pr[] )

//****************************************************************************80
//
//  Purpose:
//
//    RESID_SPIRAL computes the residual for a spiral velocity vector field.
//
//  Discussion:
//
//    Note that the continuous velocity field (U,V)(X,Y) that is discretely
//    sampled here satisfies the homogeneous continuity equation, that is,
//    it has zero divergence.  In other words:
//
//      dU/dX + dV/dY = 0.
//
//    This is by construction, since we have
//
//      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
//      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )
//
//    which guarantees zero divergence.
//
//    The underlying function PHI is defined by
//
//      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2
//
//    where C is a parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the 
//    evaluation points.
//
//    Input, double C, a parameter, typically between 0 and 2 * PI.
//
//    Output, double PR[N], the residual in the continuity equation.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double u;
  double ux;
  double v;
  double vy;

  for ( i = 0; i < n; i++ )
  {
    u =   10.0 * ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
           * pow ( 1.0 - x[i], 2 )
           * ( 
               c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
             - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
             * 2.0 * ( 1.0 - y[i] ) 
             );

    ux =   10.0 * 
      ( 
        c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
      ) 
      * 
      ( 
        c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
        * 2.0 * ( 1.0 - y[i] ) 
      );

    v = - 10.0 * ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
      * pow ( 1.0 - y[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
        );

    vy =  - 10.0 * 
      ( 
        c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
      ) 
      * 
      ( 
        c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
        * 2.0 * ( 1.0 - y[i] ) 
      );

    pr[i] = ux + vy;
  }

  return;
}
//****************************************************************************80

void spiral_gnuplot ( string header, int n, double x[], double y[], double u[], 
  double v[], double s )

//****************************************************************************80
//
//  Purpose:
//
//    SPIRAL_GNUPLOT writes the spiral vector field to files for GNUPLOT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string HEADER, a header to be used to name the files.
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the evaluation points.
//
//    Input, double U[N], V[N], the velocity components.
//
//    Input, double S, a scale factor for the velocity vectors.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  string plot_filename;
//
//  Write the data file.
//
  data_filename = header + "_data.txt";

  data_unit.open ( data_filename.c_str ( ) );

  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << x[i]
              << "  " << y[i]
              << "  " << s * u[i]
              << "  " << s * v[i] << "\n";
  }

  data_unit.close ( );

  cout << "\n";
  cout << "  Data written to '" << data_filename << "'\n";
//
//  Write the command file.
//
  command_filename = header + "_commands.txt";
  plot_filename = header + ".png";

  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "#  " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "#\n";
  command_unit << "#  Add titles and labels.\n";
  command_unit << "#\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title 'Spiral velocity flow'\n";
  command_unit << "unset key\n";
  command_unit << "#\n";
  command_unit << "#  Add grid lines.\n";
  command_unit << "#\n";
  command_unit << "set grid\n";
  command_unit << "set size ratio -1\n";
  command_unit << "#\n";
  command_unit << "#  Timestamp the plot.\n";
  command_unit << "#\n";
  command_unit << "set timestamp\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2:3:4 with vectors \\\n";
  command_unit << "  head filled lt 2 linecolor rgb 'blue'\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Commands written to '" << command_filename << "'\n";

  return;
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

void uv_spiral ( int n, double x[], double y[], double c, double u[], 
  double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    UV_SPIRAL computes a spiral velocity vector field.
//
//  Discussion:
//
//    Note that the continuous velocity field (U,V)(X,Y) that is discretely
//    sampled here satisfies the homogeneous continuity equation, that is,
//    it has zero divergence.  In other words:
//
//      dU/dX + dV/dY = 0.
//
//    This is by construction, since we have
//
//      U(X,Y) =  10 * d/dY ( PHI(X) * PHI(Y) )
//      V(X,Y) = -10 * d/dX ( PHI(X) * PHI(Y) )
//
//    which guarantees zero divergence.
//
//    The underlying function PHI is defined by
//
//      PHI(Z) = ( 1 - cos ( C * pi * Z ) ) * ( 1 - Z )^2
//
//    where C is a parameter.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 January 2015
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of evaluation points.
//
//    Input, double X[N], Y[N], the coordinates of the 
//    evaluation points.
//
//    Input, double C, a parameter, typically between 0 and 2 * PI.
//
//    Output, double U[N], V[N], the velocity components.
//
{
  int i;
  const double r8_pi = 3.141592653589793;

  for ( i = 0; i < n; i++ )
  {
    u[i] =   10.0 * ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
      * pow ( 1.0 - x[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * y[i] ) * pow ( 1.0 - y[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
          * 2.0 * ( 1.0 - y[i] ) 
        );

    v[i] = - 10.0 * ( 1.0 - cos ( c * r8_pi * y[i] ) ) 
      * pow ( 1.0 - y[i], 2 )
      * ( 
          c * r8_pi * sin ( c * r8_pi * x[i] ) * pow ( 1.0 - x[i], 2 )
        - ( 1.0 - cos ( c * r8_pi * x[i] ) ) 
        * 2.0 * ( 1.0 - x[i] ) 
        );

  }

  return;
}
