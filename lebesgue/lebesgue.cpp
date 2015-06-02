# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "lebesgue.hpp"

//****************************************************************************80

double *chebyshev1 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV1 returns the Type 1 Chebyshev points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double CHEBYSHEV1[N], the points.
//
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * i + 1 ) / ( double ) ( 2 * n );
    x[i] = cos ( angle );
  }
  return x;
}
//****************************************************************************80

double *chebyshev2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV2 returns the Type 2 Chebyshev points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double CHEBYSHEV2[N], the points.
//
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      angle = r8_pi * ( double ) ( n - i - 1 ) / ( double ) ( n - 1 );
      x[i] = cos ( angle );
    }
  }

  return x;
}
//****************************************************************************80

double *chebyshev3 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV3 returns the Type 3 Chebyshev points.
//
//  Discussion:
//
//    Note that this point set is NOT symmetric in [-1,+1].
//    It is sometimes augmented by the value -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double CHEBYSHEV3[N], the points.
//
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * n - 2 * i - 1 ) 
                  / ( double ) ( 2 * n         + 1 );
    x[i] = cos ( angle );
  }

  return x;
}
//****************************************************************************80

double *chebyshev4 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBYSHEV4 returns the Type 4 Chebyshev points.
//
//  Discussion:
//
//    Note that this point set is NOT symmetric in [-1,+1].
//    It is sometimes augmented by the value +1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double CHEBYSHEV4[N], the points.
//
{
  double angle;
  int i;
  const double r8_pi = 3.141592653589793;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    angle = r8_pi * ( double ) ( 2 * n - 2 * i )
                  / ( double ) ( 2 * n + 1 );
    x[i] = cos ( angle );
  }

  return x;
}
//****************************************************************************80

double *equidistant1 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIDISTANT1 returns the Type 1 Equidistant points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double EQUIDISTANT1[N], the points.
//
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n + 1 );
  }

  return x;
}
//****************************************************************************80

double *equidistant2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIDISTANT2 returns the Type 2 Equidistant points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double EQUIDISTANT2[N], the points.
//
{
  int i;
  double *x;

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n - 1 );
    }
  }

  return x;
}
//****************************************************************************80

double *equidistant3 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    EQUIDISTANT3 returns the Type 3 Equidistant points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double EQUIDISTANT3[N], the points.
//
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( double ) ( - n + 1 + 2 * i ) / ( double ) ( n );
  }

  return x;
}
//****************************************************************************80

double *fejer1 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER1 returns the Type 1 Fejer points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double FEJER1[N], the points.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( 2 * n - 1 - 2 * i ) 
                  / ( double ) ( 2 * n );
    x[i] = cos ( theta );
  }
  return x;
}
//****************************************************************************80

double *fejer2 ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    FEJER2 returns the Type 2 Fejer points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2018
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double FEJER2[N], the points.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( n - i ) 
                  / ( double ) ( n + 1 );
    x[i] = cos ( theta );
  }

  return x;
}
//****************************************************************************80

double *lagrange_value ( int data_num, double t_data[], int interp_num, 
  double t_interp[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE evaluates the Lagrange polynomials.
//
//  Discussion:
//
//    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
//    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
//    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
//    other abscissas.
//
//    A formal representation is:
//
//      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
//       ( T - T(J) ) / ( T(I) - T(J) )
//
//    This routine accepts a set of INTERP_NUM values at which all the Lagrange
//    polynomials should be evaluated.
//
//    Given data values P_DATA at each of the abscissas, the value of the
//    Lagrange interpolating polynomial at each of the interpolation points
//    is then simple to compute by matrix multiplication:
//
//      P_INTERP(1:INTERP_NUM) =
//        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
//
//    or, in the case where P is multidimensional:
//
//      P_INTERP(1:M,1:INTERP_NUM) =
//        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 December 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the number of data points.
//    DATA_NUM must be at least 1.
//
//    Input, double T_DATA[DATA_NUM], the data points.
//
//    Input, int INTERP_NUM, the number of
//    interpolation points.
//
//    Input, double T_INTERP[INTERP_NUM], the
//    interpolation points.
//
//    Output, double LAGRANGE_VALUE[DATA_NUM*INTERP_NUM], the values
//    of the Lagrange polynomials at the interpolation points.
//
{
  int i;
  int i1;
  int i2;
  int j;
  double *l_interp;

  l_interp = new double[data_num*interp_num];
//
//  Evaluate the polynomial.
//
  for ( j = 0; j < interp_num; j++ )
  {
    for ( i = 0; i < data_num; i++ )
    {
      l_interp[i+j*data_num] = 1.0;
    }
  }

  for ( i1 = 0; i1 < data_num; i1++ )
  {
    for ( i2 = 0; i2 < data_num; i2++ )
    {
      if ( i1 != i2 )
      {
        for ( j = 0; j < interp_num; j++ )
        {
          l_interp[i1+j*data_num] = l_interp[i1+j*data_num] 
            * ( t_interp[j] - t_data[i2] ) / ( t_data[i1] - t_data[i2] );
        }
      }
    }
  }

  return l_interp;
}
//****************************************************************************80

double lebesgue_constant ( int n, double x[], int nfun, double xfun[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_CONSTANT estimates the Lebesgue constant for a set of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
//
//  Parameters:
//
//    Input, int N, the number of interpolation points.
//
//    Input, double X[N], the interpolation points.
//
//    Input, int NFUN, the number of evaluation points.
//
//    Input, double XFUN[NFUN], the evaluation points.
//
//    Output, double LEBESGUE_CONSTANT, an estimate of the Lebesgue constant 
//    for the points.
//
{
  double *lfun;
  double lmax;

  lfun = lebesgue_function ( n, x, nfun, xfun );

  lmax = r8vec_max ( nfun, lfun );

  delete [] lfun;

  return lmax;
}
//****************************************************************************80

double *lebesgue_function ( int n, double x[], int nfun, double xfun[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_FUNCTION evaluates the Lebesgue function for a set of points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
//
//  Parameters:
//
//    Input, int N, the number of interpolation points.
//
//    Input, double X[N], the interpolation points.
//
//    Input, int NFUN, the number of evaluation points.
//
//    Input, double XFUN[NFUN], the evaluation points.
//
//    Output, double LEBESGUE_FUNCTION[NFUN], the Lebesgue function values.
//
{
  int i;
  int j;
  double *lfun;
  double *llfun;
  double t;

  lfun = new double[nfun];
//
//  Handle special case.
//
  if ( n == 1 )
  {
    for ( j = 0; j < nfun; j++ )
    {
      lfun[j] = 1.0;
    }
    return lfun;
  }

  llfun = lagrange_value ( n, x, nfun, xfun );

  for ( j = 0; j < nfun; j++ )
  {
    t = 0.0;
    for ( i = 0; i < n; i++ )
    {
      t = t + fabs ( llfun[i+j*n] );
    }
    lfun[j] = t;
  }

  delete [] llfun;

  return lfun;
}
//****************************************************************************80

void lebesgue_plot ( int n, double x[], int nfun, double xfun[], 
  string label, string filename )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_PLOT plots the Lebesgue function for a set of points.
//
//  Discussion:
//
//    The interpolation interval is assumed to be [min(XFUN), max(XFUN)].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Jean-Paul Berrut, Lloyd Trefethen,
//    Barycentric Lagrange Interpolation,
//    SIAM Review,
//    Volume 46, Number 3, September 2004, pages 501-517.
//
//  Parameters:
//
//    Input, int N, the number of interpolation points.
//
//    Input, double X[N], the interpolation points.
//
//    Input, int NFUN, the number of evaluation points.
//
//    Input, double XFUN[NFUN], the evaluation points.  
//
//    Input, string LABEL, a title for the plot.
//
//    Input, string FILENAME, a partial filename.
//    The program will create "filename_commands.txt', 'filename_data.txt',
//    and 'filename.png'.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  double *lfun;
  string png_filename;

  lfun = lebesgue_function ( n, x, nfun, xfun );
//
//  Create data file.
//
  data_filename = filename + "_data.txt";
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < nfun; i++ )
  {
    data_unit << xfun[i] << "  " 
              << lfun[i] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
//
//  Create command file.
//
  command_filename = filename + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";

  png_filename = filename + ".png";
  command_unit << "set output '" << png_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Lebesgue(X) --->'\n";
  command_unit << "set title '" << label << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "plot '" << data_filename 
               << "' using 1:2 lw 3 linecolor rgb 'red'\n";

  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] lfun;

  return;
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

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << ": " << setw(14) << a[i]  << "\n";
  }

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
