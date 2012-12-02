# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <cstring>
# include <ctime>

using namespace std;

# include "correlation.hpp"

//****************************************************************************80

double *correlation_besselj ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_BESSELJ evaluates the Bessel J correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_abs ( rho[i] ) / rho0;
    c[i] = r8_besj0 ( rhohat );
  }

  return c;
}
//****************************************************************************80

double *correlation_besselk ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_BESSELK evaluates the Bessel K correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( rho[i] == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      rhohat = r8_abs ( rho[i] ) / rho0;
      c[i] = rhohat * r8_besk1 ( rhohat );
    }
  }

  return c;
}
//****************************************************************************80

double *correlation_brownian ( int m, int n, double s[], double t[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_BROWNIAN computes the Brownian correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of arguments.
//
//    Input, double S[M], T[N], two samples.
//    0 <= S(*), T(*).
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[M*N], the correlations.
//
{
  double *c;
  int i;
  int j;

  c = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( 0.0 < r8_max ( s[i], t[j] ) )
      {
        c[i+j*m] = sqrt ( r8_min ( s[i], t[j] ) / r8_max ( s[i], t[j] ) );
      }
      else
      {
        c[i+j*m] = 1.0;
      }
    }
  }

  return c;
}
//****************************************************************************80

void correlation_brownian_display ( )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_BROWNIAN_DISPLAY displays 4 slices of the Brownian Correlation.
//
//  Discussion:
//
//    The correlation function is C(S,T) = sqrt ( min ( s, t ) / max ( s, t ) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
{
  double *c;
  string command_filename = "brownian_plots_commands.txt";
  ofstream command_unit;
  string data_filename = "brownian_plots_data.txt";
  ofstream data_unit;
  int i;
  int j;
  int n = 101;
  int n2 = 4;
  double *s;
  double t[4] = { 0.25, 1.50, 2.50, 3.75 };

  s = r8vec_linspace_new ( n, 0.0, 5.0 );

  c = new double[n*n2];

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < n2; j++ )
    {
      c[i+j*n] = sqrt ( r8_min ( s[i], t[j] ) / r8_max ( s[i], t[j] ) );
    }
  }

  data_unit.open ( data_filename.c_str() );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << s[i];
    for ( j = 0; j < n2; j++ )
    {
      data_unit << "  " << c[i+j*n];
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command_unit.open ( command_filename.c_str() );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set key off\n";
  command_unit << "set output \"brownian_plots.png\"\n";
  command_unit << "set title 'Brownian correlation C(S,T), S = 0.25, 1.5, 2.5, 3.75'\n";
  command_unit << "set xlabel 'S'\n";
  command_unit << "set ylabel 'C(s,t)'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot \"" << data_filename << "\" using 1:2 lw 3 linecolor rgb 'blue',\\\n";
  command_unit << "     \"" << data_filename << "\" using 1:3 lw 3 linecolor rgb 'blue',\\\n";
  command_unit << "     \"" << data_filename << "\" using 1:4 lw 3 linecolor rgb 'blue',\\\n";
  command_unit << "     \"" << data_filename << "\" using 1:5 lw 3 linecolor rgb 'blue'\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file \"" << command_filename << "\".\n";

  delete [] c;
  delete [] s;

  return;
}
//****************************************************************************80

double *correlation_circular ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_CIRCULAR evaluates the circular correlation function.
//
//  Discussion:
//
//    This correlation is based on the area of overlap of two circles
//    of radius RHO0 and separation RHO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double pi = 3.141592653589793;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_min ( r8_abs ( rho[i] ) / rho0, 1.0 );

    c[i] = ( 1.0 - ( 2.0 / pi ) 
      * ( rhohat * sqrt ( 1.0 - rhohat * rhohat ) + asin ( rhohat ) ) );
  }

  return c;
}
//****************************************************************************80

double *correlation_constant ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_CONSTANT evaluates the constant correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = 1.0;
  }
  return c;
}
//****************************************************************************80

double *correlation_cubic ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_CUBIC evaluates the cubic correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_min ( r8_abs ( rho[i] ) / rho0, 1.0 );

    c[i] = 1.0 
         - 7.0  * pow ( rhohat, 2 ) 
         + 8.75 * pow ( rhohat, 3 )
         - 3.5  * pow ( rhohat, 5 )
         + 0.75 * pow ( rhohat, 7 );
  }

  return c;
}
//****************************************************************************80

double *correlation_damped_cosine ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_DAMPED_COSINE evaluates the damped cosine correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = exp ( - r8_abs ( rho[i] ) / rho0 ) * cos ( r8_abs ( rho[i] ) / rho0 );
  }
  return c;
}
//****************************************************************************80

double *correlation_damped_sine ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( rho[i] == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      rhohat = r8_abs ( rho[i] ) / rho0;
      c[i] = sin ( rhohat ) / rhohat;
    }
  }
  return c;
}
//****************************************************************************80

double *correlation_exponential ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_EXPONENTIAL evaluates the exponential correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];
  for ( i = 0; i < n; i++ )
  {
    c[i] = exp ( - r8_abs ( rho[i] ) / rho0 );
  }
  return c;
}
//****************************************************************************80

double *correlation_gaussian ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_GAUSSIAN evaluates the Gaussian correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = exp ( - pow ( rho[i] / rho0, 2 ) );
  }
  return c;
}
//****************************************************************************80

double *correlation_hole ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_HOLE evaluates the hole correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = ( 1.0 - r8_abs ( rho[i] ) / rho0 ) 
      * exp ( - r8_abs ( rho[i] ) / rho0 );
  }
  return c;
}
//****************************************************************************80

double *correlation_linear ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_LINEAR evaluates the linear correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( rho0 < r8_abs ( rho[i] ) )
    {
      c[i] = 0.0;
    }
    else
    {
      c[i] = ( rho0 - r8_abs ( rho[i] ) ) / rho0;
    }
  }
  return c;
}
//****************************************************************************80

double *correlation_matern ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_MATERN evaluates the Matern correlation function.
//
//  Discussion:
//
//    In order to call this routine under a dummy name, I had to drop NU from
//    the parameter list.
//
//    The Matern correlation is
//
//      rho1 = 2 * sqrt ( nu ) * rho / rho0
//
//      c(rho) = ( rho1 )^nu * BesselK ( nu, rho1 ) 
//               / gamma ( nu ) / 2 ^ ( nu - 1 )
//
//    The Matern covariance has the form:
//
//      K(rho) = sigma^2 * c(rho)
//
//    A Gaussian process with Matern covariance has sample paths that are
//    differentiable (nu - 1) times.
//
//    When nu = 0.5, the Matern covariance is the exponential covariance.
//
//    As nu goes to +oo, the correlation converges to exp ( - (rho/rho0)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//    0.0 <= RHO.
//
//    Input, double RHO0, the correlation length.
//    0.0 < RHO0.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double nu;
  double rho1;

  nu = 2.5;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rho1 = 2.0 * sqrt ( nu ) * r8_abs ( rho[i] ) / rho0;

    if ( rho1 == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      c[i] = pow ( rho1, nu ) * r8_besk ( nu, rho1 ) / r8_gamma ( nu ) 
        / pow ( 2.0, nu - 1.0 );
    }
  }
  return c;
}
//****************************************************************************80

double *correlation_pentaspherical ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_PENTASPHERICAL evaluates the pentaspherical correlation function.
//
//  Discussion:
//
//    This correlation is based on the volume of overlap of two spheres
//    of radius RHO0 and separation RHO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_min ( r8_abs ( rho[i] ) / rho0, 1.0 );

    c[i] = 1.0 - 1.875 * rhohat + 1.25 * pow ( rhohat, 3 )
      - 0.375 * pow ( rhohat, 5 );
  }

  return c;
}
//****************************************************************************80

void correlation_plot ( int n, double rho[], double c[], string header, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_PLOT makes a plot of a correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double C[N], the correlations.
//
//    Input, string HEADER, an identifier for the files.
//
//    Input, string TITLE, a title for the plot.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  double rho0;

  data_filename = header + "_data.txt";

  data_unit.open ( data_filename.c_str() );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << rho[i]
              << "  " << c[i] << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";
  command_filename = header + "_commands.txt";

  command_unit.open ( command_filename.c_str() );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output \"" << header << "_plot.png\"\n";
  command_unit << "set xlabel 'Distance Rho'\n";
  command_unit << "set ylabel 'Correlation C(Rho)'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'blue'\n";
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void correlation_plots ( int n, int n2, double rho[], double rho0[], double c[], 
  string header, string title )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_PLOTS plots correlations for a range of correlation lengths.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values of RHO.
//
//    Input, int N2, the number of values of RHO0.
//
//    Input, double RHO[N], the independent value.
//
//    Input, double RHO0[N2], the correlation lengths.
//
//    Input, double C[N*N2], the correlations.
//
//    Input, string HEADER, an identifier for the files.
//
//    Input, string TITLE, a title for the plot.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;

  data_filename = header + "_plots_data.txt";

  data_unit.open ( data_filename.c_str() );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << rho[i];
    for ( j = 0; j < n2; j++ )
    {
      data_unit << "  " << c[i+j*n];
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command_filename = header + "_plots_commands.txt";

  command_unit.open ( command_filename.c_str() );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output \"" << header << "_plots.png\"\n";
  command_unit << "set xlabel 'Rho'\n";
  command_unit << "set ylabel 'Correlation(Rho)'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set key off\n";
  if ( n2 == 1 )
  {
    command_unit << "plot '" << data_filename << "' using 1:2 lw 3\n";
  }
  else
  {
    command_unit << "plot '" << data_filename << "' using 1:2 lw 3, \\\n";
    for ( i = 2; i < n2; i++ )
    {
      command_unit << "     '" << data_filename << "' using 1:" << i + 1 << " lw 3, \\\n";
    }
    command_unit << "     '" << data_filename << "' using 1:" << n2 + 1 << " lw 3\n";
  }
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

double *correlation_power ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_POWER evaluates the power correlation function.
//
//  Discussion:
//
//    In order to be able to call this routine under a dummy name, I had
//    to drop E from the argument list.
//
//    The power correlation is
//
//      C(rho) = ( 1 - |rho| )^e  if 0 <= |rho| <= 1
//             = 0                otherwise
//
//      The constraint on the exponent is 2 <= e.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//    0.0 <= RHO.
//
//    Input, double RHO0, the correlation length.
//    0.0 < RHO0.
//
//    Input, double E, the exponent.
//    E has a default value of 2.0;
//    2.0 <= E.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  double e;
  int i;
  double rhohat;

  e = 2.0;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_abs ( rho[i] ) / rho0;
    if ( rhohat <= 1.0 )
    {
      c[i] = pow ( 1.0 - rhohat, e );
    }
    else
    {
      c[i] = 0.0;
    }
  }
  return c;
}
//****************************************************************************80

double *correlation_rational_quadratic ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_RATIONAL_QUADRATIC: rational quadratic correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    c[i] = 1.0 / ( 1.0 + pow ( rho[i] / rho0, 2 ) );
  }

  return c;
}
//****************************************************************************80

double *correlation_spherical ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_SPHERICAL evaluates the spherical correlation function.
//
//  Discussion:
//
//    This correlation is based on the volume of overlap of two spheres
//    of radius RHO0 and separation RHO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;
  double rhohat;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    rhohat = r8_min ( r8_abs ( rho[i] ) / rho0, 1.0 );
    c[i] = 1.0 - 1.5 * rhohat + 0.5 * pow ( rhohat, 3 );
  }

  return c;
}
//****************************************************************************80

double *correlation_to_covariance ( int n, double c[], double sigma[] )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_TO_COVARIANCE: covariance matrix from a correlation matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double C[N*N], the correlation matrix.
//
//    Input, double SIGMA[N], the standard deviations.
//
//    Output, double K[N*N], the covariance matrix.
//
{
  double c_max;
  double c_min;
  double e;
  double error_frobenius;
  int i;
  int j;
  double *k;
  double tol;

  tol = sqrt ( r8_epsilon ( ) );
//
//  C must be symmetric.
//
  error_frobenius = r8mat_is_symmetric ( n, n, c );

  if ( tol < error_frobenius )
  {
    cerr << "\n";
    cerr << "CORRELATION_TO_COVARIANCE - Fatal error!\n";
    cerr << "  Input matrix C fails symmetry test with error " << error_frobenius << "\n";
    exit ( 1 );
  }
//
//  The diagonal must be 1.
//
  for ( i = 0; i < n; i++ )
  {
    e = r8_abs ( c[i+i*n] - 1.0 );
    if ( tol < e )
    {
      cerr << "\n";
      cerr << "CORRELATION_TO_COVARIANCE - Fatal error!\n";
      cerr << "  Input matrix C has non-unit diagonal entries.\n";
      cerr << "  Error on row " << i << " is " << e << "\n";
      exit ( 1 );
    }
  }
//
//  Off-diagonals must be between -1 and 1.
//
  c_min = r8mat_min ( n, n, c );

  if ( c_min < - 1.0 - tol )
  {
    cerr << "\n";
    cerr << "CORRELATION_TO_COVARIANCE - Fatal error!\n";
    cerr << "  Input matrix C has entries less than -1.0\n";
    exit ( 1 );
  }

  c_max = r8mat_max ( n, n, c );

  if ( 1.0 + tol < c_max )
  {
    cerr << "\n";
    cerr << "CORRELATION_TO_COVARIANCE - Fatal error!\n";
    cerr << "  Input matrix C has entries greater than +1.0\n";
    exit ( 1 );
  }
//
//  Form K.
//
  k = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      k[i+j*n] = sigma[i] * c[i+j*n] * sigma[j];
    }
  }

  return k;
}
//****************************************************************************80

double *correlation_white_noise ( int n, double rho[], double rho0 )

//****************************************************************************80
//
//  Purpose:
//
//    CORRELATION_WHITE_NOISE evaluates the white noise correlation function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Petter Abrahamsen,
//    A Review of Gaussian Random Fields and Correlation Functions,
//    Norwegian Computing Center, 1997.
//
//  Parameters:
//
//    Input, int N, the number of arguments.
//
//    Input, double RHO[N], the arguments.
//
//    Input, double RHO0, the correlation length.
//
//    Output, double C[N], the correlations.
//
{
  double *c;
  int i;

  c = new double[n];

  for ( i = 0; i < n; i++ )
  {
    if ( rho[i] == 0.0 )
    {
      c[i] = 1.0;
    }
    else
    {
      c[i] = 0.0;
    }
  }

  return c;
}
//****************************************************************************80

void covariance_to_correlation ( int n, double k[], double c[], double sigma[] )

//****************************************************************************80
//
//  Purpose:
//
//    COVARIANCE_TO_CORRELATION: correlation matrix from a covariance matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double K[N*N], the covariance matrix.
//
//    Output, double C[N*N], the correlation matrix.
//
//    Output, double SIGMA[N], the standard deviations.
//
{
  double e;
  double error_frobenius;
  int i;
  int j;
  double sigma_min;
  double tol;

  tol = sqrt ( r8_epsilon ( ) );
//
//  K must be symmetric.
//
  error_frobenius = r8mat_is_symmetric ( n, n, k );

  if ( tol < error_frobenius )
  {
    cerr << "\n";
    cerr << "COVARIANCE_TO_CORRELATION - Fatal error\n";
    cerr << "  Input matrix K fails symmetry test with error " << error_frobenius << "\n";
    exit ( 1 );
  }
//
//  It must be the case that K(I,J)^2 <= K(I,I) * K(J,J).
//
  e = 0.0;
  for ( i = 0; i < n; i++ )
  {
    for ( j = i + 1; j < n; j++ )
    {
      e = r8_max ( e, k[i+j*n] * k[i+j*n] - k[i+i*n] * k[j+j*n] );
    }
  }

  if ( tol < e )
  {
    cerr << "\n";
    cerr << "COVARIANCE_TO_CORRELATION - Fatal error\n";
    cerr << "  Input matrix K fails K(I,J)^2 <= K(I,I)*K(J,J)\n";
    exit ( 1 );
  }
//
//  Get the diagonal.
//
  for ( i = 0; i < n; i++ )
  {
    sigma[i] = k[i+i*n];
  }
//
//  Ensure the diagonal is positive.
//
  sigma_min = r8vec_min ( n, sigma );

  if ( sigma_min <= 0.0 )
  {
    cerr << "\n";
    cerr << "COVARIANCE_TO_CORRELATION - Fatal error!\n";
    cerr << "  Input matrix K has nonpositive diagonal entry = " << sigma_min << "\n";
    exit ( 1 );
  }
//
//  Convert from variance to standard deviation.
//
  for ( i = 0; i < n; i++ )
  {
    sigma[i] = sqrt ( sigma[i] );
  }
//
//  Form C.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = k[i+j*n] / sigma[i] / sigma[j];
    }
  }

  return;
}
//****************************************************************************80

int i4_abs ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_ABS returns the absolute value of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, an integer.
//
//    Output, int I4_ABS, the absolute value of the integer.
//
{
  int value;

  if ( 0 <= i )
  {
    value = i;
  }
  else
  {
    value = - i;
  }
  return value;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

double *minij ( int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    MINIJ returns the MINIJ matrix.
//
//  Discussion:
//
//    A(I,J) = min ( I, J )
//
//  Example:
//
//    N = 5
//
//    1 1 1 1 1
//    1 2 2 2 2
//    1 2 3 3 3
//    1 2 3 4 4
//    1 2 3 4 5
//
//  Properties:
//
//    A is integral, therefore det ( A ) is integral, and 
//    det ( A ) * inverse ( A ) is integral.
//
//    A is positive definite.
//
//    A is symmetric: A' = A.
//
//    Because A is symmetric, it is normal.
//
//    Because A is normal, it is diagonalizable.
//
//    The inverse of A is tridiagonal.
//
//    The eigenvalues of A are
//
//      LAMBDA[i] = 0.5 / ( 1 - cos ( ( 2 * I - 1 ) * pi / ( 2 * N + 1 ) ) ),
//
//    (N+1)*ONES[N] - A also has a tridiagonal inverse.
//
//    Gregory and Karney consider the matrix defined by
//
//      B(I,J) = N + 1 - MAX(I,J)
//
//    which is equal to the MINIJ matrix, but with the rows and
//    columns reversed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Robert Gregory, David Karney,
//    Example 3.12, Example 4.14,
//    A Collection of Matrices for Testing Computational Algorithms,
//    Wiley, 1969, page 41, page 74, 
//    LC: QA263.G68.
//
//    Daniel Rutherford,
//    Some continuant determinants arising in physics and chemistry II,
//    Proceedings of the Royal Society Edinburgh,
//    Volume 63, A, 1952, pages 232-241.
//
//    John Todd,
//    Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
//    Academic Press, 1977, page 158.
//
//    Joan Westlake,
//    A Handbook of Numerical Matrix Inversion and Solution of 
//    Linear Equations,
//    John Wiley, 1968,
//    ISBN13: 978-0471936756,
//    LC: QA263.W47.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns 
//    of the matrix.
//
//    Output, double MINIJ[M*N], the matrix.
//
{
  double *a;
  int i;
  int j;

  a = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = ( double ) ( i4_min ( i + 1, j + 1 ) );
    }
  }

  return a;
}
//****************************************************************************80

void paths_plot ( int n, int n2, double rho[], double x[], string header, 
  string title )

//****************************************************************************80
//
//  Purpose:
//
//    PATHS_PLOT plots a sequence of paths or simulations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points in each path.
//
//    Input, int N2, the number of paths.
//
//    Input, double RHO[N], the independent value.
//
//    Input, double X[N*N2], the path values.
//
//    Input, string HEADER, an identifier for the files.
//
//    Input, string TITLE, a title for the plot.
//
{
  string command_filename;
  ofstream command_unit;
  string data_filename;
  ofstream data_unit;
  int i;
  int j;
  double rho0;

  data_filename = header + "_path_data.txt";

  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < n; i++ )
  {
    data_unit << "  " << rho[i];
    for ( j = 0; j < n2; j++ )
    {
      data_unit << "  " << x[i+j*n];
    }
    data_unit << "\n";
  }
  data_unit.close ( );
  cout << "  Created data file \"" << data_filename << "\".\n";

  command_filename = header + "_path_commands.txt";

  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output \"" << header << "_paths.png\"\n";
  command_unit << "set xlabel 'Rho'\n";
  command_unit << "set ylabel 'X(Rho)'\n";
  command_unit << "set title '" << title << "'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "set key off\n";
  if ( n2 == 1 )
  {
    command_unit << "plot '" << data_filename << "' using 1:2 lw 3\n";
  }
  else
  {
    command_unit << "plot '" << data_filename << "' using 1:2, \\\n";
    for ( i = 2; i < n2; i++ )
    {
      command_unit << "     '" << data_filename << "' using 1:" << i + 1 << ", \\\n";
    }
    command_unit << "     '" << data_filename << "' using 1:" << n2 + 1 << "\n";
  }
  command_unit << "quit\n";
  command_unit.close ( );
  cout << "  Created command file \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

double pythag ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    PYTHAG computes SQRT ( A * A + B * B ) carefully.
//
//  Discussion:
//
//    The formula
//
//      PYTHAG = sqrt ( A * A + B * B )
//
//    is reasonably accurate, but can fail if, for example, A^2 is larger
//    than the machine overflow.  The formula can lose most of its accuracy
//    if the sum of the squares is very large or very small.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Modified:
//
//    08 November 2012
//
//  Parameters:
//
//    Input, double A, B, the two legs of a right triangle.
//
//    Output, double PYTHAG, the length of the hypotenuse.
//
{
  double p;
  double r;
  double s;
  double t;
  double u;

  p = r8_max ( r8_abs ( a ), r8_abs ( b ) );

  if ( p != 0.0 )
  {
    r = r8_min ( r8_abs ( a ), r8_abs ( b ) ) / p;
    r = r * r;

    while ( 1 )
    {
      t = 4.0 + r;

      if ( t == 4.0 )
      {
        break;
      }

      s = r / t;
      u = 1.0 + 2.0 * s;
      p = u * p;
      r = ( s / u ) * ( s / u ) * r;
    }
  }
  return p;
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_aint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_AINT truncates an R8 argument to an integer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    1 September 2011
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_AINT, the truncated version of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = - ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }
  else
  {
    value =   ( double ) ( ( int ) ( r8_abs ( x ) ) );
  }

  return value;
}
//****************************************************************************80

void r8_b0mp ( double x, double &ampl, double &theta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double &AMPL, &THETA, the modulus and phase.
//
{
  static double bm0cs[37] = {
    +0.9211656246827742712573767730182E-01,
    -0.1050590997271905102480716371755E-02,
    +0.1470159840768759754056392850952E-04,
    -0.5058557606038554223347929327702E-06,
    +0.2787254538632444176630356137881E-07,
    -0.2062363611780914802618841018973E-08,
    +0.1870214313138879675138172596261E-09,
    -0.1969330971135636200241730777825E-10,
    +0.2325973793999275444012508818052E-11,
    -0.3009520344938250272851224734482E-12,
    +0.4194521333850669181471206768646E-13,
    -0.6219449312188445825973267429564E-14,
    +0.9718260411336068469601765885269E-15,
    -0.1588478585701075207366635966937E-15,
    +0.2700072193671308890086217324458E-16,
    -0.4750092365234008992477504786773E-17,
    +0.8615128162604370873191703746560E-18,
    -0.1605608686956144815745602703359E-18,
    +0.3066513987314482975188539801599E-19,
    -0.5987764223193956430696505617066E-20,
    +0.1192971253748248306489069841066E-20,
    -0.2420969142044805489484682581333E-21,
    +0.4996751760510616453371002879999E-22,
    -0.1047493639351158510095040511999E-22,
    +0.2227786843797468101048183466666E-23,
    -0.4801813239398162862370542933333E-24,
    +0.1047962723470959956476996266666E-24,
    -0.2313858165678615325101260800000E-25,
    +0.5164823088462674211635199999999E-26,
    -0.1164691191850065389525401599999E-26,
    +0.2651788486043319282958336000000E-27,
    -0.6092559503825728497691306666666E-28,
    +0.1411804686144259308038826666666E-28,
    -0.3298094961231737245750613333333E-29,
    +0.7763931143074065031714133333333E-30,
    -0.1841031343661458478421333333333E-30,
    +0.4395880138594310737100799999999E-31 };
  static double bm02cs[40] = {
    +0.9500415145228381369330861335560E-01,
    -0.3801864682365670991748081566851E-03,
    +0.2258339301031481192951829927224E-05,
    -0.3895725802372228764730621412605E-07,
    +0.1246886416512081697930990529725E-08,
    -0.6065949022102503779803835058387E-10,
    +0.4008461651421746991015275971045E-11,
    -0.3350998183398094218467298794574E-12,
    +0.3377119716517417367063264341996E-13,
    -0.3964585901635012700569356295823E-14,
    +0.5286111503883857217387939744735E-15,
    -0.7852519083450852313654640243493E-16,
    +0.1280300573386682201011634073449E-16,
    -0.2263996296391429776287099244884E-17,
    +0.4300496929656790388646410290477E-18,
    -0.8705749805132587079747535451455E-19,
    +0.1865862713962095141181442772050E-19,
    -0.4210482486093065457345086972301E-20,
    +0.9956676964228400991581627417842E-21,
    -0.2457357442805313359605921478547E-21,
    +0.6307692160762031568087353707059E-22,
    -0.1678773691440740142693331172388E-22,
    +0.4620259064673904433770878136087E-23,
    -0.1311782266860308732237693402496E-23,
    +0.3834087564116302827747922440276E-24,
    -0.1151459324077741271072613293576E-24,
    +0.3547210007523338523076971345213E-25,
    -0.1119218385815004646264355942176E-25,
    +0.3611879427629837831698404994257E-26,
    -0.1190687765913333150092641762463E-26,
    +0.4005094059403968131802476449536E-27,
    -0.1373169422452212390595193916017E-27,
    +0.4794199088742531585996491526437E-28,
    -0.1702965627624109584006994476452E-28,
    +0.6149512428936330071503575161324E-29,
    -0.2255766896581828349944300237242E-29,
    +0.8399707509294299486061658353200E-30,
    -0.3172997595562602355567423936152E-30,
    +0.1215205298881298554583333026514E-30,
    -0.4715852749754438693013210568045E-31 };
  static double bt02cs[39] = {
    -0.24548295213424597462050467249324,
    +0.12544121039084615780785331778299E-02,
    -0.31253950414871522854973446709571E-04,
    +0.14709778249940831164453426969314E-05,
    -0.99543488937950033643468850351158E-07,
    +0.85493166733203041247578711397751E-08,
    -0.86989759526554334557985512179192E-09,
    +0.10052099533559791084540101082153E-09,
    -0.12828230601708892903483623685544E-10,
    +0.17731700781805131705655750451023E-11,
    -0.26174574569485577488636284180925E-12,
    +0.40828351389972059621966481221103E-13,
    -0.66751668239742720054606749554261E-14,
    +0.11365761393071629448392469549951E-14,
    -0.20051189620647160250559266412117E-15,
    +0.36497978794766269635720591464106E-16,
    -0.68309637564582303169355843788800E-17,
    +0.13107583145670756620057104267946E-17,
    -0.25723363101850607778757130649599E-18,
    +0.51521657441863959925267780949333E-19,
    -0.10513017563758802637940741461333E-19,
    +0.21820381991194813847301084501333E-20,
    -0.46004701210362160577225905493333E-21,
    +0.98407006925466818520953651199999E-22,
    -0.21334038035728375844735986346666E-22,
    +0.46831036423973365296066286933333E-23,
    -0.10400213691985747236513382399999E-23,
    +0.23349105677301510051777740800000E-24,
    -0.52956825323318615788049749333333E-25,
    +0.12126341952959756829196287999999E-25,
    -0.28018897082289428760275626666666E-26,
    +0.65292678987012873342593706666666E-27,
    -0.15337980061873346427835733333333E-27,
    +0.36305884306364536682359466666666E-28,
    -0.86560755713629122479172266666666E-29,
    +0.20779909972536284571238399999999E-29,
    -0.50211170221417221674325333333333E-30,
    +0.12208360279441714184191999999999E-30,
    -0.29860056267039913454250666666666E-31 };
  static double bth0cs[44] = {
    -0.24901780862128936717709793789967,
    +0.48550299609623749241048615535485E-03,
    -0.54511837345017204950656273563505E-05,
    +0.13558673059405964054377445929903E-06,
    -0.55691398902227626227583218414920E-08,
    +0.32609031824994335304004205719468E-09,
    -0.24918807862461341125237903877993E-10,
    +0.23449377420882520554352413564891E-11,
    -0.26096534444310387762177574766136E-12,
    +0.33353140420097395105869955014923E-13,
    -0.47890000440572684646750770557409E-14,
    +0.75956178436192215972642568545248E-15,
    -0.13131556016891440382773397487633E-15,
    +0.24483618345240857495426820738355E-16,
    -0.48805729810618777683256761918331E-17,
    +0.10327285029786316149223756361204E-17,
    -0.23057633815057217157004744527025E-18,
    +0.54044443001892693993017108483765E-19,
    -0.13240695194366572724155032882385E-19,
    +0.33780795621371970203424792124722E-20,
    -0.89457629157111779003026926292299E-21,
    +0.24519906889219317090899908651405E-21,
    -0.69388422876866318680139933157657E-22,
    +0.20228278714890138392946303337791E-22,
    -0.60628500002335483105794195371764E-23,
    +0.18649748964037635381823788396270E-23,
    -0.58783732384849894560245036530867E-24,
    +0.18958591447999563485531179503513E-24,
    -0.62481979372258858959291620728565E-25,
    +0.21017901684551024686638633529074E-25,
    -0.72084300935209253690813933992446E-26,
    +0.25181363892474240867156405976746E-26,
    -0.89518042258785778806143945953643E-27,
    +0.32357237479762298533256235868587E-27,
    -0.11883010519855353657047144113796E-27,
    +0.44306286907358104820579231941731E-28,
    -0.16761009648834829495792010135681E-28,
    +0.64292946921207466972532393966088E-29,
    -0.24992261166978652421207213682763E-29,
    +0.98399794299521955672828260355318E-30,
    -0.39220375242408016397989131626158E-30,
    +0.15818107030056522138590618845692E-30,
    -0.64525506144890715944344098365426E-31,
    +0.26611111369199356137177018346367E-31 };
  double eta;
  static int nbm0 = 0;
  static int nbm02 = 0;
  static int nbt02 = 0;
  static int nbth0 = 0;
  static double pi4 = 0.785398163397448309615660845819876;
  static double xmax = 0.0;
  double z;

  if ( nbm0 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nbm0 = r8_inits ( bm0cs, 37, eta );
    nbt02 = r8_inits ( bt02cs, 39, eta );
    nbm02 = r8_inits ( bm02cs, 40, eta );
    nbth0 = r8_inits ( bth0cs, 44, eta );
    xmax = 1.0 / r8_mach ( 4 );
  }

  if ( x < 4.0 )
  {
    cerr << "\n";
    cerr << "R8_B0MP - Fatal error!\n";
    cerr << "  X < 4.\n";
    exit ( 1 );
  }
  else if ( x <= 8.0 )
  {
    z = ( 128.0 / x / x - 5.0 ) / 3.0;
    ampl = ( 0.75 + r8_csevl ( z, bm0cs, nbm0 ) ) / sqrt ( x );
    theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x;
  }
  else
  {
    z = 128.0 / x / x - 1.0;
    ampl = ( 0.75 + r8_csevl ( z, bm02cs, nbm02) ) / sqrt ( x );
    theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x;
  }
  return;
}
//****************************************************************************80

double r8_besi1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI1, the Bessel function I of order 1 of X.
//
{
  static double bi1cs[17] = {
    -0.19717132610998597316138503218149E-02,
    +0.40734887667546480608155393652014,
    +0.34838994299959455866245037783787E-01,
    +0.15453945563001236038598401058489E-02,
    +0.41888521098377784129458832004120E-04,
    +0.76490267648362114741959703966069E-06,
    +0.10042493924741178689179808037238E-07,
    +0.99322077919238106481371298054863E-10,
    +0.76638017918447637275200171681349E-12,
    +0.47414189238167394980388091948160E-14,
    +0.24041144040745181799863172032000E-16,
    +0.10171505007093713649121100799999E-18,
    +0.36450935657866949458491733333333E-21,
    +0.11205749502562039344810666666666E-23,
    +0.29875441934468088832000000000000E-26,
    +0.69732310939194709333333333333333E-29,
    +0.14367948220620800000000000000000E-31 };
  static int nti1 = 0;
  double  value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nti1 == 0 )
  {
    nti1 = r8_inits ( bi1cs, 17, 0.1 * r8_mach ( 3 ) );
    xmin = 2.0 * r8_mach ( 1 );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
    xmax = log ( r8_mach ( 2 ) );
  }

  y = r8_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0;
  }
  else if ( y <= xsml )
  {
    value = 0.5 * x;
  }
  else if ( y <= 3.0 )
  {
    value = x * ( 0.875 + r8_csevl ( y * y / 4.5 - 1.0, bi1cs, nti1 ) );
  }
  else if ( y <= xmax )
  {
    value = exp ( y ) * r8_besi1e ( x );
  }
  else
  {
    cerr << "\n";
    cerr << "R8_BESI1 - Fatal error!\n";
    cerr << "  Result overflows.\n";
    exit ( 1 );
  }
  return value;
}
//****************************************************************************80

double r8_besi1e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESI1E, the exponentially scaled Bessel 
//    function I1(X).
//
{
  static double ai12cs[69] = {
    +0.2857623501828012047449845948469E-01,
    -0.9761097491361468407765164457302E-02,
    -0.1105889387626237162912569212775E-03,
    -0.3882564808877690393456544776274E-05,
    -0.2512236237870208925294520022121E-06,
    -0.2631468846889519506837052365232E-07,
    -0.3835380385964237022045006787968E-08,
    -0.5589743462196583806868112522229E-09,
    -0.1897495812350541234498925033238E-10,
    +0.3252603583015488238555080679949E-10,
    +0.1412580743661378133163366332846E-10,
    +0.2035628544147089507224526136840E-11,
    -0.7198551776245908512092589890446E-12,
    -0.4083551111092197318228499639691E-12,
    -0.2101541842772664313019845727462E-13,
    +0.4272440016711951354297788336997E-13,
    +0.1042027698412880276417414499948E-13,
    -0.3814403072437007804767072535396E-14,
    -0.1880354775510782448512734533963E-14,
    +0.3308202310920928282731903352405E-15,
    +0.2962628997645950139068546542052E-15,
    -0.3209525921993423958778373532887E-16,
    -0.4650305368489358325571282818979E-16,
    +0.4414348323071707949946113759641E-17,
    +0.7517296310842104805425458080295E-17,
    -0.9314178867326883375684847845157E-18,
    -0.1242193275194890956116784488697E-17,
    +0.2414276719454848469005153902176E-18,
    +0.2026944384053285178971922860692E-18,
    -0.6394267188269097787043919886811E-19,
    -0.3049812452373095896084884503571E-19,
    +0.1612841851651480225134622307691E-19,
    +0.3560913964309925054510270904620E-20,
    -0.3752017947936439079666828003246E-20,
    -0.5787037427074799345951982310741E-22,
    +0.7759997511648161961982369632092E-21,
    -0.1452790897202233394064459874085E-21,
    -0.1318225286739036702121922753374E-21,
    +0.6116654862903070701879991331717E-22,
    +0.1376279762427126427730243383634E-22,
    -0.1690837689959347884919839382306E-22,
    +0.1430596088595433153987201085385E-23,
    +0.3409557828090594020405367729902E-23,
    -0.1309457666270760227845738726424E-23,
    -0.3940706411240257436093521417557E-24,
    +0.4277137426980876580806166797352E-24,
    -0.4424634830982606881900283123029E-25,
    -0.8734113196230714972115309788747E-25,
    +0.4045401335683533392143404142428E-25,
    +0.7067100658094689465651607717806E-26,
    -0.1249463344565105223002864518605E-25,
    +0.2867392244403437032979483391426E-26,
    +0.2044292892504292670281779574210E-26,
    -0.1518636633820462568371346802911E-26,
    +0.8110181098187575886132279107037E-28,
    +0.3580379354773586091127173703270E-27,
    -0.1692929018927902509593057175448E-27,
    -0.2222902499702427639067758527774E-28,
    +0.5424535127145969655048600401128E-28,
    -0.1787068401578018688764912993304E-28,
    -0.6565479068722814938823929437880E-29,
    +0.7807013165061145280922067706839E-29,
    -0.1816595260668979717379333152221E-29,
    -0.1287704952660084820376875598959E-29,
    +0.1114548172988164547413709273694E-29,
    -0.1808343145039336939159368876687E-30,
    -0.2231677718203771952232448228939E-30,
    +0.1619029596080341510617909803614E-30,
    -0.1834079908804941413901308439210E-31 };
  static double ai1cs[46] = {
    -0.2846744181881478674100372468307E-01,
    -0.1922953231443220651044448774979E-01,
    -0.6115185857943788982256249917785E-03,
    -0.2069971253350227708882823777979E-04,
    +0.8585619145810725565536944673138E-05,
    +0.1049498246711590862517453997860E-05,
    -0.2918338918447902202093432326697E-06,
    -0.1559378146631739000160680969077E-07,
    +0.1318012367144944705525302873909E-07,
    -0.1448423418183078317639134467815E-08,
    -0.2908512243993142094825040993010E-09,
    +0.1266388917875382387311159690403E-09,
    -0.1664947772919220670624178398580E-10,
    -0.1666653644609432976095937154999E-11,
    +0.1242602414290768265232168472017E-11,
    -0.2731549379672432397251461428633E-12,
    +0.2023947881645803780700262688981E-13,
    +0.7307950018116883636198698126123E-14,
    -0.3332905634404674943813778617133E-14,
    +0.7175346558512953743542254665670E-15,
    -0.6982530324796256355850629223656E-16,
    -0.1299944201562760760060446080587E-16,
    +0.8120942864242798892054678342860E-17,
    -0.2194016207410736898156266643783E-17,
    +0.3630516170029654848279860932334E-18,
    -0.1695139772439104166306866790399E-19,
    -0.1288184829897907807116882538222E-19,
    +0.5694428604967052780109991073109E-20,
    -0.1459597009090480056545509900287E-20,
    +0.2514546010675717314084691334485E-21,
    -0.1844758883139124818160400029013E-22,
    -0.6339760596227948641928609791999E-23,
    +0.3461441102031011111108146626560E-23,
    -0.1017062335371393547596541023573E-23,
    +0.2149877147090431445962500778666E-24,
    -0.3045252425238676401746206173866E-25,
    +0.5238082144721285982177634986666E-27,
    +0.1443583107089382446416789503999E-26,
    -0.6121302074890042733200670719999E-27,
    +0.1700011117467818418349189802666E-27,
    -0.3596589107984244158535215786666E-28,
    +0.5448178578948418576650513066666E-29,
    -0.2731831789689084989162564266666E-30,
    -0.1858905021708600715771903999999E-30,
    +0.9212682974513933441127765333333E-31,
    -0.2813835155653561106370833066666E-31 };
  static double bi1cs[17] = {
    -0.19717132610998597316138503218149E-02,
    +0.40734887667546480608155393652014,
    +0.34838994299959455866245037783787E-01,
    +0.15453945563001236038598401058489E-02,
    +0.41888521098377784129458832004120E-04,
    +0.76490267648362114741959703966069E-06,
    +0.10042493924741178689179808037238E-07,
    +0.99322077919238106481371298054863E-10,
    +0.76638017918447637275200171681349E-12,
    +0.47414189238167394980388091948160E-14,
    +0.24041144040745181799863172032000E-16,
    +0.10171505007093713649121100799999E-18,
    +0.36450935657866949458491733333333E-21,
    +0.11205749502562039344810666666666E-23,
    +0.29875441934468088832000000000000E-26,
    +0.69732310939194709333333333333333E-29,
    +0.14367948220620800000000000000000E-31 };
  double eta;
  static int ntai1 = 0;
  static int ntai12 = 0;
  static int nti1 = 0;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( nti1 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    nti1 = r8_inits ( bi1cs, 17, eta );
    ntai1 = r8_inits ( ai1cs, 46, eta );
    ntai12 = r8_inits ( ai12cs, 69, eta );
    xmin = 2.0 * r8_mach ( 1 );
    xsml = sqrt ( 8.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= xmin )
  {
    value = 0.0;
  }
  else if ( y <= xsml )
  {
    value = 0.5 * x * exp ( - y );
  }
  else if ( y <= 3.0 )
  {
    value = x * ( 0.875 + r8_csevl ( y * y / 4.5 - 1.0, bi1cs, nti1 ) )
      * exp ( - y );
  }
  else if ( y <= 8.0 )
  {
    value = ( 0.375 + r8_csevl ( ( 48.0 / y - 11.0) / 5.0, 
      ai1cs, ntai1 ) ) / sqrt ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  else
  {
    value = ( 0.375 + r8_csevl ( 16.0 / y - 1.0, ai12cs, ntai12 ) ) / sqrt ( y );
    if ( x < 0.0 )
    {
      value = - value;
    }
  }
  return value;
}
//****************************************************************************80

double r8_besj0 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESJ0, the Bessel function J of order 0 of X.
//
{
  double ampl;
  static double bj0cs[19] = {
    +0.10025416196893913701073127264074,
    -0.66522300776440513177678757831124,
    +0.24898370349828131370460468726680,
    -0.33252723170035769653884341503854E-01,
    +0.23114179304694015462904924117729E-02,
    -0.99112774199508092339048519336549E-04,
    +0.28916708643998808884733903747078E-05,
    -0.61210858663032635057818407481516E-07,
    +0.98386507938567841324768748636415E-09,
    -0.12423551597301765145515897006836E-10,
    +0.12654336302559045797915827210363E-12,
    -0.10619456495287244546914817512959E-14,
    +0.74706210758024567437098915584000E-17,
    -0.44697032274412780547627007999999E-19,
    +0.23024281584337436200523093333333E-21,
    -0.10319144794166698148522666666666E-23,
    +0.40608178274873322700800000000000E-26,
    -0.14143836005240913919999999999999E-28,
    +0.43910905496698880000000000000000E-31 };
  static int ntj0 = 0;
  double theta;
  double value;
  static double xsml = 0.0;
  double y;

  if ( ntj0 == 0 )
  {
    ntj0 = r8_inits ( bj0cs, 19, 0.1 * r8_mach ( 3 ) );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  y = r8_abs ( x );

  if ( y <= xsml )
  {
    value = 1.0;
  }
  else if ( y <= 4.0 )
  {
    value = r8_csevl ( 0.125 * y * y - 1.0, bj0cs, ntj0 );
  }
  else
  {
    r8_b0mp ( y, ampl, theta );
    value = ampl * cos ( theta );
  }
  return value;
}
//****************************************************************************80

double r8_besk ( double nu, double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK evaluates the Bessel function K of order NU of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double NU, the order.
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK, the Bessel function K of order NU at X.
//
{
  double *bke;
  int nin;
  double value;
  double xnu;

  xnu = nu - ( int ) ( nu );
  nin = ( int ) ( nu ) + 1;
  bke = r8_besks ( xnu, x, nin );

  value = bke[nin-1];

  delete [] bke;

  return value;
}
//****************************************************************************80

double r8_besk1 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK1, the Bessel function K of order 1 of X.
//
{
  static double bk1cs[16] = {
    +0.25300227338947770532531120868533E-01,
    -0.35315596077654487566723831691801,
    -0.12261118082265714823479067930042,
    -0.69757238596398643501812920296083E-02,
    -0.17302889575130520630176507368979E-03,
    -0.24334061415659682349600735030164E-05,
    -0.22133876307347258558315252545126E-07,
    -0.14114883926335277610958330212608E-09,
    -0.66669016941993290060853751264373E-12,
    -0.24274498505193659339263196864853E-14,
    -0.70238634793862875971783797120000E-17,
    -0.16543275155100994675491029333333E-19,
    -0.32338347459944491991893333333333E-22,
    -0.53312750529265274999466666666666E-25,
    -0.75130407162157226666666666666666E-28,
    -0.91550857176541866666666666666666E-31 };
  static int ntk1 = 0;
  double value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntk1 == 0 )
  {
    ntk1 = r8_inits ( bk1cs, 16, 0.1 * r8_mach ( 3 ) );
    xmin = exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax - 0.5 * xmax * log ( xmax ) 
      / ( xmax + 0.5 ) - 0.01;
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK1 = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = log ( 0.5 * x ) * r8_besi1 ( x ) + ( 0.75 
      + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = log ( 0.5 * x ) * r8_besi1 ( x ) + ( 0.75 
      + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x;
  }
  else if ( x <= xmax )
  {
    value = exp ( - x ) * r8_besk1e ( x );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_besk1e ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_BESK1E, the exponentially scaled Bessel 
//    function K1(X).
//
{
  static double ak12cs[33] = {
    +0.6379308343739001036600488534102E-01,
    +0.2832887813049720935835030284708E-01,
    -0.2475370673905250345414545566732E-03,
    +0.5771972451607248820470976625763E-05,
    -0.2068939219536548302745533196552E-06,
    +0.9739983441381804180309213097887E-08,
    -0.5585336140380624984688895511129E-09,
    +0.3732996634046185240221212854731E-10,
    -0.2825051961023225445135065754928E-11,
    +0.2372019002484144173643496955486E-12,
    -0.2176677387991753979268301667938E-13,
    +0.2157914161616032453939562689706E-14,
    -0.2290196930718269275991551338154E-15,
    +0.2582885729823274961919939565226E-16,
    -0.3076752641268463187621098173440E-17,
    +0.3851487721280491597094896844799E-18,
    -0.5044794897641528977117282508800E-19,
    +0.6888673850418544237018292223999E-20,
    -0.9775041541950118303002132480000E-21,
    +0.1437416218523836461001659733333E-21,
    -0.2185059497344347373499733333333E-22,
    +0.3426245621809220631645388800000E-23,
    -0.5531064394246408232501248000000E-24,
    +0.9176601505685995403782826666666E-25,
    -0.1562287203618024911448746666666E-25,
    +0.2725419375484333132349439999999E-26,
    -0.4865674910074827992378026666666E-27,
    +0.8879388552723502587357866666666E-28,
    -0.1654585918039257548936533333333E-28,
    +0.3145111321357848674303999999999E-29,
    -0.6092998312193127612416000000000E-30,
    +0.1202021939369815834623999999999E-30,
    -0.2412930801459408841386666666666E-31 };
  static double ak1cs[38] = {
    +0.27443134069738829695257666227266,
    +0.75719899531993678170892378149290E-01,
    -0.14410515564754061229853116175625E-02,
    +0.66501169551257479394251385477036E-04,
    -0.43699847095201407660580845089167E-05,
    +0.35402774997630526799417139008534E-06,
    -0.33111637792932920208982688245704E-07,
    +0.34459775819010534532311499770992E-08,
    -0.38989323474754271048981937492758E-09,
    +0.47208197504658356400947449339005E-10,
    -0.60478356628753562345373591562890E-11,
    +0.81284948748658747888193837985663E-12,
    -0.11386945747147891428923915951042E-12,
    +0.16540358408462282325972948205090E-13,
    -0.24809025677068848221516010440533E-14,
    +0.38292378907024096948429227299157E-15,
    -0.60647341040012418187768210377386E-16,
    +0.98324256232648616038194004650666E-17,
    -0.16284168738284380035666620115626E-17,
    +0.27501536496752623718284120337066E-18,
    -0.47289666463953250924281069568000E-19,
    +0.82681500028109932722392050346666E-20,
    -0.14681405136624956337193964885333E-20,
    +0.26447639269208245978085894826666E-21,
    -0.48290157564856387897969868800000E-22,
    +0.89293020743610130180656332799999E-23,
    -0.16708397168972517176997751466666E-23,
    +0.31616456034040694931368618666666E-24,
    -0.60462055312274989106506410666666E-25,
    +0.11678798942042732700718421333333E-25,
    -0.22773741582653996232867840000000E-26,
    +0.44811097300773675795305813333333E-27,
    -0.88932884769020194062336000000000E-28,
    +0.17794680018850275131392000000000E-28,
    -0.35884555967329095821994666666666E-29,
    +0.72906290492694257991679999999999E-30,
    -0.14918449845546227073024000000000E-30,
    +0.30736573872934276300799999999999E-31 };
  static double bk1cs[16] = {
    +0.25300227338947770532531120868533E-01,
    -0.35315596077654487566723831691801,
    -0.12261118082265714823479067930042,
    -0.69757238596398643501812920296083E-02,
    -0.17302889575130520630176507368979E-03,
    -0.24334061415659682349600735030164E-05,
    -0.22133876307347258558315252545126E-07,
    -0.14114883926335277610958330212608E-09,
    -0.66669016941993290060853751264373E-12,
    -0.24274498505193659339263196864853E-14,
    -0.70238634793862875971783797120000E-17,
    -0.16543275155100994675491029333333E-19,
    -0.32338347459944491991893333333333E-22,
    -0.53312750529265274999466666666666E-25,
    -0.75130407162157226666666666666666E-28,
    -0.91550857176541866666666666666666E-31 };
  double eta;
  static int ntak1 = 0;
  static int ntak12 = 0;
  static int ntk1 = 0;
  double value;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ntk1 == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntk1 = r8_inits ( bk1cs, 16, eta );
    ntak1 = r8_inits ( ak1cs, 38, eta );
    ntak12 = r8_inits ( ak12cs, 33, eta );
    xmin = exp ( r8_max ( log ( r8_mach ( 1 ) ), 
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    xsml = sqrt ( 4.0 * r8_mach ( 3 ) );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESK1E = Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }
  else if ( x <= xsml )
  {
    y = 0.0;
    value = exp ( x ) * ( log ( 0.5 * x ) * r8_besi1 ( x )
      + ( 0.75 + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 2.0 )
  {
    y = x * x;
    value = exp ( x ) * ( log ( 0.5 * x ) * r8_besi1 ( x )
      + ( 0.75 + r8_csevl ( 0.5 * y - 1.0, bk1cs, ntk1 ) ) / x );
  }
  else if ( x <= 8.0 )
  {
    value = ( 1.25 
      + r8_csevl ( ( 16.0 / x - 5.0 ) / 3.0, ak1cs, ntak1 ) ) / sqrt ( x );
  }
  else
  {
    value = ( 1.25 +
      r8_csevl ( 16.0 / x - 1.0, ak12cs, ntak12 ) ) / sqrt ( x );
  }
  return value;
}
//****************************************************************************80

double *r8_beskes ( double xnu, double x, int nin )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESKES: a sequence of exponentially scaled K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, ?
//    |XNU| < 1.
//
//    Input, double X, the argument.
//
//    Input, int NIN, indicates the number of terms to compute.
//
//    Output, double BESKES(abs(NIN)), the exponentially scaled 
//    K Bessel functions.
//
{
  double *bke;
  double bknu1;
  double direct;
  int i;
  int iswtch;
  int n;
  double v;
  double vend;
  double vincr;

  v = r8_abs ( xnu );
  n = i4_abs ( nin );

  if ( 1.0 <= v )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  |XNU| must be less than 1.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  if ( n == 0 )
  {
    cerr << "\n";
    cerr << "R8_BESKES - Fatal error!\n";
    cerr << "  N = 0.\n";
    exit ( 1 );
  }

  bke = new double[abs(nin)];

  r8_knus ( v, x, bke[0], bknu1, iswtch );

  if ( n == 1 )
  {
    return bke;
  }

  if ( nin < 0 )
  {
    vincr = - 1.0;
  }
  else
  {
    vincr = + 1.0;
  }

  if ( xnu < 0.0 )
  {
    direct = - vincr;
  }
  else
  {
    direct = vincr;
  }

  bke[1] = bknu1;

  if ( direct < 0.0 )
  {
    r8_knus ( r8_abs ( xnu + vincr ), x, bke[1], bknu1, iswtch );
  }

  if ( n == 2 )
  {
    return bke;
  }

  vend = r8_abs ( xnu + ( double ) ( nin ) ) - 1.0;

  v = xnu;
  for ( i = 3; i <= n; i++ )
  {
    v = v + vincr;
    bke[i-1] = 2.0 * v * bke[i-2] / x + bke[i-3];
  }
  return bke;
}
//****************************************************************************80

double *r8_besks ( double xnu, double x, int nin )

//****************************************************************************80
//
//  Purpose:
//
//    R8_BESKS evaluates a sequence of K Bessel functions at X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, ?
//    |XNU| < 1.
//
//    Input, double X, the argument.
//
//    Input, int NIN, indicates the number of terms to compute.
//
//    Output, double BESKS(abs(NIN)), the K Bessel functions.
//
{
  double *bk;
  double expxi;
  int i;
  int n;
  static double xmax = 0.0;

  if ( xmax == 0.0 )
  {
    xmax = - log ( r8_mach ( 1 ) );
    xmax = xmax + 0.5 * log ( 3.14 * 0.5 / xmax );
  }

  bk = r8_beskes ( xnu, x, nin );

  expxi = exp ( - x );
  n = i4_abs ( nin );

  for ( i = 0; i < n; i++ )
  {
    bk[i] = expxi * bk[i];
  }
  return bk;
}
//****************************************************************************80

 double r8_csevl ( double x, double a[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CSEVL evaluates a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, double CS[N], the Chebyshev coefficients.
//
//    Input, int N, the number of Chebyshev coefficients.
//
//    Output, double R8_CSEVL, the Chebyshev series evaluated at X.
//
{
  double b0;
  double b1;
  double b2;
  int i;
  double twox;
  double value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms <= 0.\n";
    exit ( 1 );
  }

  if ( 1000 < n )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  Number of terms greater than 1000.\n";
    exit ( 1 );
 }

  if ( x < -1.1 || 1.1 < x )
  {
    cerr << "\n";
    cerr << "R8_CSEVL - Fatal error!\n";
    cerr << "  X outside (-1,+1).\n";
    exit ( 1 );
  }

  twox = 2.0 * x;
  b1 = 0.0;
  b0 = 0.0;

  for ( i = n - 1; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[i];
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  static double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

void r8_gaml ( double &xmin, double &xmax )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAML evaluates bounds for an R8 argument of the gamma function.
//
//  Discussion:
//
//    This function calculates the minimum and maximum legal bounds 
//    for X in the evaluation of GAMMA ( X ).
//
//    XMIN and XMAX are not the only bounds, but they are the only 
//    non-trivial ones to calculate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Output, double &XMIN, &XMAX, the bounds.
//
{
  double alnbig;
  double alnsml;
  int i;
  int j;
  double xln;
  double xold;

  alnsml = log ( r8_mach ( 1 ) );
  xmin = - alnsml;

  for ( i = 1; i <= 10; i++ )
  {
    xold = xmin;
    xln = log ( xmin );
    xmin = xmin - xmin * ( ( xmin + 0.5 ) * xln - xmin 
      - 0.2258 + alnsml ) / ( xmin * xln + 0.5 );

    if ( r8_abs ( xmin - xold ) < 0.005 )
    {
      xmin = - xmin + 0.01;

      alnbig = log ( r8_mach ( 2 ) );
      xmax = alnbig;

      for ( j = 1; j <= 10; j++ )
      {
        xold = xmax;
        xln = log ( xmax );
        xmax = xmax - xmax * ( ( xmax - 0.5 ) * xln - xmax 
          + 0.9189 - alnbig ) / ( xmax * xln - 0.5 );

        if ( r8_abs ( xmax - xold ) < 0.005 )
        {
          xmax = xmax - 0.01;
          xmin = r8_max ( xmin, - xmax + 1.0 );
          return;
        }
      }
      cerr << "\n";
      cerr << "R8_GAML - Fatal error!\n";
      cerr << "  Unable to find XMAX.\n";
      exit ( 1 );
    }
  }
  cerr << "\n";
  cerr << "R8_GAML - Fatal error!\n";
  cerr << "  Unable to find XMIN.\n";
  exit ( 1 );
}
//****************************************************************************80

double r8_gamma ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_GAMMA evaluates the gamma function of an R8 argument.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_GAMMA, the gamma function of X.
//
{
  static double dxrel = 0.0;
  static double gcs[42] = {
    +0.8571195590989331421920062399942E-02,
    +0.4415381324841006757191315771652E-02,
    +0.5685043681599363378632664588789E-01,
    -0.4219835396418560501012500186624E-02,
    +0.1326808181212460220584006796352E-02,
    -0.1893024529798880432523947023886E-03,
    +0.3606925327441245256578082217225E-04,
    -0.6056761904460864218485548290365E-05,
    +0.1055829546302283344731823509093E-05,
    -0.1811967365542384048291855891166E-06,
    +0.3117724964715322277790254593169E-07,
    -0.5354219639019687140874081024347E-08,
    +0.9193275519859588946887786825940E-09,
    -0.1577941280288339761767423273953E-09,
    +0.2707980622934954543266540433089E-10,
    -0.4646818653825730144081661058933E-11,
    +0.7973350192007419656460767175359E-12,
    -0.1368078209830916025799499172309E-12,
    +0.2347319486563800657233471771688E-13,
    -0.4027432614949066932766570534699E-14,
    +0.6910051747372100912138336975257E-15,
    -0.1185584500221992907052387126192E-15,
    +0.2034148542496373955201026051932E-16,
    -0.3490054341717405849274012949108E-17,
    +0.5987993856485305567135051066026E-18,
    -0.1027378057872228074490069778431E-18,
    +0.1762702816060529824942759660748E-19,
    -0.3024320653735306260958772112042E-20,
    +0.5188914660218397839717833550506E-21,
    -0.8902770842456576692449251601066E-22,
    +0.1527474068493342602274596891306E-22,
    -0.2620731256187362900257328332799E-23,
    +0.4496464047830538670331046570666E-24,
    -0.7714712731336877911703901525333E-25,
    +0.1323635453126044036486572714666E-25,
    -0.2270999412942928816702313813333E-26,
    +0.3896418998003991449320816639999E-27,
    -0.6685198115125953327792127999999E-28,
    +0.1146998663140024384347613866666E-28,
    -0.1967938586345134677295103999999E-29,
    +0.3376448816585338090334890666666E-30,
    -0.5793070335782135784625493333333E-31 };
  int i;
  int n;
  static int ngcs = 0;
  static double pi = 3.14159265358979323846264338327950;
  double sinpiy;
  static double sq2pil = 0.91893853320467274178032973640562;
  double value;
  static double xmax = 0.0;
  static double xmin = 0.0;
  static double xsml = 0.0;
  double y;

  if ( ngcs == 0 )
  {
    ngcs = r8_inits ( gcs, 42, 0.1 * r8_mach ( 3 ) );
    r8_gaml ( xmin, xmax );
    xsml = exp ( r8_max ( log ( r8_mach ( 1 ) ),
      - log ( r8_mach ( 2 ) ) ) + 0.01 );
    dxrel = sqrt ( r8_mach ( 4 ) );
  }

  y = r8_abs ( x );

  if ( y <= 10.0 )
  {
    n = ( int ) ( x );
    if ( x < 0.0 )
    {
      n = n - 1;
    }
    y = x - ( double ) ( n );
    n = n - 1;
    value = 0.9375 + r8_csevl ( 2.0 * y - 1.0, gcs, ngcs );

    if ( n == 0 )
    {
      return value;
    }
    else if ( n < 0 )
    {
      n = - n;

      if ( x == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is 0.\n";
        exit ( 1 );
      }

      if ( x < 0.0 && x + ( double ) ( n - 2 ) == 0.0 )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is a negative int.\n";
        exit ( 1 );
      }

      if ( x < - 0.5 && r8_abs ( ( x - r8_aint ( x - 0.5 ) ) / x ) < dxrel )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Warning!\n";
        cerr << "  X too near a negative int,\n";
        cerr << "  answer is half precision.\n";
      }

      if ( y < xsml )
      {
        cerr << "\n";
        cerr << "R8_GAMMA - Fatal error!\n";
        cerr << "  X is so close to zero that Gamma overflows.\n";
        exit ( 1 );
      }

      for ( i = 1; i <= n; i++ )
      {
        value = value / ( x + ( double ) ( i - 1 ) );
      }

    }
    else if ( n == 0 )
    {
    }
    else
    {
      for ( i = 1; i <= n; i++ )
      {
        value = ( y + ( double ) ( i ) ) * value;
      }
    }
  }
  else
  {
    if ( xmax < x )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Fatal error!\n";
      cerr << "  X so big that Gamma overflows.\n";
      exit ( 1 );
    }
//
//  Underflow.
//
    if ( x < xmin )
    {
      value = 0.0;
      return value;
    }

    value = exp ( ( y - 0.5 ) * log ( y ) - y + sq2pil + r8_lgmc ( y ) );

    if ( 0.0 < x )
    {
      return value;
    }

    if ( r8_abs ( ( x - r8_aint ( x - 0.5 ) ) / x ) < dxrel )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Warning!\n";
      cerr << "  X too near a negative int,\n";
      cerr << "  answer is half precision.\n";
    }

    sinpiy = sin ( pi * y );

    if ( sinpiy == 0.0 )
    {
      cerr << "\n";
      cerr << "R8_GAMMA - Fatal error!\n";
      cerr << "  X is a negative int.\n";
      exit ( 1 );
    }
    value = - pi / ( y * sinpiy * value );
  }
  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

int r8_inits ( double dos[], int nos, double eta )

//****************************************************************************80
//
//  Purpose:
//
//    R8_INITS initializes a Chebyshev series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Roger Broucke,
//    Algorithm 446:
//    Ten Subroutines for the Manipulation of Chebyshev Series,
//    Communications of the ACM,
//    Volume 16, Number 4, April 1973, pages 254-256.
//
//  Parameters:
//
//    Input, double DOS[NOS], the Chebyshev coefficients.
//
//    Input, int NOS, the number of coefficients.
//
//    Input, double ETA, the desired accuracy.
//
//    Output, int R8_INITS, the number of terms of the series needed
//    to ensure the requested accuracy.
//
{
  double err;
  int i;
  int value;

  if ( nos < 1 )
  {
    cerr << "\n";
    cerr << "R8_INITS - Fatal error!\n";
    cerr << "  Number of coefficients < 1.\n";
    exit ( 1 );
  }

  err = 0.0;

  for ( i = nos - 1; 0 <= i; i-- )
  {
    err = err + r8_abs ( dos[i] );
    if ( eta < err )
    {
      value = i + 1;
      return value;
    }
  }

  value = i;
  cerr << "\n";
  cerr << "R8_INITS - Warning!\n";
  cerr << "  ETA may be too small.\n";

  return value;
}
//****************************************************************************80

void r8_knus ( double xnu, double x, double &bknu, double &bknu1, int &iswtch )

//****************************************************************************80
//
//  Purpose:
//
//    R8_KNUS computes a sequence of K Bessel functions.
//
//  Discussion:
//
//    This routine computes Bessel functions 
//      exp(x) * k-sub-xnu (x)  
//    and
//      exp(x) * k-sub-xnu+1 (x) 
//    for 0.0 <= xnu < 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double XNU, the order parameter.
//
//    Input, double X, the argument.
//
//    Output, double &BKNU, &BKNU1, the two K Bessel functions.
//
//    Output, int &ISWTCH, ?
//
{
  double a[32];
  double a0;
  static double aln2 = 0.69314718055994530941723212145818;
  static double alnbig = 0;
  static double alneps = 0;
  static double alnsml = 0;
  double alnz;
  double alpha[32];
  double an;
  double b0;
  double beta[32];
  double bknu0;
  double bknud;
  double bn;
  double c0;
  static double c0kcs[29] = {
    +0.60183057242626108387577445180329E-01,
    -0.15364871433017286092959755943124,
    -0.11751176008210492040068229226213E-01,
    -0.85248788891979509827048401550987E-03,
    -0.61329838767496791874098176922111E-04,
    -0.44052281245510444562679889548505E-05,
    -0.31631246728384488192915445892199E-06,
    -0.22710719382899588330673771793396E-07,
    -0.16305644608077609552274620515360E-08,
    -0.11706939299414776568756044043130E-09,
    -0.84052063786464437174546593413792E-11,
    -0.60346670118979991487096050737198E-12,
    -0.43326960335681371952045997366903E-13,
    -0.31107358030203546214634697772237E-14,
    -0.22334078226736982254486133409840E-15,
    -0.16035146716864226300635791528610E-16,
    -0.11512717363666556196035697705305E-17,
    -0.82657591746836959105169479089258E-19,
    -0.59345480806383948172333436695984E-20,
    -0.42608138196467143926499613023976E-21,
    -0.30591266864812876299263698370542E-22,
    -0.21963541426734575224975501815516E-23,
    -0.15769113261495836071105750684760E-24,
    -0.11321713935950320948757731048056E-25,
    -0.81286248834598404082792349714433E-27,
    -0.58360900893453226552829349315949E-28,
    -0.41901241623610922519452337780905E-29,
    -0.30083737960206435069530504212862E-30,
    -0.21599152067808647728342168089832E-31 };
  double eta;
  static double euler = 0.57721566490153286060651209008240;
  double expx;
  int i;
  int ii;
  int inu;
  int n;
  static int ntc0k = 0;
  int nterms;
  static int ntznu1 = 0;
  double p1;
  double p2;
  double p3;
  double qq;
  double result;
  static double sqpi2 = +1.2533141373155002512078826424055;
  double sqrtx;
  double v;
  double vlnz;
  double x2n;
  double x2tov;
  double xi;
  double xmu;
  static double xnusml = 0.0;
  static double xsml = 0.0;
  double z;
  static double znu1cs[20] = {
    +0.203306756994191729674444001216911,
    +0.140077933413219771062943670790563,
    +0.791679696100161352840972241972320E-02,
    +0.339801182532104045352930092205750E-03,
    +0.117419756889893366664507228352690E-04,
    +0.339357570612261680333825865475121E-06,
    +0.842594176976219910194629891264803E-08,
    +0.183336677024850089184748150900090E-09,
    +0.354969844704416310863007064469557E-11,
    +0.619032496469887332205244342078407E-13,
    +0.981964535680439424960346115456527E-15,
    +0.142851314396490474211473563005985E-16,
    +0.191894921887825298966162467488436E-18,
    +0.239430979739498914162313140597128E-20,
    +0.278890246815347354835870465474995E-22,
    +0.304606650633033442582845214092865E-24,
    +0.313173237042191815771564260932089E-26,
    +0.304133098987854951645174908005034E-28,
    +0.279840384636833084343185097659733E-30,
    +0.244637186274497596485238794922666E-32 };
  double ztov;

  if ( ntc0k == 0 )
  {
    eta = 0.1 * r8_mach ( 3 );
    ntc0k = r8_inits ( c0kcs, 29, eta );
    ntznu1 = r8_inits ( znu1cs, 20, eta );
    xnusml = sqrt ( r8_mach ( 3 ) / 8.0 );
    xsml = 0.1 * r8_mach ( 3 );
    alnsml = log ( r8_mach ( 1 ) );
    alnbig = log ( r8_mach ( 2 ) );
    alneps = log ( 0.1 * r8_mach ( 3 ) );
  }

  if ( xnu < 0.0 || 1.0 <= xnu )
  {
    cerr << "\n";
    cerr << "R8_KNUS - Fatal error!\n";
    cerr << "  XNU < 0 or 1 <= XNU.\n";
    exit ( 1 );
  }

  if ( x <= 0.0 )
  {
    cerr << "\n";
    cerr << "R8_KNUS - Fatal error!\n";
    cerr << "  X <= 0.\n";
    exit ( 1 );
  }

  iswtch = 0;
//
//  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
//  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
//  then to (0., .5), because k of negative order (-nu) = k of positive
//  order (+nu).
//
  if ( x <= 2.0 )
  {
    if ( xnu <= 0.5 )
    {
      v = xnu;
    }
    else
    {
      v = 1.0 - xnu;
    }
//
//  carefully find (x/2)^xnu and z^xnu where z = x*x/4.
//
    alnz = 2.0 * ( log ( x ) - aln2 );

    if ( x <= xnu )
    {
      if ( alnbig < - 0.5 * xnu * alnz - aln2 - log ( xnu ) )
      {
        cerr << "\n";
        cerr << "R8_KNUS - Fatal error!\n";
        cerr << "  Small X causing overflow.\n";
        exit ( 1 );
      }
    }

    vlnz = v * alnz;
    x2tov = exp ( 0.5 * vlnz );

    if ( vlnz <= alnsml )
    {
      ztov = 0.0;
    }
    else
    {
      ztov = x2tov * x2tov;
    }

    a0 = 0.5 * r8_gamma ( 1.0 + v );
    b0 = 0.5 * r8_gamma ( 1.0 - v );
    c0 = - euler;

    if ( 0.5 <= ztov && xnusml < v )
    {
      c0 = - 0.75 + r8_csevl ( ( 8.0 * v ) * v - 1.0, c0kcs, ntc0k );
    }

    if ( ztov <= 0.5 )
    {
      alpha[0] = ( a0 - ztov * b0 ) / v;
    }
    else
    {
      alpha[0] = c0 - alnz * ( 0.75 +
        r8_csevl ( vlnz / 0.35 + 1.0, znu1cs, ntznu1 ) ) * b0;
    }

    beta[0] = - 0.5 * ( a0 + ztov * b0 );

    if ( x <= xsml )
    {
      z = 0.0;
    }
    else
    {
      z = 0.25 * x * x;
    }

    nterms = i4_max ( 2, ( int ) ( 11.0 
      + ( 8.0 * alnz - 25.19 - alneps ) / ( 4.28 - alnz ) ) );

    for ( i = 2; i <= nterms; i++ )
    {
      xi = ( double ) ( i - 1 );
      a0 = a0 / ( xi * ( xi - v ) );
      b0 = b0 / ( xi * ( xi + v ) );
      alpha[i-1] = ( alpha[i-2] + 2.0 * xi * a0 ) 
        / ( xi * ( xi + v ) );
      beta[i-1] = ( xi - 0.5 * v ) * alpha[i-1] - ztov * b0;
    }

    bknu = alpha[nterms-1];
    bknud = beta[nterms-1];
    for ( ii = 2; ii <= nterms; ii++ )
    {
      i = nterms + 1 - ii;
      bknu = alpha[i-1] + bknu * z;
      bknud = beta[i-1] + bknud * z;
    }

    expx = exp ( x );
    bknu = expx * bknu / x2tov;

    if ( alnbig < - 0.5 * ( xnu + 1.0 ) * alnz - 2.0 * aln2 )
    {
      iswtch = 1;
      return;
    }

    bknud = expx * bknud * 2.0 / ( x2tov * x );

    if ( xnu <= 0.5 )
    {
      bknu1 = v * bknu / x - bknud;
      return;
    }
    bknu0 = bknu;
    bknu = - v * bknu / x - bknud;
    bknu1 = 2.0 * xnu * bknu / x + bknu0;
  }
//
//  x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
//  rational expansion.
//
  else
  {
    sqrtx = sqrt ( x );

    if ( 1.0 / xsml < x )
    {
      bknu = sqpi2 / sqrtx;
      bknu1 = bknu;
      return;
    }

    an = - 0.60 - 1.02 / x;
    bn = - 0.27 - 0.53 / x;
    nterms = i4_min ( 32, i4_max ( 3, ( int ) ( an + bn * alneps ) ) );

    for ( inu = 1; inu <= 2; inu++ )
    {
      if ( inu == 1 )
      {
        if ( xnu <= xnusml )
        {
          xmu = 0.0;
        }
        else
        {
          xmu = ( 4.0 * xnu ) * xnu;
        }
      }
      else
      {
        xmu = 4.0 * ( r8_abs ( xnu ) + 1.0 ) * ( r8_abs ( xnu ) + 1.0 );
      }

      a[0] = 1.0 - xmu;
      a[1] = 9.0 - xmu;
      a[2] = 25.0 - xmu;

      if ( a[1] == 0.0 )
      {
        result = sqpi2 * ( 16.0 * x + xmu + 7.0 ) / ( 16.0 * x * sqrtx );
      }
      else
      {
        alpha[0] = 1.0;
        alpha[1] = ( 16.0 * x + a[1] ) / a[1];
        alpha[2] = ( ( 768.0 * x + 48.0 * a[2] ) * x 
          + a[1] * a[2] ) / ( a[1] * a[2] );

        beta[0] = 1.0;
        beta[1] = ( 16.0 * x + ( xmu + 7.0 ) ) / a[1];
        beta[2] = ( ( 768.0 * x + 48.0 * ( xmu + 23.0 ) ) * x +
          ( ( xmu + 62.0 ) * xmu + 129.0 ) ) / ( a[1] * a[2] );

        for ( i = 4; i <= nterms; i++ )
        {
          n = i - 1;
          x2n = ( double ) ( 2 * n - 1 );

          a[i-1] = ( x2n + 2.0 ) * ( x2n + 2.0 ) - xmu;
          qq = 16.0 * x2n / a[i-1];
          p1 = - x2n * ( ( double ) ( 12 * n * n - 20 * n ) - a[0] ) 
            / ( ( x2n - 2.0 ) * a[i-1] ) - qq * x;
          p2 = ( ( double ) ( 12 * n * n - 28 * n + 8 ) - a[0] ) 
            / a[i-1] - qq * x;
          p3 = - x2n * a[i-4] / ( ( x2n - 2.0 ) * a[i-1] );

          alpha[i-1] = - p1 * alpha[i-2]
                       - p2 * alpha[i-3] 
                       - p3 * alpha[i-4];

          beta[i-1] =  - p1 * beta[i-2]
                       - p2 * beta[i-3] 
                       - p3 * beta[i-4];

        }
        result = sqpi2 * beta[nterms-1] / ( sqrtx * alpha[nterms-1] );
      }

      if ( inu == 1 )
      {
        bknu = result;
      }
      else
      {
        bknu1 = result;
      }
    }
  }
  return;
}
//****************************************************************************80

double r8_lgmc ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_LGMC evaluates the log gamma correction factor for an R8 argument.
//
//  Discussion:
//
//    For 10 <= X, compute the log gamma correction factor so that
//
//      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
//                          + ( x - 0.5 ) * log ( x ) - x 
//                          + r8_lgmc ( x )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 September 2011
//
//  Author:
//
//    Original FORTRAN77 version by Wayne Fullerton.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Wayne Fullerton,
//    Portable Special Function Routines,
//    in Portability of Numerical Software,
//    edited by Wayne Cowell,
//    Lecture Notes in Computer Science, Volume 57,
//    Springer 1977,
//    ISBN: 978-3-540-08446-4,
//    LC: QA297.W65.
//
//  Parameters:
//
//    Input, double X, the argument.
//
//    Output, double R8_LGMC, the correction factor.
//
{
  static double algmcs[15] = {
    +0.1666389480451863247205729650822,
    -0.1384948176067563840732986059135E-04,
    +0.9810825646924729426157171547487E-08,
    -0.1809129475572494194263306266719E-10,
    +0.6221098041892605227126015543416E-13,
    -0.3399615005417721944303330599666E-15,
    +0.2683181998482698748957538846666E-17,
    -0.2868042435334643284144622399999E-19,
    +0.3962837061046434803679306666666E-21,
    -0.6831888753985766870111999999999E-23,
    +0.1429227355942498147573333333333E-24,
    -0.3547598158101070547199999999999E-26,
    +0.1025680058010470912000000000000E-27,
    -0.3401102254316748799999999999999E-29,
    +0.1276642195630062933333333333333E-30 };
  static int nalgm = 0;
  double value;
  static double xbig = 0.0;
  static double xmax = 0.0;

  if ( nalgm == 0 )
  {
    nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) );
    xbig = 1.0 / sqrt ( r8_mach ( 3 ) );
    xmax = exp ( r8_min ( log ( r8_mach ( 2 ) / 12.0 ), 
      - log ( 12.0 * r8_mach ( 1 ) ) ) );
  }

  if ( x < 10.0 )
  {
    cerr << "\n";
    cerr << "R8_LGMC - Fatal error!\n";
    cerr << "  X must be at least 10.\n";
    exit ( 1 );
  }
  else if ( x < xbig )
  {
    value = r8_csevl ( 2.0 * ( 10.0 / x ) 
      * ( 10.0 / x ) - 1.0, algmcs, nalgm ) / x;
  }
  else if ( x < xmax )
  {
    value = 1.0 / ( 12.0 * x );
  }
  else
  {
    value = 0.0;
  }
  return value;
}
//****************************************************************************80

double r8_mach ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MACH returns double precision real machine constants.
//
//  Discussion:
//
//    Assuming that the internal representation of a double precision real
//    number is in base B, with T the number of base-B digits in the mantissa,
//    and EMIN the smallest possible exponent and EMAX the largest possible 
//    exponent, then
//
//      R8_MACH(1) = B^(EMIN-1), the smallest positive magnitude.
//      R8_MACH(2) = B^EMAX*(1-B^(-T)), the largest magnitude.
//      R8_MACH(3) = B^(-T), the smallest relative spacing.
//      R8_MACH(4) = B^(1-T), the largest relative spacing.
//      R8_MACH(5) = log10(B).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 April 2007
//
//  Author:
//
//    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Phyllis Fox, Andrew Hall, Norman Schryer,
//    Algorithm 528:
//    Framework for a Portable Library,
//    ACM Transactions on Mathematical Software,
//    Volume 4, Number 2, June 1978, page 176-188.
//
//  Parameters:
//
//    Input, int I, chooses the parameter to be returned.
//    1 <= I <= 5.
//
//    Output, double R8_MACH, the value of the chosen parameter.
//
{
  double value;

  if ( i == 1 )
  {
    value = 4.450147717014403E-308;
  }
  else if ( i == 2 )
  {
    value = 8.988465674311579E+307;
  }
  else if ( i == 3 )
  {
    value = 1.110223024625157E-016;
  }
  else if ( i == 4 )
  {
    value = 2.220446049250313E-016;
  }
  else if ( i == 5 )
  {
    value = 0.301029995663981E+000;
  }
  else
  {
    cerr << "\n";
    cerr << "R8_MACH - Fatal error!\n";
    cerr << "  The input argument I is out of bounds.\n";
    cerr << "  Legal values satisfy 1 <= I <= 5.\n";
    cerr << "  I = " << i << "\n";
    value = 0.0;
    exit ( 1 );
  }

  return value;
}
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

double r8_sign ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SIGN returns the sign of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the number whose sign is desired.
//
//    Output, double R8_SIGN, the sign of X.
//
{
  double value;

  if ( x < 0.0 )
  {
    value = -1.0;
  }
  else
  {
    value = 1.0;
  }
  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double *r8mat_cholesky_factor ( int n, double a[], int &flag )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The matrix must be symmetric and positive semidefinite.
//
//    For a positive semidefinite symmetric matrix A, the Cholesky factorization
//    is a lower triangular matrix L such that:
//
//      A = L * L'
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix A.
//
//    Input, double A[N*N], the N by N matrix.
//
//    Output, int &FLAG, an error flag.
//    0, no error occurred.
//    1, the matrix is not positive definite.
//    2, the matrix is not nonnegative definite.
//
//    Output, double R8MAT_CHOLESKY_FACTOR[N*N], the N by N lower triangular
//    Cholesky factor.
//
{
  double *c;
  int i;
  int j;
  int k;
  double sum2;
  double tol;

  flag = 0;
  tol = sqrt ( r8_epsilon ( ) );

  c = r8mat_copy_new ( n, n, a );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      c[i+j*n] = 0.0;
    }
    for ( i = j; i < n; i++ )
    {
      sum2 = c[j+i*n];
      for ( k = 0; k < j; k++ )
      {
        sum2 = sum2 - c[j+k*n] * c[i+k*n];
      }
      if ( i == j )
      {
        if ( 0.0 < sum2 )
        {
          c[i+j*n] = sqrt ( sum2 );
        }
        else if ( sum2 < - tol )
        {
          flag = 2;
          cerr << "\n";
          cerr << "R8MAT_CHOLESKY_FACTOR - Fatal error!\n";
          cerr << "  Matrix is not nonnegative definite.\n";
          cerr << "  Diagonal I = " << i << "\n";
          cerr << "  SUM2 = " << sum2 << "\n";
          exit ( 1 );
        }
        else
        {
          flag = 1;
          c[i+j*n] = 0.0;
        }
      }
      else
      {

        if ( c[j+j*n] != 0.0 )
        {
          c[i+j*n] = sum2 / c[j+j*n];
        }
        else
        {
          c[i+j*n] = 0.0;
        }
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8mat_copy_new ( int m, int n, double a1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's, which
//    may be stored as a vector in column-major order.
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
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A1[M*N], the matrix to be copied.
//
//    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
//
{
  double *a2;
  int i;
  int j;

  a2 = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      a2[i+j*m] = a1[i+j*m];
    }
  }
  return a2;
}
//****************************************************************************80

double r8mat_is_symmetric ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
//
//  Discussion:
//
//    An R8MAT is a matrix of double precision real values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Input, double A[M*N], the matrix.
//
//    Output, double RMAT_IS_SYMMETRIC, measures the 
//    Frobenius norm of ( A - A' ), which would be zero if the matrix
//    were exactly symmetric.
//
{
  int i;
  int j;
  double value;

  if ( m != n )
  {
    value = r8_huge ( );
    return value;
  }

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a[i+j*m] - a[j+i*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
}
//****************************************************************************80

double r8mat_max ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MAX returns the maximum entry of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Output, double R8MAT_MAX, the maximum entry of A.
//
{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( value < a[i+j*m] )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

double r8mat_min ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MIN returns the minimum entry of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Output, double DMIN_MAX, the minimum entry of A.
//
{
  int i;
  int j;
  double value;

  value = a[0+0*m];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < value )
      {
        value = a[i+j*m];
      }
    }
  }
  return value;
}
//****************************************************************************80

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

double *r8mat_normal_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.
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
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
//
{
  double *r;

  r = r8vec_normal_01_new ( m * n, seed );

  return r;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j - 1 << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i - 1 << ": ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
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

double *r8vec_normal_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, double Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  double pi = 3.141592653589793;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
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

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
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
//    19 August 2004
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
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

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *sample_paths_cholesky ( int n, int n2, double rhomax, double rho0, 
  double *correlation ( int n, double rho_vec[], double rho0 ), int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_PATHS_CHOLESKY: sample paths for stationary correlation functions.
//
//  Discussion:
//
//    This method uses the Cholesky factorization of the correlation matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points on each path.
//
//    Input, int N2, the number of paths.
//
//    Input, double RHOMAX, the maximum value of RHO.
//
//    Input, double RHO0, the correlation length.
//
//    Input, double *CORRELATION ( int n, double rho_vec[], double rho0), 
//    the name of the function which evaluates the correlation.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double X[N*N2], the sample paths.
//
{
  double *cor;
  double *cor_vec;
  int flag;
  int i;
  int j;
  int k;
  double *l;
  double *r;
  double *rho_vec;
  double rhomin;
  double *x;
//
//  Choose N equally spaced sample points from 0 to RHOMAX.
//
  rhomin = 0.0;
  rho_vec = r8vec_linspace_new ( n, rhomin, rhomax );
//
//  Evaluate the correlation function.
//
  cor_vec = correlation ( n, rho_vec, rho0 );
//
//  Construct the correlation matrix;
//
//  From the vector 
//    [ C(0), C(1), C(2), ... C(N-1) ]
//  construct the vector
//    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
//  Every row of the correlation matrix can be constructed by a subvector
//  of this vector.
//
  cor = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      k = i4_wrap ( j - i, 0, n - 1 );
      cor[i+j*n] = cor_vec[k];
    }
  }
//
//  Get the Cholesky factorization of COR:
//
//    COR = L * L'.
//
  l = r8mat_cholesky_factor ( n, cor, flag );
//
//  The matrix might not be nonnegative definite.
//
  if ( flag == 2 )
  {
    cerr << "\n";
    cerr << "SAMPLE_PATHS_CHOLESKY - Fatal error!\n";
    cerr << "  The correlation matrix is not\n";
    cerr << "  symmetric nonnegative definite.\n";
    exit ( 1 );
  }
//
//  Compute a matrix of N by N2 normally distributed values.
//
  r = r8mat_normal_01_new ( n, n2, seed );
//
//  Compute the sample path.
//
  x = r8mat_mm_new ( n, n, n2, l, r );

  delete [] cor;
  delete [] cor_vec;
  delete [] l;
  delete [] r;
  delete [] rho_vec;

  return x;
}
//****************************************************************************80

double *sample_paths_eigen ( int n, int n2, double rhomax, double rho0, 
  double *correlation ( int n, double rho_vec[], double rho0 ), int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_PATHS_EIGEN: sample paths for stationary correlation functions.
//
//  Discussion:
//
//    This method uses the eigen-decomposition of the correlation matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points on each path.
//
//    Input, int N2, the number of paths.
//
//    Input, double RHOMAX, the maximum value of RHO.
//
//    Input, double RHO0, the correlation length.
//
//    Input, double *CORRELATION ( int n, double rho_vec[], double rho0), 
//    the name of the function which evaluates the correlation.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double X[N*N2], the sample paths.
//
{
  double *c;
  double *cor;
  double *cor_vec;
  double *d;
  double dmin;
  int i;
  int ierr;
  int j;
  int k;
  double *r;
  double *rho_vec;
  double rhomin;
  double *v;
  double *w;
  double *x;
//
//  Choose N equally spaced sample points from 0 to RHOMAX.
//
  rhomin = 0.0;
  rho_vec = r8vec_linspace_new ( n, rhomin, rhomax );
//
//  Evaluate the correlation function.
//
  cor_vec = correlation ( n, rho_vec, rho0 );
//
//  Construct the correlation matrix;
//
//  From the vector 
//    [ C(0), C(1), C(2), ... C(N-1) ]
//  construct the vector
//    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
//  Every row of the correlation matrix can be constructed by a subvector
//  of this vector.
//
  cor = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      k = i4_wrap ( abs ( i - j ), 0, n - 1 );
      cor[i+j*n] = cor_vec[k];
    }
  }
//
//  Get the eigendecomposition of COR:
//
//    COR = V * D * V'.
//
//  Because COR is symmetric, V is orthogonal.
//
  d = new double[n];
  w = new double[n];
  v = new double[n*n];

  tred2 ( n, cor, d, w, v );

  ierr = tql2 ( n, d, w, v );
//
//  We assume COR is non-negative definite, and hence that there
//  are no negative eigenvalues.
//
  dmin = r8vec_min ( n, d );

  if ( dmin < - sqrt ( r8_epsilon ( ) ) )
  {
    cout << "\n";
    cout << "SAMPLE_PATHS_EIGEN - Warning!\n";
    cout << "  Negative eigenvalues observed as low as " << dmin << "\n";
  }

  for ( i = 0; i < n; i++ )
  {
    d[i] = r8_max ( d[i], 0.0 );
  }
//
//  Compute the eigenvalues of the factor C.
//
  for ( i = 0; i < n; i++ )
  {
    d[i] = sqrt ( d[i] );
  }
//
//  Compute C, such that C' * C = COR.
//
  c = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i+j*n] = c[i+j*n] + d[k] * v[i+k*n] * v[j+k*n];
      }
    }
  }
//
//  Compute N by N2 independent random normal values.
//
  r = r8mat_normal_01_new ( n, n2, seed );
//
//  Multiply to get the variables X which have correlation COR.
//
  x = r8mat_mm_new ( n, n, n2, c, r );

  delete [] c;
  delete [] cor;
  delete [] cor_vec;
  delete [] d;
  delete [] r;
  delete [] rho_vec;
  delete [] v;
  delete [] w;

  return x;
}
//****************************************************************************80

double *sample_paths2_cholesky ( int n, int n2, double rhomin, double rhomax, 
  double rho0, double *correlation2 ( int m, int n, double s[], double t[], 
  double rho0 ), int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_PATHS2_CHOLESKY: sample paths for stationary correlation functions.
//
//  Discussion:
//
//    This method uses the Cholesky factorization of the correlation matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points on each path.
//
//    Input, int N2, the number of paths.
//
//    Input, double RHOMIN, RHOMAX, the range of RHO.
//
//    Input, double RHO0, the correlation length.
//
//    Input, double *CORRELATION2 ( int m, int n, double s[], double t[], 
//    double rho0 ), the name of the function which evaluates the correlation.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double X[N*N2], the sample paths.
//
{
  double *cor;
  int flag;
  int i;
  int j;
  int k;
  double *l;
  double *r;
  double *s;
  double *x;
//
//  Choose N equally spaced sample points from RHOMIN to RHOMAX.
//
  s = r8vec_linspace_new ( n, rhomin, rhomax );
//
//  Evaluate the correlation function.
//
  cor = correlation2 ( n, n, s, s, rho0 );
//
//  Get the Cholesky factorization of COR:
//
//    COR = L * L'.
//
  l = r8mat_cholesky_factor ( n, cor, flag );
//
//  The matrix might not be nonnegative definite.
//
  if ( flag == 2 )
  {
    cerr << "\n";
    cerr << "SAMPLE_PATHS2_CHOLESKY - Fatal error!\n";
    cerr << "  The correlation matrix is not\n";
    cerr << "  symmetric nonnegative definite.\n";
    exit ( 1 );
  }
//
//  Compute a matrix of N by N2 normally distributed values.
//
  r = r8mat_normal_01_new ( n, n2, seed );
//
//  Compute the sample path.
//
  x = r8mat_mm_new ( n, n, n2, l, r );

  delete [] cor;
  delete [] l;
  delete [] r;
  delete [] s;

  return x;
}
//****************************************************************************80

double *sample_paths2_eigen ( int n, int n2, double rhomin, double rhomax, 
  double rho0, double *correlation2 ( int m, int n, double s[], double t[], 
  double rho0 ), int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_PATHS2_EIGEN: sample paths for stationary correlation functions.
//
//  Discussion:
//
//    This method uses the eigen-decomposition of the correlation matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points on each path.
//
//    Input, int N2, the number of paths.
//
//    Input, double RHOMIN, RHOMAX, the range of RHO.
//
//    Input, double RHO0, the correlation length.
//
//    Input, double *CORRELATION2 ( int m, int n, double s[], double t[], 
//    double rho0 ), the name of the function which evaluates the correlation.
//
//    Input/output, int &SEED, a seed for the random number
//    generator.
//
//    Output, double X[N*N2], the sample paths.
//
{
  double *c;
  double *cor;
  double *d;
  double dmin;
  int i;
  int ierr;
  int j;
  int k;
  double *r;
  double *s;
  double *v;
  double *w;
  double *x;
//
//  Choose N equally spaced sample points from RHOMIN to RHOMAX.
//
  s = r8vec_linspace_new ( n, rhomin, rhomax );
//
//  Evaluate the correlation function.
//
  cor = correlation2 ( n, n, s, s, rho0 );
//
//  Get the eigendecomposition of COR:
//
//    COR = V * D * V'.
//
//  Because COR is symmetric, V is orthogonal.
//
  d = new double[n];
  w = new double[n];
  v = new double[n*n];

  tred2 ( n, cor, d, w, v );

  ierr = tql2 ( n, d, w, v );
//
//  We assume COR is non-negative definite, and hence that there
//  are no negative eigenvalues.
//
  dmin = r8vec_min ( n, d );

  if ( dmin < - sqrt ( r8_epsilon ( ) ) )
  {
    cout << "\n";
    cout << "SAMPLE_PATHS2_EIGEN - Warning!\n";
    cout << "  Negative eigenvalues observed as low as " << dmin << "\n";
  }

  for ( i = 0; i < n; i++ )
  {
    d[i] = r8_max ( d[i], 0.0 );
  }
//
//  Compute the eigenvalues of the factor C.
//
  for ( i = 0; i < n; i++ )
  {
    d[i] = sqrt ( d[i] );
  }
//
//  Compute C, such that C' * C = COR.
//
  c = new double[n*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      c[i+j*n] = 0.0;
      for ( k = 0; k < n; k++ )
      {
        c[i+j*n] = c[i+j*n] + d[k] * v[i+k*n] * v[j+k*n];
      }
    }
  }
//
//  Compute N by N2 independent random normal values.
//
  r = r8mat_normal_01_new ( n, n2, seed );
//
//  Multiply to get the variables X which have correlation COR.
//
  x = r8mat_mm_new ( n, n, n2, c, r );

  delete [] c;
  delete [] cor;
  delete [] d;
  delete [] r;
  delete [] s;
  delete [] v;
  delete [] w;

  return x;
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

int tql2 ( int n, double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
//
//  Discussion:
//
//    This subroutine finds the eigenvalues and eigenvectors of a symmetric
//    tridiagonal matrix by the QL method.  The eigenvectors of a full
//    symmetric matrix can also be found if TRED2 has been used to reduce this
//    full matrix to tridiagonal form.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Bowdler, Martin, Reinsch, Wilkinson,
//    TQL2,
//    Numerische Mathematik,
//    Volume 11, pages 293-306, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double D[N].  On input, the diagonal elements of
//    the matrix.  On output, the eigenvalues in ascending order.  If an error
//    exit is made, the eigenvalues are correct but unordered for indices
//    1,2,...,IERR-1.
//
//    Input/output, double E[N].  On input, E(2:N) contains the
//    subdiagonal elements of the input matrix, and E(1) is arbitrary.
//    On output, E has been destroyed.
//
//    Input, double Z[N*N].  On input, the transformation matrix
//    produced in the reduction by TRED2, if performed.  If the eigenvectors of
//    the tridiagonal matrix are desired, Z must contain the identity matrix.
//    On output, Z contains the orthonormal eigenvectors of the symmetric
//    tridiagonal (or full) matrix.  If an error exit is made, Z contains
//    the eigenvectors associated with the stored eigenvalues.
//
//    Output, int TQL2, error flag.
//    0, normal return,
//    J, if the J-th eigenvalue has not been determined after
//    30 iterations.
//
{
  double c;
  double c2;
  double c3;
  double dl1;
  double el1;
  double f;
  double g;
  double h;
  int i;
  int ierr;
  int ii;
  int j;
  int k;
  int l;
  int l1;
  int l2;
  int m;
  int mml;
  double p;
  double r;
  double s;
  double s2;
  double t;
  double tst1;
  double tst2;

  ierr = 0;

  if ( n == 1 )
  {
    return ierr;
  }

  for ( i = 1; i < n; i++ )
  {
    e[i-1] = e[i];
  }

  f = 0.0;
  tst1 = 0.0;
  e[n-1] = 0.0;

  for ( l = 0; l < n; l++ )
  {
    j = 0;
    h = r8_abs ( d[l] ) + r8_abs ( e[l] );
    tst1 = r8_max ( tst1, h );
//
//  Look for a small sub-diagonal element.
//
    for ( m = l; m < n; m++ )
    {
      tst2 = tst1 + r8_abs ( e[m] );
      if ( tst2 == tst1 )
      {
        break;
      }
    }

    if ( m != l )
    {
      for ( ; ; )
      {
        if ( 30 <= j )
        {
          ierr = l + 1;
          return ierr;
        }

        j = j + 1;
//
//  Form shift.
//
        l1 = l + 1;
        l2 = l1 + 1;
        g = d[l];
        p = ( d[l1] - g ) / ( 2.0 * e[l] );
        r = pythag ( p, 1.0 );
        d[l] = e[l] / ( p + r8_sign ( p ) * r8_abs ( r ) );
        d[l1] = e[l] * ( p + r8_sign ( p ) * r8_abs ( r ) );
        dl1 = d[l1];
        h = g - d[l];
        for ( i = l2; i < n; i++ )
        {
          d[i] = d[i] - h;
        }
        f = f + h;
//
//  QL transformation.
//
        p = d[m];
        c = 1.0;
        c2 = c;
        el1 = e[l1];
        s = 0.0;
        mml = m - l;

        for ( ii = 1; ii <= mml; ii++ )
        {
          c3 = c2;
          c2 = c;
          s2 = s;
          i = m - ii;
          g = c * e[i];
          h = c * p;
          r = pythag ( p, e[i] );
          e[i+1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i+1] = h + s * ( c * g + s * d[i] );
//
//  Form vector.
//
          for ( k = 0; k < n; k++ )
          {
            h = z[k+(i+1)*n];
            z[k+(i+1)*n] = s * z[k+i*n] + c * h;
            z[k+i*n] = c * z[k+i*n] - s * h;
          }
        }
        p = - s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
        tst2 = tst1 + r8_abs ( e[l] );

        if ( tst2 <= tst1 )
        {
          break;
        }
      }
    }
    d[l] = d[l] + f;
  }
//
//  Order eigenvalues and eigenvectors.
//
  for ( ii = 1; ii < n; ii++ )
  {
    i = ii - 1;
    k = i;
    p = d[i];
    for ( j = ii; j < n; j++ )
    {
      if ( d[j] < p )
      {
        k = j;
        p = d[j];
      }
    }

    if ( k != i )
    {
      d[k] = d[i];
      d[i] = p;
      for ( j = 0; j < n; j++ )
      {
        t        = z[j+i*n];
        z[j+i*n] = z[j+k*n];
        z[j+k*n] = t;
      }
    }
  }
  return ierr;
}
//****************************************************************************80

void tred2 ( int n, double a[], double d[], double e[], double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
//
//  Discussion:
//
//    This subroutine reduces a real symmetric matrix to a
//    symmetric tridiagonal matrix using and accumulating
//    orthogonal similarity transformations.
//
//    A and Z may coincide, in which case a single storage area is used
//    for the input of A and the output of Z.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 November 2012
//
//  Author:
//
//    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
//    Klema, Moler.
//    C version by John Burkardt.
//
//  Reference:
//
//    Martin, Reinsch, Wilkinson,
//    TRED2,
//    Numerische Mathematik,
//    Volume 11, pages 181-195, 1968.
//
//    James Wilkinson, Christian Reinsch,
//    Handbook for Automatic Computation,
//    Volume II, Linear Algebra, Part 2,
//    Springer, 1971,
//    ISBN: 0387054146,
//    LC: QA251.W67.
//
//    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
//    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
//    Matrix Eigensystem Routines, EISPACK Guide,
//    Lecture Notes in Computer Science, Volume 6,
//    Springer Verlag, 1976,
//    ISBN13: 978-3540075462,
//    LC: QA193.M37.
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the real symmetric input matrix.  Only the
//    lower triangle of the matrix need be supplied.
//
//    Output, double D[N], the diagonal elements of the tridiagonal
//    matrix.
//
//    Output, double E[N], contains the subdiagonal elements of the
//    tridiagonal matrix in E(2:N).  E(1) is set to zero.
//
//    Output, double Z[N*N], the orthogonal transformation matrix
//    produced in the reduction.
//
{
  double f;
  double g;
  double h;
  double hh;
  int i;
  int ii;
  int j;
  int k;
  int l;
  double scale;

  for ( j = 0; j < n; j++ )
  {
    for ( i = j; i < n; i++ )
    {
      z[i+j*n] = a[i+j*n];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = a[n-1+j*n];
  }

  for ( i = n - 1; 1 <= i; i-- )
  {
    l = i - 1;
    h = 0.0;
//
//  Scale row.
//
    scale = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      scale = scale + r8_abs ( d[k] );
    }

    if ( scale == 0.0 )
    {
      e[i] = d[l];

      for ( j = 0; j <= l; j++ )
      {
        d[j]     = z[l+j*n];
        z[i+j*n] = 0.0;
        z[j+i*n] = 0.0;
      }
      d[i] = 0.0;
      continue;
    }

    for ( k = 0; k <= l; k++ )
    {
      d[k] = d[k] / scale;
    }

    h = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      h = h + d[k] * d[k];
    }

    f = d[l];
    g = - sqrt ( h ) * r8_sign ( f );
    e[i] = scale * g;
    h = h - f * g;
    d[l] = f - g;
//
//  Form A*U.
//
    for ( k = 0; k <= l; k++ )
    {
      e[k] = 0.0;
    }

    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      z[j+i*n] = f;
      g = e[j] + z[j+j*n] * f;

      for ( k = j + 1; k <= l; k++ )
      {
        g = g + z[k+j*n] * d[k];
        e[k] = e[k] + z[k+j*n] * f;
      }
      e[j] = g;
    }
//
//  Form P.
//
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] / h;
    }
    f = 0.0;
    for ( k = 0; k <= l; k++ )
    {
      f = f + e[k] * d[k];
    }
    hh = 0.5 * f / h;
//
//  Form Q.
//
    for ( k = 0; k <= l; k++ )
    {
      e[k] = e[k] - hh * d[k];
    }
//
//  Form reduced A.
//
    for ( j = 0; j <= l; j++ )
    {
      f = d[j];
      g = e[j];

      for ( k = j; k <= l; k++ )
      {
        z[k+j*n] = z[k+j*n] - f * e[k] - g * d[k];
      }
      d[j] = z[l+j*n];
      z[i+j*n] = 0.0;
    }
    d[i] = h;
  }
//
//  Accumulation of transformation matrices.
//
  for ( i = 1; i < n; i++ )
  {
    l = i - 1;
    z[n-1+l*n] = z[l+l*n];
    z[l+l*n] = 1.0;
    h = d[i];

    if ( h != 0.0 )
    {
      for ( k = 0; k <= l; k++ )
      {
        d[k] = z[k+i*n] / h;
      }
      for ( j = 0; j <= l; j++ )
      {
        g = 0.0;
        for ( k = 0; k <= l; k++ )
        {
          g = g + z[k+i*n] * z[k+j*n];
        }
        for ( k = 0; k <= l; k++ )
        {
          z[k+j*n] = z[k+j*n] - g * d[k];
        }
      }
    }
    for ( k = 0; k <= l; k++ )
    {
      z[k+i*n] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    d[j] = z[n-1+j*n];
  }

  for ( j = 0; j < n - 1; j++ )
  {
    z[n-1+j*n] = 0.0;
  }
  z[n-1+(n-1)*n] = 1.0;

  e[0] = 0.0;

  return;
}

