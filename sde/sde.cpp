# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sde.hpp"
# include "qr_solve.hpp"

//****************************************************************************80

double *bpath ( int &seed, int n )

//****************************************************************************80
//
//  Purpose:
//
//    BPATH performs a Brownian path simulation.
//
//  Discussion:
//
//    This routine computes one simulation of discretized Brownian 
//    motion over the time interval [0,1] using N time steps.
//    The user specifies a random number seed.  Different values of
//    the seed will result in different realizations of the path.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random number 
//    generator.
//
//    Input, int N, the number of steps.
//
//    Output, double BPATH[N+1], the Brownian path.
//
{
  double dt;
  double *dw;
  int j;
  double tmax;
  double *w;

  tmax = 1.0;
  dt = tmax / ( double ) ( n );
//
//  Define the increments dW.
//
  dw = r8vec_normal_01_new ( n, seed );

  for ( j = 0; j < n; j++ )
  {
    dw[j] = sqrt ( dt ) * dw[j];
  }
//
//  W is the sum of the previous increments.
//
  w = new double[n+1];

  w[0] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    w[j] = w[j-1] + dw[j];
  }

  delete [] dw;

  return w;
}
//****************************************************************************80

void bpath_gnuplot ( int n, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    BPATH_GNUPLOT writes a GNUPLOT input file to plot BPATH data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int N, the number of steps.
//
//    Input, double W[N+1], the Brownian path.
//
{
  string command_filename = "bpath_commands.txt";
  ofstream command;
  string data_filename = "bpath_data.txt";
  ofstream data;
  int i;
  double t;
//
//  Create the data file.
//
  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i <= n; i++ )
  {
    t = ( double ) ( i ) / ( double ) ( n );
    data << "  " << t
         << "  " << w[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  BPATH data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# bpath_commands.txt\n";
  command << "# created by sde::bpath_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < bpath_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'bpath.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel 'W(t)'\n";
  command << "set title 'Brownian motion by BPATH'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot 'bpath_data.txt' using 1:2\n";
  command << "quit\n";

  command.close ( );

  cout << "  BPATH plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void bpath_average ( int &seed, int m, int n, double u[], double umean[], 
  double &error )

//****************************************************************************80
//
//  Purpose:
//
//    BPATH_AVERAGE: displays the average of 1000 Brownian paths.
//
//  Discussion:
//
//    This routine computes M simulations of discretized Brownian 
//    motion W(t) over the time interval [0,1] using N time steps.
//    The user specifies a random number seed.  Different values of
//    the seed will result in a different set of realizations of the path.
//
//    Actually, we are interested in a function u(W(t)):
//
//      u(W(t)) = exp ( t + W(t)/2 )
//
//    The routine plots 5 of the simulations, as well as the average
//    of all the simulations.  
//
//    The plot of the average should be quite smooth.  Its expected
//    value is exp ( 9 * t / 8 ), and we compute the 'error', that is,
//    the difference between the averaged value and this expected
//    value.  This 'error' should decrease as the number of simulation
//    is increased.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Input, int M, the number of simulations to compute 
//    and average.  A typical value is 1000.
//
//    Input, int N, the number of steps.  A typical value
//    is 500.
//
//    Output, double U[M*(N+1)], the M paths.
//
//    Output, double UMEAN[N+1], the averaged path.
//
//    Output, double &ERROR, the maximum difference between the
//    averaged path and the exact expected value.
//
{
  double dt;
  double *dw;
  int i;
  int j;
  double *t;
  double tmax;
  double *w;

  tmax = 1.0;
  dt = tmax / ( double ) ( n );

  t = new double[n+1];
  for ( j = 0; j <= n; j++ )
  {
    t[j] = ( double ) ( j ) * tmax / ( double ) ( n );
  }

  w = new double[n+1];

  for ( i = 0; i < m; i++ )
  {
//
//  Define the increments dW.
//
    dw = r8vec_normal_01_new ( n, seed );

    for ( j = 0; j < n; j++ )
    {
      dw[j] = sqrt ( dt ) * dw[j];
    }
//
//  W is the sum of the previous increments.
//
    w[0] = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      w[j] = w[j-1] + dw[j-1];
    }

    for ( j = 0; j <= n; j++ )
    {
      u[i+j*m] = exp ( t[j] + 0.5 * w[j] );
    }
    delete [] dw;
  }
//
//  Average the M estimates of the path.
//
  for ( j = 0; j <= n; j++ )
  {
    umean[j] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      umean[j] = umean[j] + u[i+j*m];
    }
    umean[j] = umean[j] / ( double ) ( m );
  }

  error = 0.0;
  for ( j = 0; j <= n; j++ )
  {
    error = r8_max ( error, r8_abs ( umean[j] - exp ( 9.0 * t[j] / 8.0 ) ) );
  }

  delete [] t;
  delete [] w;

  return;
}
//****************************************************************************80

void bpath_average_gnuplot ( int m, int n, double u[], double umean[] )

//****************************************************************************80
//
//  Purpose:
//
//    BPATH_AVERAGE_GNUPLOT writes a GNUPLOT input file to plot BPATH_AVERAGE data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int M, the number of simulations.
//
//    Input, int N, the number of steps. 
//
//    Input, double U[M*(N+1)], the M paths.
//
//    Input, double UMEAN[N+1], the averaged path.
//
{
  string command_filename = "bpath_average_commands.txt";
  ofstream command;
  string data_filename = "bpath_average_data.txt";
  ofstream data;
  int i;
  int j;
  double t;
//
//  Create the data file.
//
  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i <= n; i++ )
  {
    t = ( double ) ( i ) / ( double ) ( n );
    data << "  " << t;
    for ( j = 0; j < 5; j++ )
    {
      data << "  " << u[j+i*m];
    }
    data << "  " << umean[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  BPATH_AVERAGE data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# bpath_average_commands.txt\n";
  command << "# created by sde::bpath_average_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < bpath_average_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'bpath_average.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel 'W(t)'\n";
  command << "set title 'Averaged Brownian paths'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot 'bpath_average_data.txt' using 1:2 title 'sample 1', \\\n";
  command << "     'bpath_average_data.txt' using 1:3 title 'sample 2', \\\n";
  command << "     'bpath_average_data.txt' using 1:4 title 'sample 3', \\\n";
  command << "     'bpath_average_data.txt' using 1:5 title 'sample 4', \\\n";
  command << "     'bpath_average_data.txt' using 1:6 title 'sample 5', \\\n";
  command << "     'bpath_average_data.txt' using 1:7 title 'average' lw 3\n";
  command << "quit\n";

  command.close ( );

  cout << "  BPATH_AVERAGE plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void chain ( int &seed, int n, double xem[], double vem[], double &diff )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN tests the stochastic Chain Rule.
//
//  Discussion:
//
//    This function solves a stochastic differential equation for 
//
//      V = sqrt(X) 
//
//    where X satisfies the stochastic differential equation:
// 
//      dX = ( alpha - X ) * dt + beta * sqrt(X) dW,
//      X(0) = Xzero,
//
//    with 
//
//      alpha = 2,
//      beta = 1,
//      Xzero = 1.
//
//    From the stochastic Chain Rule, the SDE for V is therefore:
//
//      dV = ( ( 4 * alpha - beta^2 ) / ( 8 * V ) - 1/2 V ) dt + 1/2 beta dW
//      V(0) = sqrt ( Xzero ).
//
//    Xem is the Euler-Maruyama solution for X. 
//
//    Vem is the Euler-Maruyama solution of the SDE for V from
//    the stochastic Chain Rule.
//
//    Hence, we compare sqrt(Xem) and Vem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int &SEED, a seed for the random number generator.
//
//    Input, int N, the number of time steps.
//
//    Output, double XEM[N+1], the computed value of X.
//
//    Output, double VEM[N+1], the computed value of V.
//
//    Output, double &DIFF, the maximum value of |sqrt(XEM)-V|.
//
{
  double alpha;
  double beta;
  double dt;
  double dt2;
  double *dw;
  int i;
  int j;
  double tmax;
//
//  Set problem parameters.
//
  alpha = 2.0;
  beta = 1.0;
//
//  Stepping parameters.
//  dt2 is the size of the Euler-Maruyama steps.
//
  tmax = 1.0;
  dt = tmax / ( double ) ( n );
  dt2 = dt;
//
//  Define the increments dW.
//
  dw = r8vec_normal_01_new ( n, seed );

  for ( i = 0; i < n; i++ )
  {
    dw[i] = sqrt ( dt ) * dw[i];
  }
//
//  Solve for X(t).
//
  xem[0] = 1.0;
  for ( j = 1; j <= n; j++ )
  {
    xem[j] = xem[j-1] + ( alpha - xem[j-1] ) * dt2 
                      + beta * sqrt ( xem[j-1] ) * dw[j-1];
  }
//
//  Solve for V(t).
//
  vem[0] = sqrt ( xem[0] );
  for ( j = 1; j <= n; j++ )
  {
    vem[j] = vem[j-1] 
      + ( ( 4.0 * alpha - beta * beta ) / ( 8.0 * vem[j-1] ) 
      - 0.5 * vem[j-1] ) * dt2 
      + 0.5 * beta * dw[j-1];
  }
//
//  Compare sqrt(X) and V.
//
  diff = 0.0;
  for ( i = 0; i <= n; i++ )
  {
    diff = r8_max ( diff, r8_abs ( sqrt ( xem[i] ) - vem[i] ) );
  }

  delete [] dw;

  return;
}
//****************************************************************************80

void chain_gnuplot ( int n, double x[], double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    CHAIN_GNUPLOT writes a GNUPLOT input file to plot CHAIN data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int N, the number of steps.
//
//    Input, double X[N+1], the value of X.
//
//    Input, double V[N+1], the value of V.
//
{
  string command_filename = "chain_commands.txt";
  ofstream command;
  string data_filename = "chain_data.txt";
  ofstream data;
  int i;
  double t;
//
//  Create the data file.
//

  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i <= n; i++ )
  {
    t = ( double ) ( i ) / ( double ) ( n );
    data << "  " << t
         << "  " << sqrt ( x[i] )
         << "  " << v[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  CHAIN data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# chain_commands.txt\n";
  command << "# created by sde::chain_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < chain_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'chain.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel 'Sqrt(X(t)) vs V(X(t))'\n";
  command << "set title 'V(X(t)) from X(t) and from Chain Rule'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot 'chain_data.txt' using 1:2 title 'Sqrt(X(t))', \\\n";
  command << "     'chain_data.txt' using 1:3 title 'V(X(t))'\n";
  command << "quit\n";

  command.close ( );

  cout << "  CHAIN plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void em ( int &seed, int n, double t[], double xtrue[], double t2[], 
  double xem[], double &error )

//****************************************************************************80
//
//  Purpose:
//
//    EM applies the Euler-Maruyama method to a linear SDE.
//
//  Discussion:
//
//    The SDE is 
//
//      dX = lambda * X dt + mu * X dW,   
//      X(0) = Xzero,
//
//    where 
//
//      lambda = 2,
//      mu = 1,
//      Xzero = 1.
//
//    The discretized Brownian path over [0,1] uses
//    a stepsize dt = 2^(-8).
//
//    The Euler-Maruyama method uses a larger timestep Dt = R*dt,
//    where R is an integer.  For an SDE of the form
//
//      dX = f(X(t)) dt + g(X(t)) dW(t)
//
//    it has the form
//
//      X(j) = X(j-1) + f(X(j-1)) * Dt + g(X(j-1)) * ( W(j*Dt) - W((j-1)*Dt) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Input, int N, the number of time steps.  A typical
//    value is 2^8.  N should be a multiple of 4.
//
//    Output, double T[N+1], the time values for the exact solution.
//
//    Output, double XTRUE[N+1], the exact solution.
//
//    Output, double T2[N/4+1], the time values for the 
//    Euler-Maruyama solution.
//
//    Output, double XEM[N/4+1], the Euler-Maruyama solution.
//
//    Output, double &ERROR, the value of | XEM(T) - XTRUE(T) |.
//
{
  double dt;
  double dt2;
  double *dw;
  double dw2;
  int i;
  int j;
  int l;
  double lambda;
  double mu;
  int r;
  double tmax;
  double *w;
  double xzero;
//
//  Set problem parameters.
//
  lambda = 2.0;
  mu = 1.0;
  xzero = 1.0;
//
//  Set stepping parameters.
//
  tmax = 1.0;
  dt = tmax / ( double ) ( n );
//
//  Define the increments dW.
//
  dw = r8vec_normal_01_new ( n, seed );

  for ( j = 0; j < n; j++ )
  {
    dw[j] = sqrt ( dt ) * dw[j];
  }
//
//  Sum the Brownian increments.
//
  w = new double[n+1];
  w[0] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    w[j] = w[j-1] + dw[j-1];
  }

  for ( j = 0; j <= n; j++)
  {
    t[j] = ( double ) ( j ) * tmax / ( double ) ( n );
  }
//
//  Compute the discretized Brownian path.
//
  for ( j = 0; j <= n; j++ )
  {
    xtrue[j] = xzero * exp ( ( lambda - 0.5 * mu * mu ) * ( t[j] + mu * w[j] ) );
  }
//
//  Set:
//  R, the multiplier for the EM step, 
//  Dt, the EM stepsize,
//  L, the number of EM steps (we need N to be a multiple of R!)
//
  r = 4;
  dt2 = ( double ) ( r ) * dt;
  l = n / r;

  for ( j = 0; j <= l; j++ )
  {
    t2[j] = ( double ) ( j ) * tmax / ( double ) ( l );
  }
//
//  Compute XEM.
//
  xem[0] = xzero;
  for ( j = 1; j <= l; j++ )
  {
    dw2 = 0.0;
    for ( i = r * ( j - 1 ); i < r * j; i++ )
    {
      dw2 = dw2 + dw[i];
    }
    xem[j] = xem[j-1] + dt2 * lambda * xem[j-1] + mu * xem[j-1] * dw2;
  }

  error = r8_abs ( xem[l] - xtrue[n] );

  delete [] dw;
  delete [] w;

  return;
}
//****************************************************************************80

void em_gnuplot ( int n, double t[], double xtrue[], double t2[], double xem[] )

//****************************************************************************80
//
//  Purpose:
//
//    EM_GNUPLOT writes a GNUPLOT input file to plot EM data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int N, the number of steps.
//
//    Input, double T[N+1], the time values for the exact solution.
//
//    Input, double XTRUE[N+1], the exact solution.
//
//    Input, double T2[N/4+1], the time values for the 
//    Euler-Maruyama solution.
//
//    Input, double XEM[N/4+1], the Euler-Maruyama solution.
//
{
  string command_filename = "em_commands.txt";
  ofstream command;
  string data1_filename = "em1_data.txt";
  string data2_filename = "em2_data.txt";
  ofstream data;
  int i;
//
//  Create data file #1.
//
  data.open ( data1_filename.c_str ( ) );

  for ( i = 0; i <= n; i++ )
  {
    data << "  " << t[i]
         << "  " << xtrue[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  EM data #1 stored in \"" << data1_filename << "\".\n";
//
//  Create data file #2.
//
  data.open ( data2_filename.c_str ( ) );

  for ( i = 0; i <= n / 4; i++ )
  {
    data << "  " << t2[i]
         << "  " << xem[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  EM data #2 stored in \"" << data2_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# em_commands.txt\n";
  command << "# created by sde::em_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < em_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'em.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel 'X(t)'\n";
  command << "set title 'Exact X(t) and Euler-Maruyama Estimate'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot 'em1_data.txt' using 1:2 title 'Exact X(t))', \\\n";
  command << "     'em2_data.txt' using 1:2 title 'EM X(t)'\n";
  command << "quit\n";

  command.close ( );

  cout << "  EM plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void emstrong ( int &seed, int m, int n, int p_max, double dtvals[], 
  double xerr[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMSTRONG tests the strong convergence of the EM method.
//
//  Discussion:
//
//    The SDE is 
//
//      dX = lambda * X dt + mu * X dW,   
//      X(0) = Xzero,
//
//    where 
//
//      lambda = 2,
//      mu = 1,
//      Xzero = 1.
//
//    The discretized Brownian path over [0,1] has dt = 2^(-9).
//
//    The Euler-Maruyama method uses 5 different timesteps: 
//      16*dt, 8*dt, 4*dt, 2*dt, dt.
//
//    We are interested in examining strong convergence at T=1,
//    that is
//
//      E | X_L - X(T) |.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Input, int M, the number of simulations to perform.
//    A typical value is M = 1000.
//
//    Input, int N, the number of time steps to take.
//    A typical value is N = 512.
//
//    Input, int P_MAX, the number of time step sizes to use.
//    A typical value is 5.
//
//    Output, double DTVALS[P_MAX], the time steps used.
//
//    Output, double XERR[P_MAX], the averaged absolute error in the
//    solution estimate at the final time.
//
{
  double *a;
  double dt;
  double dt2;
  double *dw;
  double e;
  int i;
  int j;
  int k;
  int l;
  double lambda;
  double mu;
  int p;
  double q;
  int r;
  double resid;
  double *rhs;
  int s;
  double *sol;
  double tmax;
  double *w;
  double winc;
  double xtemp;
  double xtrue;
  double xzero;
//
//  Set problem parameters.
//
  lambda = 2.0;
  mu = 1.0;
  xzero = 1.0;
//
//  Set stepping parameters.
//
  tmax = 1.0;
  dt = tmax / ( double ) ( n );

  for ( p = 0; p < p_max; p++ )
  {
    dtvals[p] = dt * pow ( 2.0, p );
  }
//
//  Sample over discrete Brownian paths.
//
  for ( p = 0; p < p_max; p++ )
  {
    xerr[p] = 0.0;
  }

  for ( s = 0; s < m; s++ )
  {
//
//  Define the increments dW.
//
    dw = r8vec_normal_01_new ( n, seed );

    for ( j = 0; j < n; j++ )
    {
      dw[j] = sqrt ( dt ) * dw[j];
    }
//
//  Sum the increments to get the Brownian path.
//
    w = new double[n+1];
    w[0] = 0.0;
    for ( j = 1; j <= n; j++ )
    {
      w[j] = w[j-1] + dw[j-1];
    }
//
//  Determine the true solution.
//
    xtrue = xzero * exp ( ( lambda - 0.5 * mu * mu ) + mu * w[n] );
//
//  Use the Euler-Maruyama method with 5 different time steps dt2 = r * dt
//  to estimate the solution value at time TMAX.
//
    for ( p = 0; p < p_max; p++ )
    {
      dt2 = dtvals[p];
      r = i4_power ( 2, p );
      l = n / r;
      xtemp = xzero;
      for ( j = 0; j < l; j++ )
      {
        winc = 0.0;
        for ( k = r * j; k < r * ( j + 1 ); k++ )
        {
          winc = winc + dw[k];
        }
        xtemp = xtemp + dt2 * lambda * xtemp + mu * xtemp * winc;
      }
      xerr[p] = xerr[p] + r8_abs ( xtemp - xtrue );
    }
    delete [] dw;
    delete [] w;
  }

  for ( p = 0; p < p_max; p++ )
  {
    xerr[p] = xerr[p] / ( double ) ( m );
  }
//
//  Least squares fit of error = c * dt^q.
//
  a = new double[p_max*2];
  rhs = new double[p_max];

  for ( i = 0; i < p_max; i++ )
  {
    a[i+0*p_max] = 1.0;
    a[i+1*p_max] = log ( dtvals[i] );
    rhs[i] = log ( xerr[i] );
  }
  sol = qr_solve ( p_max, 2, a, rhs );

  cout << "\n";
  cout << "EMSTRONG:\n";
  cout << "  Least squares solution to Error = c * dt ^ q\n";
  cout << "  (Expecting Q to be about 1/2.)\n";
  cout << "  Computed Q = " << sol[1] << "\n";

  resid = 0.0;
  for ( i = 0; i < p_max; i++ )
  {
    e = a[i+0*p_max] * sol[0] + a[i+1*p_max] * sol[1] - rhs[i];
    resid = resid + e * e;
  }
  resid = sqrt ( resid );
  cout << "  Residual is " << resid << "\n";

  delete [] a;
  delete [] rhs;
  delete [] sol;

  return;
}
//****************************************************************************80

void emstrong_gnuplot ( int p_max, double dtvals[], double xerr[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMSTRONG_GNUPLOT writes a GNUPLOT input file to plot EMSTRONG data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int P_MAX, the number of time step sizes to use.
//
//    Input, double DTVALS(P_MAX), the time steps used.
//
//    Input, double XERR(P_MAX), the averaged absolute error in the
//    solution estimate at the final time.
//
{
  string command_filename = "emstrong_commands.txt";
  ofstream command;
  string data_filename = "emstrong_data.txt";
  ofstream data;
  int i;
//
//  Create data file.
//
  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i < p_max; i++ )
  {
    data << "  " << dtvals[i]
         << "  " << xerr[i]
         << "  " << sqrt ( dtvals[i] ) << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  EMSTRONG data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# emstrong_commands.txt\n";
  command << "# created by sde::emstrong_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < emstrong_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'emstrong.png'\n";
  command << "set xlabel 'Log(dt)'\n";
  command << "set ylabel 'Log(Averaged Error at final T)'\n";
  command << "set logscale xy 10\n";
  command << "set title 'Euler-Maruyama Error as function of DT'\n";
  command << "set grid\n";
  command << "set style data linespoints\n";
  command << "plot 'emstrong_data.txt' using 1:2 title 'Error', \\\n";
  command << "     'emstrong_data.txt' using 1:3 title 'Slope = 1/2'\n";
  command << "quit\n";

  command.close ( );

  cout << "  EMSTRONG plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void emweak ( int &seed, int method, int m, int p_max, double dtvals[], 
  double xerr[] )

//****************************************************************************80
//
//  Purpose:
//
//    EMWEAK tests the weak convergence of the Euler-Maruyama method.
//
//  Discussion:
//
//    The SDE is 
//
//      dX = lambda * X dt + mu * X dW,   
//      X(0) = Xzero,
//
//    where 
//
//      lambda = 2,
//      mu = 1,
//      Xzero = 1.
//
//    The discretized Brownian path over [0,1] has dt = 2^(-9).
//
//    The Euler-Maruyama method will use 5 different timesteps:
//
//      2^(p-10),  p = 1,2,3,4,5.
//
//    We examine weak convergence at T=1:
//
//      | E (X_L) - E (X(T)) |.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 September 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546
//
//  Parameters:
//
//    Input, int &SEED, a seed for the random number generator.
//
//    Input, int METHOD.
//    0, use the standard Euler-Maruyama method;
//    1, use the weak Euler-Maruyama method.
//
//    Input, int M, the number of simulations to perform.
//    A typical value is M = 1000.
//
//    Input, int P_MAX, the number of time step sizes to use.
//    A typical value is 5.
//
//    Output, double DTVALS[P_MAX], the time steps used.
//
//    Output, double XERR[P_MAX], the averaged absolute error in the
//    solution estimate at the final time.
//
{
  double *a;
  double dt;
  double dt2;
  double e;
  int i;
  int j;
  int l;
  double lambda;
  double mu;
  int p;
  double q;
  int r;
  double resid;
  double *rhs;
  int s;
  double *sol;
  double tmax;
  double *winc;
  double *xem;
  double *xtemp;
  double xtrue;
  double xzero;
//
//  Problem parameters;
//
  lambda = 2.0;
  mu = 0.1;
  xzero = 1.0;
//
//  Stepping parameters.
//
  tmax = 1.0;
  for ( p = 0; p < p_max; p++ )
  {
    dtvals[p] = pow ( 2.0, p - 9 );
  }
//
//  Take various Euler timesteps.
//  For stepsize dt, we will need to take L Euler steps to reach time TMAX.
//
  xtemp = new double[m];
  xem = new double[p_max];

  for ( p = 0; p < p_max; p++ )
  {
    l = i4_power ( 2, 9 - p );
    dt = dtvals[p];

    for ( i = 0; i < m; i++ )
    {
      xtemp[i] = xzero;
    }

    for ( j = 0; j < l; j++ )
    {
      winc = r8vec_normal_01_new ( m, seed );
      if ( method == 0 )
      {
        for ( i = 0; i < m; i++ )
        {
          winc[i] = sqrt ( dt ) * winc[i];
        }
      }
      else
      {
        for ( i = 0; i < m; i++ )
        {
          winc[i] = sqrt ( dt ) * r8_sign ( winc[i] );
        }
      }

      for ( i = 0; i < m; i++ )
      {
        xtemp[i] = xtemp[i] + dt * lambda * xtemp[i] 
          + mu * xtemp[i] * winc[i];
      }
      delete [] winc;
    }
//
//  Average the M results for this stepsize.
//
    xem[p] = r8vec_mean ( m, xtemp );
  }
  delete [] xtemp;
//
//  Compute the error in the estimates for each stepsize.
//
  for ( p = 0; p < p_max; p++ )
  { 
    xerr[p] = r8_abs ( xem[p] - exp ( lambda ) );
  }
  delete [] xem;
//
//  Least squares fit of error = c * dt^q.
//
  a = new double[p_max*2];
  rhs = new double[p_max];

  for ( i = 0; i < p_max; i++ )
  {
    a[i+0*p_max] = 1.0;
    a[i+1*p_max] = log ( dtvals[i] );
    rhs[i] = log ( xerr[i] );
  }
  sol = qr_solve ( p_max, 2, a, rhs );

  cout << "\n";
  cout << "EMWEAK:\n";
  if ( method == 0 )
  {
    cout << "  Using standard Euler-Maruyama method.\n";
  }
  else
  {
    cout << "  Using weak Euler-Maruyama method.\n";
  }
  cout << "  Least squares solution to Error = c * dt ^ q\n";
  cout << "  (Expecting Q to be about 1.)\n";
  cout << "  Computed Q = " << sol[1] << "\n";

  resid = 0.0;
  for ( i = 0; i < p_max; i++ )
  {
    e = a[i+0*p_max] * sol[0] + a[i+1*p_max] * sol[1] - rhs[i];
    resid = resid + e * e;
  }
  resid = sqrt ( resid );
  cout << "  Residual is " << resid << "\n";

  delete [] a;
  delete [] rhs;
  delete [] sol;

  return;
}
//****************************************************************************80

void emweak_gnuplot ( int p_max, double dtvals[], double xerr[], int method )

//****************************************************************************80
//
//  Purpose:
//
//    EMWEAK_GNUPLOT writes a GNUPLOT input file to plot EMWEAK data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int P_MAX, the number of time step sizes to use.
//
//    Input, double DTVALS[P_MAX], the time steps used.
//
//    Input, double XERR[P_MAX], the averaged absolute error in the
//    solution estimate at the final time.
//
//    Input, int METHOD.
//    0, use the standard Euler-Maruyama method;
//    1, use the weak Euler-Maruyama method.
//
{
  string command_filename;
  ofstream command;
  string data_filename;
  ofstream data;
  int i;
//
//  Create data file.
//
  if ( method == 0 )
  {
    data_filename = "emweak0_data.txt";
  }
  else 
  {
    data_filename = "emweak1_data.txt";
  }

  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i < p_max; i++ )
  {
    data << "  " << dtvals[i]
         << "  " << xerr[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  EMWEAK data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  if ( method == 0 )
  {
    command_filename = "emweak0_commands.txt";
  }
  else
  {
    command_filename = "emweak1_commands.txt";
  }

  command.open ( command_filename.c_str ( ) );

  command << "# " << command_filename << "\n";
  command << "# created by sde::emweak_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  if ( method == 0 )
  {
    command << "set output 'emweak0.png'\n";
  }
  else
  {
    command << "set output 'emweak1.png'\n";
  }
  command << "set xlabel 'Log(dt)'\n";
  command << "set ylabel 'Log(Averaged Error at final T)'\n";
  command << "set logscale xy 10\n";
  if ( method == 0 )
  {
    command << "set title 'Standard Euler-Maruyama Error as function of DT'\n";
  }
  else
  {
    command << "set title 'Weak Euler-Maruyama Error as function of DT'\n";
  }
  command << "set grid\n";
  command << "set style data linespoints\n";
  if ( method == 0 )
  {
    command << "plot 'emweak0_data.txt' using 1:2 title 'Error', \\\n";
    command << "     'emweak0_data.txt' using 1:1 title 'Slope = 1'\n";
  }
  else
  {
    command << "plot 'emweak1_data.txt' using 1:2 title 'Error', \\\n";
    command << "     'emweak1_data.txt' using 1:1 title 'Slope = 1'\n";
  }

  command << "quit\n";

  command.close ( );

  cout << "  EMWEAK plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void filename_inc ( string *filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILENAME_INC increments a partially numeric file name.
//
//  Discussion:
//
//    It is assumed that the digits in the name, whether scattered or
//    connected, represent a number that is to be increased by 1 on
//    each call.  If this number is all 9's on input, the output number
//    is all 0's.  Non-numeric letters of the name are unaffected.
//
//    If the name is empty, then the routine stops.
//
//    If the name contains no digits, the empty string is returned.
//
//  Example:
//
//      Input            Output
//      -----            ------
//      "a7to11.txt"     "a7to12.txt"  (typical case.  Last digit incremented)
//      "a7to99.txt"     "a8to00.txt"  (last digit incremented, with carry.)
//      "a9to99.txt"     "a0to00.txt"  (wrap around)
//      "cat.txt"        " "           (no digits to increment)
//      " "              STOP!         (error)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, string *FILENAME, the filename to be incremented.
//
{
  char c;
  int change;
  int i;
  int lens;

  lens = (*filename).length ( );

  if ( lens <= 0 )
  {
    cerr << "\n";
    cerr << "FILENAME_INC - Fatal error!\n";
    cerr << "  The input string is empty.\n";
    exit ( 1 );
  }

  change = 0;

  for ( i = lens - 1; 0 <= i; i-- )
  {
    c = (*filename)[i];

    if ( '0' <= c && c <= '9' )
    {
      change = change + 1;

      if ( c == '9' )
      {
        c = '0';
        (*filename)[i] = c;
      }
      else
      {
        c = c + 1;
        (*filename)[i] = c;
        return;
      }
    }
  }
//
//  No digits were found.  Return blank.
//
  if ( change == 0 )
  {
    for ( i = lens - 1; 0 <= i; i-- )
    {
      (*filename)[i] = ' ';
    }
  }

  return;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cerr << "\n";
      cerr << "I4_POWER - Fatal error!\n";
      cerr << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void milstrong ( int &seed, int p_max, double dtvals[], double xerr[] )

//****************************************************************************80
//
//  Purpose:
//
//    MILSTRONG tests the strong convergence of the Milstein method.
//
//  Discussion:
//
//    This function solves the stochastic differential equation
//
//      dX = sigma * X * ( k - X ) dt + beta * X dW,  
//      X(0) = Xzero,
//
//    where 
//
//       sigma = 2, 
//       k = 1, 
//       beta = 1,
//       Xzero = 0.5.
//
//    The discretized Brownian path over [0,1] has dt = 2^(-11).
//
//    The Milstein method uses timesteps 128*dt, 64*dt, 32*dt, 16*dt 
//    (also dt for reference).
//
//    We examine strong convergence at T=1:  
//
//      E | X_L - X(T) |.
//
//    The code is vectorized: all paths computed simultaneously.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original MATLAB version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Input, int P_MAX, the number of time step sizes to use.
//    A typical value is 4.
//
//    Output, double DTVALS[P_MAX], the time steps used.
//
//    Output, double XERR[P_MAX], the averaged absolute error in the
//    solution estimate at the final time.
//
{
  double *a;
  double beta;
  double dt;
  double dtp;
  double *dw;
  double e;
  int i;
  int i2;
  int j;
  double k;
  int l;
  int m;
  int n;
  int p;
  int r;
  double resid;
  double *rhs;
  double sigma;
  double *sol;
  double tmax;
  double winc;
  double *xref;
  double *xtemp;
  double xzero;
//
//  Set problem parameters.
//
  sigma = 2.0;
  k = 1.0;
  beta = 0.25;
  xzero = 0.5;
//
//  Set stepping parameters.
//
  tmax = 1.0;
  n = i4_power ( 2, 11 );
  dt = tmax / ( double ) ( n );
//
//  Number of paths sampled.
//
  m = 500;
//
//  Define the increments dW.
//
  dw = r8mat_normal_01_new ( m, n, seed );
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      dw[i+j*m] = sqrt ( dt ) * dw[i+j*m];
    }
  }
//
//  Estimate the reference solution at time T M times.
//
  xref = new double[m];

  for ( i = 0; i < m; i++ )
  {
    xref[i] = xzero;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      xref[i] = xref[i] 
        + dt * sigma * xref[i] * ( k - xref[i] ) 
        + beta * xref[i] * dw[i+j*m] 
        + 0.5 * beta * beta * xref[i] * ( dw[i+j*m] * dw[i+j*m] - dt );
    }
  }
//
//  Now compute M Milstein approximations at each of 4 timesteps,
//  and record the average errors.
//
  for ( p = 0; p < p_max; p++ )
  {
    dtvals[p] = dt * 8.0 * pow ( 2.0, p + 1 );
  }
  for ( p = 0; p < p_max; p++ )
  {
    xerr[p] = 0.0;
  }
  xtemp = new double[m];
  for ( p = 0; p < p_max; p++ )
  {
    r = 8 * i4_power ( 2, p + 1 );
    dtp = dtvals[p];
    l = n / r;
    for ( i = 0; i < m; i++ )
    {
      xtemp[i] = xzero;
    }
    for ( j = 0; j < l; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        winc = 0.0;
        for ( i2 = r * j; i2 < r * ( j + 1 ); i2++ )
        {
          winc = winc + dw[i+i2*m];
        }
        xtemp[i] = xtemp[i] 
          + dtp * sigma * xtemp[i] * ( k - xtemp[i] ) 
          + beta * xtemp[i] * winc
          + 0.5 * beta * beta * xtemp[i] * ( winc * winc - dtp );
      }
    }

    xerr[p] = 0.0;
    for ( i = 0; i < m; i++ )
    {
      xerr[p] = xerr[p] + r8_abs ( xtemp[i] - xref[i] );
    }
    xerr[p] = xerr[p] / ( double ) ( m );
  }
//
//  Least squares fit of error = C * dt^q
//
  a = new double[p_max*2];
  rhs = new double[p_max];
  for ( p = 0; p < p_max; p++ )
  {
    a[p+0*p_max] = 1.0;
    a[p+1*p_max] = log ( dtvals[p] );
    rhs[p] = log ( xerr[p] );
  }
  sol = qr_solve ( p_max, 2, a, rhs );

  cout << "\n";
  cout << "MILSTEIN:\n";
  cout << "  Least squares solution to Error = c * dt ^ q\n";
  cout << "  Expecting Q to be about 1.\n";
  cout << "  Computed Q = " << sol[1] << "\n";

  resid = 0.0;
  for ( i = 0; i < p_max; i++ )
  {
    e = a[i+0*p_max] * sol[0] + a[i+1*p_max] * sol[1] - rhs[i];
    resid = resid + e * e;
  }
  resid = sqrt ( resid );
  cout << "  Residual is " << resid << "\n";

  delete [] a;
  delete [] dw;
  delete [] rhs;
  delete [] sol;
  delete [] xref;
  delete [] xtemp;

  return;
}
//****************************************************************************80

void milstrong_gnuplot ( int p_max, double dtvals[], double xerr[] )

//****************************************************************************80
//
//  Purpose:
//
//    MILSTRONG_GNUPLOT writes a GNUPLOT input file to plot MILSTRONG data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int P_MAX, the number of time step sizes to use.
//
//    Input, double DTVALS(P_MAX), the time steps used.
//
//    Input, double XERR(P_MAX), the averaged absolute error in the
//    solution estimate at the final time.
//
{
  string command_filename = "milstrong_commands.txt";
  ofstream command;
  string data_filename = "milstrong_data.txt";
  ofstream data;
  int i;
//
//  Create data file.
//
  data.open ( data_filename.c_str ( ) );

  for ( i = 0; i < p_max; i++ )
  {
    data << "  " << dtvals[i]
         << "  " << xerr[i] << "\n";
  }
  data.close ( );

  cout << "\n";
  cout << "  MILSTRONG data stored in \"" << data_filename << "\".\n";
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# milstrong_commands.txt\n";
  command << "# created by sde::milstrong_gnuplot.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < milstrong_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'milstrong.png'\n";
  command << "set xlabel 'Log(dt)'\n";
  command << "set ylabel 'Log(Averaged Error at final T)'\n";
  command << "set logscale xy 10\n";
  command << "set title 'Milstein Error as function of DT'\n";
  command << "set grid\n";
  command << "set style data linespoints\n";
  command << "plot 'milstrong_data.txt' using 1:2 title 'Error', \\\n";
  command << "     'milstrong_data.txt' using 1:1 title 'Slope = 1'\n";
  command << "quit\n";

  command.close ( );

  cout << "  MILSTRONG plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

double r8_normal_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = r8_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
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

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
//****************************************************************************80

double r8vec_mean ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEAN returns the mean of an R8VEC.
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
//    02 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double X[N], the vector whose mean is desired.
//
//    Output, double R8VEC_MEAN, the mean, or average, of the vector entries.
//
{
  int i;
  double mean;

  mean = 0.0;
  for ( i = 0; i < n; i++ )
  {
    mean = mean + x[i];
  }

  mean = mean / ( double ) n;

  return mean;
}
//****************************************************************************80

void r8vec_normal_01 ( int n, int &seed, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
//
//  Discussion:
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
//    09 April 2012
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
//    Output, double X[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, double R(N+1), is used to store some uniform random values.
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
# define R8_PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
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
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
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

    x[x_hi-1] = sqrt ( - 2.0 * log ( r[0] ) ) * cos ( 2.0 * R8_PI * r[1] );
    y =         sqrt ( - 2.0 * log ( r[0] ) ) * sin ( 2.0 * R8_PI * r[1] );

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

    for ( i = 0; i <= 2 * m - 2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
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

    for ( i = 0; i <= 2 * m - 4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( - 2.0 * log ( r[i] ) ) * cos ( 2.0 * R8_PI * r[i+1] );
    y           = sqrt ( - 2.0 * log ( r[i] ) ) * sin ( 2.0 * R8_PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return;
# undef R8_PI
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

void r8vec_uniform_01 ( int n, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

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

void stab_asymptotic ( int &seed, int n, int p_max )

//****************************************************************************80
//
//  Purpose:
//
//    STAB_ASYMPTOTIC examines asymptotic stability.
//
//  Discussion:
//
//    The function tests the asymptotic stability
//    of the Euler-Maruyama method applied to a stochastic differential
//    equation (SDE).
//
//    The SDE is
//
//      dX = lambda*X dt + mu*X dW,
//      X(0) = Xzero,
//
//    where 
//
//      lambda is a constant,
//      mu is a constant,
//      Xzero = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input/output, int &SEED, a seed for the random 
//    number generator.
//
//    Input, int N, the number of time steps for the
//    first solution.
//
//    Input, int P_MAX, the number of time step sizes.
//
{
  string command_filename = "stab_asymptotic_commands.txt";
  ofstream command;
  string data_filename;
  string data_filename0 = "stab_asymptotic0_data.txt";
  ofstream data;
  double dt;
  double *dtvals;
  int i;
  int j;
  double lambda;
  double mu;
  int nval;
  int p;
  double t;
  double test;
  double tmax;
  double *u;
  double winc;
  double *xemabs;
  double xmin;
  double xtemp;
  double xzero;

  cout << "\n";
  cout << "STAB_ASYMPTOTIC:\n";
  cout << "  Investigate asymptotic stability of Euler-Maruyama\n";
  cout << "  solution with stepsize DT and MU.\n";
  cout << "\n";
  cout << "  SDE is asymptotically stable if\n";
  cout << "    Real ( lambda - 1/2 mu^2 ) < 0.\n";
  cout << "\n";
  cout << "  EM with DT is asymptotically stable if\n";
  cout << "    E log ( | 1 + lambda dt - sqrt(dt) mu n(0,1) | ) < 0.\n";
  cout << "  where n(0,1) is a normal random value.\n";
//
//  Problem parameters.
//
  lambda = 0.5;
  mu = sqrt ( 6.0 );
  xzero = 1.0;
//
//  Test the SDE.
//
  cout << "\n";
  cout << "  Lambda = " << lambda << "\n";
  cout << "  Mu =     " << mu << "\n";
  test = lambda - 0.5 * mu * mu;
  cout << "  SDE asymptotic stability test = " << test << "\n";
//
//  Step parameters.
//
  tmax = 500.0;
//
//  For each stepsize, compute the Euler-Maruyama solution.
//
  data_filename = data_filename0;
  dtvals = new double[p_max];

  for ( p = 0; p < p_max; p++ )
  {
    nval = n * i4_power ( 2, p );
    dt = tmax / ( double ) ( nval );
    dtvals[p] = dt;
//
//  Test the EM for this DT.
//
    cout << "\n";
    cout << "  dt = " << dt << "\n";
    u = r8vec_normal_01_new ( 1000, seed );
    for ( i = 0; i < 1000; i++ )
    {
      u[i] = log ( r8_abs ( 1.0 + lambda * dt - sqrt ( dt ) * mu * u[i] ) );
    }
    test = r8vec_mean ( 1000, u );
    cout << "  EM asymptotic test = " << test << "\n";
    delete [] u;

    xtemp = xzero;
    xemabs = new double[nval+1];
    xemabs[0] = xtemp;

    for ( j = 1; j <= nval; j++ )
    {
      winc = sqrt ( dt ) * r8_normal_01 ( seed );
      xtemp = xtemp + dt * lambda * xtemp + mu * xtemp * winc;
      xemabs[j] = r8_abs ( xtemp );
    }
//
//  Write this data to a file.
//
    filename_inc ( &data_filename );

    data.open ( data_filename.c_str ( ) );
//
//  We have to impose a tiny lower bound on the values because we
//  will end up plotting their logs.
//
    xmin = exp ( -200.0 );
    for ( i = 0; i <= nval; i++ )
    {
      t = tmax * ( double ) ( i ) / ( double ) ( nval );
      data << "  " << t
           << "  " << r8_max ( xemabs[i], xmin ) << "\n";
    }
    data.close ( );

    cout << "\n";
    cout << "  Data for DT = " << dt << " stored in \"" << data_filename << "\"\n";

    delete [] xemabs;
  }
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# stab_asymptotic_commands.txt\n";
  command << "# created by sde::stab_asymptotic.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < stab_asymptotic_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'stab_asymptotic.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel '|X(t)|'\n";
  command << "set title 'Absolute value of EM Solution'\n";
  command << "set grid\n";
  command << "set logscale y 10\n";
  command << "set style data lines\n";

  data_filename = data_filename0;

  filename_inc ( &data_filename );
  command << "plot '" << data_filename << "' using 1:2, \\\n";

  for ( p = 1; p < p_max - 1; p++ )
  {
    filename_inc ( &data_filename );
    command << "     '" << data_filename << "' using 1:2, \\\n";
  }
  filename_inc ( &data_filename );
  command << "     '" << data_filename << "' using 1:2\n";

  command << "quit\n";

  command.close ( );

  cout << "  STAB_ASYMPTOTIC plot stored in \"" << command_filename << "\".\n";

  delete [] dtvals;

  return;
}
//****************************************************************************80

void stab_meansquare ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    STAB_MEANSQUARE examines mean-square stability.
//
//  Discussion:
//
//    The function tests the mean-square stability
//    of the Euler-Maruyama method applied to a stochastic differential
//    equation (SDE).
//
//    The SDE is
//
//      dX = lambda*X dt + mu*X dW,
//      X(0) = Xzero,
//
//    where 
//
//      lambda is a constant,
//      mu is a constant,
//      Xzero = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int &SEED, a seed for the random number generator.
//    In the reference, this value is set to 100.
//
{
  string command_filename = "stab_meansquare_commands.txt";
  ofstream command;
  string data_filename0 = "stab_meansquare0_data.txt";
  string data_filename;
  ofstream data;
  double dt;
  int i;
  int j;
  int k;
  double lambda;
  int m;
  double mu;
  int n;
  double t;
  double test;
  double tmax;
  double *winc;
  double *xms;
  double *xtemp;
  double xzero;

  cout << "\n";
  cout << "STAB_MEANSQUARE:\n";
  cout << "  Investigate mean square stability of Euler-Maruyama\n";
  cout << "  solution with stepsize DT and MU.\n";
  cout << "\n";
  cout << "  SDE is mean square stable if\n";
  cout << "    Real ( lambda + 1/2 |mu|^2 ) < 0.\n";
  cout << "\n";
  cout << "  EM with DT is mean square stable if\n";
  cout << "    |1+dt^2| + dt * |mu|^2 - 1.0 < 0.\n";
//
//  Set problem parameters.
//
  tmax = 20.0;
  m = 50000;
  xzero = 1.0 ;
//
//  Problem parameters.
//
  lambda = -3.0;
  mu = sqrt ( 3.0 );
//
//  Test the SDE.
//
  cout << "\n";
  cout << "  Lambda = " << lambda << "\n";
  cout << "  Mu =     " << mu << "\n";
  test = lambda + 0.5 * mu * mu;
  cout << "  SDE mean square stability test = " << test << "\n";
//
//  XMS is the mean square estimate of M paths.
//
  data_filename = data_filename0;

  for ( k = 0; k < 3; k++ )
  {
    dt = pow ( 2.0, - k );
    n = 20 * i4_power ( 2, k );
//
//  Test the EM for this DT.
//
    cout << "\n";           
    cout << "  dt = " << dt << "\n";
    test = pow ( 1.0 + dt * lambda, 2 ) + dt * mu * mu - 1.0;
    cout << "  EM mean square stability test = " << test << "\n";

    xms = new double[n+1];
    xtemp = new double[m];
    for ( i = 0; i < m; i++ )
    {
      xtemp[i] = xzero;
    }
    xms[0] = xzero;

    for ( j = 0; j <= n; j++ )
    {
      winc = r8vec_normal_01_new ( m, seed );
      for ( i = 0; i < m; i++ )
      {
        winc[i] = sqrt ( dt ) * winc[i];
      }
      for ( i = 0; i < m; i++ )
      {
        xtemp[i] = xtemp[i] 
          + dt * lambda * xtemp[i] 
          + mu * xtemp[i] * winc[i];
      }
      xms[j] = 0.0;
      for ( i = 0; i < m; i++ )
      {
        xms[j] = xms[j] + xtemp[i] * xtemp[i];
      }
      xms[j] = xms[j] / ( double ) ( m );
      delete [] winc;
    }
//
//  Write this data to a file.
//
    filename_inc ( &data_filename );

    data.open ( data_filename.c_str ( ) );

    for ( j = 0; j <= n; j++ )
    {
      t = tmax * ( double ) ( j ) / ( double ) ( n );
      data << "  " << t
           << "  " << xms[j] << "\n";
    }
    data.close ( );

    cout << "\n";
    cout << "  Data for DT = " << dt << " stored in \"" << data_filename << "\".\n";

    delete [] xtemp;
    delete [] xms;
  }
//
//  Create the command file.
//
  command.open ( command_filename.c_str ( ) );

  command << "# stab_meansquare_commands.txt\n";
  command << "# created by sde::stab_meansquare.\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < stab_meansquare_commands.txt\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output 'stab_meansquare.png'\n";
  command << "set xlabel 't'\n";
  command << "set ylabel 'E|X^2(t)|'\n";
  command << "set title 'Mean Square of EM Solution'\n";
  command << "set grid\n";
  command << "set logscale y 10\n";
  command << "set style data lines\n";

  data_filename = data_filename0;

  filename_inc ( &data_filename );
  command << "plot '" << data_filename << "' using 1:2, \\\n";

  for ( k = 1; k <= 1; k++ )
  {
    filename_inc ( &data_filename );
    command << "     '" << data_filename << "' using 1:2, \\\n";
  }
  filename_inc ( &data_filename );
  command << "     '" << data_filename << "' using 1:2\n";

  command << "quit\n";

  command.close ( );

  cout << "  STAB_MEANSQUARE plot commands stored in \"" << command_filename << "\".\n";

  return;
}
//****************************************************************************80

void stochastic_integral_ito ( int n, int &seed, double &estimate, 
  double &exact, double &error )

//****************************************************************************80
//
//  Purpose:
//
//    STOCHASTIC_INTEGRAL_ITO approximates the Ito integral of W(t) dW.
//
//  Discussion:
//
//    This function estimates the Ito integral of W(t) dW over 
//    the interval [0,1].
//
//    The estimates is made by taking N steps.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int N, the number of steps to take.
//
//    Input, int &SEED, a seed for the random number generator.
//
//    Output, double &ESTIMATE, the estimate of the integral.
//
//    Output, double &EXACT, the exact value of the integral.
//
//    Output, double &ERROR, the error in the integral estimate.
//
{
  double dt;
  double *dw;
  int j;
  double tmax;
  double *w;
//
//  Set step parameters.
//
  tmax = 1.0;
  dt = tmax / ( double ) ( n );
//
//  Define the increments dW.
//
  dw = r8vec_normal_01_new ( n, seed );
  for ( j = 0; j < n; j++ )
  {
    dw[j] = sqrt ( dt ) * dw[j];
  }
//
//  Sum the increments to get the Brownian path.
//
  w = new double[n+1];
  w[0] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    w[j] = w[j-1] + dw[j-1];
  }
//
//  Approximate the Ito integral.
//
  estimate = r8vec_dot_product ( n, w, dw );
//
//  Compare with the exact solution.
//
  exact = 0.5 * ( w[n] * w[n] - tmax );
  error = r8_abs ( estimate - exact );

  delete [] dw;
  delete [] w;

  return;
}
//****************************************************************************80

void stochastic_integral_strat ( int n, int &seed, double &estimate, 
  double &exact, double &error )

//****************************************************************************80
//
//  Purpose:
//
//    STOCHASTIC_INTEGRAL_STRAT approximates the Stratonovich integral of W(t) dW.
//
//  Discussion:
//
//    This function estimates the Stratonovich integral of W(t) dW over 
//    the interval [0,1].
//
//    The estimates is made by taking N steps.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2012
//
//  Author:
//
//    Original Matlab version by Desmond Higham.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Desmond Higham,
//    An Algorithmic Introduction to Numerical Simulation of
//    Stochastic Differential Equations,
//    SIAM Review,
//    Volume 43, Number 3, September 2001, pages 525-546.
//
//  Parameters:
//
//    Input, int N, the number of steps to take.
//
//    Input, int &SEED, a seed for the random number generator.
//
//    Output, double &ESTIMATE, the estimate of the integral.
//
//    Output, double &EXACT, the exact value of the integral.
//
//    Output, double &ERROR, the error in the integral estimate.
//
{
  double dt;
  double *dw;
  int j;
  double tmax;
  double *u;
  double *v;
  double *w;
//
//  Set step parameters.
//
  tmax = 1.0;
  dt = tmax / ( double ) ( n );
//
//  Define the increments dW.
//
  dw = r8vec_normal_01_new ( n, seed );

  for ( j = 0; j < n; j++ )
  {
    dw[j] = sqrt ( dt ) * dw[j];
  }
//
//  Sum the increments to get the Brownian path.
//
  w = new double[n+1];
  w[0] = 0.0;
  for ( j = 1; j <= n; j++ )
  {
    w[j] = w[j-1] + dw[j-1];
  }
//
//  Approximate the Stratonovich integral.
//
  u = r8vec_normal_01_new ( n, seed );

  v = new double[n];
  for ( j = 0; j < n; j++ )
  {
    v[j] = 0.5 * ( w[j] + w[j+1] ) + 0.5 * sqrt ( dt ) * u[j];
  }

  estimate = r8vec_dot_product ( n, v, dw );
//
//  Compare with the exact solution.
//
  exact = 0.5 * w[n-1] * w[n-1];
  error = r8_abs ( estimate - exact );

  delete [] dw;
  delete [] u;
  delete [] v;
  delete [] w;

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
