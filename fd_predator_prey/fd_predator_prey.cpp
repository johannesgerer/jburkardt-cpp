# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>

using namespace std;

int main ( int argc, char *argv[] );
void r8mat_write ( string output_filename, int m, int n, double table[] );
int s_len_trim ( char *s );
int s_to_i4 ( char *s, int *last, bool *error );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    FD_PREDATOR_PREY solves a pair of predator-prey ODE's.
//
//  Discussion:
//
//    The physical system under consideration is a pair of animal populations.
//
//    The PREY reproduce rapidly; for each animal alive at the beginning of the
//    year, two more will be born by the end of the year.  The prey do not have
//    a natural death rate; instead, they only die by being eaten by the predator.
//    Every prey animal has 1 chance in 1000 of being eaten in a given year by
//    a given predator.
//
//    The PREDATORS only die of starvation, but this happens very quickly.
//    If unfed, a predator will tend to starve in about 1/10 of a year.
//    On the other hand, the predator reproduction rate is dependent on
//    eating prey, and the chances of this depend on the number of available prey.
//
//    The resulting differential equations can be written:
//
//      PREY(0) = 5000
//      PRED(0) =  100
//
//      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
//      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
//
//    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
//
//    The pair of ordinary differential equations that result have an interesting
//    behavior.  For certain choices of the interaction coefficients (such as
//    those given here), the populations of predator and prey will tend to
//    a periodic oscillation.  The two populations will be out of phase; the number
//    of prey will rise, then after a delay, the predators will rise as the prey
//    begins to fall, causing the predator population to crash again.
//
//    In this program, the pair of ODE's is solved with a simple finite difference
//    approximation using a fixed step size.  In general, this is NOT an efficient
//    or reliable way of solving differential equations.  However, this program is
//    intended to illustrate the ideas of finite difference approximation.
//
//    In particular, if we choose a fixed time step size DT, then a derivative
//    such as dPREY/dT is approximated by:
//
//      d PREY / dT = approximately ( PREY(T+DT) - PREY(T) ) / DT
//
//    which means that the first differential equation can be written as
//
//      PREY(T+DT) = PREY(T) + DT * ( 2 * PREY(T) - 0.001 * PREY(T) * PRED(T) ).
//
//    We can rewrite the equation for PRED as well.  Then, since we know the
//    values of PRED and PREY at time 0, we can use these finite difference
//    equations to estimate the values of PRED and PREY at time DT.  These values
//    can be used to get estimates at time 2*DT, and so on.  To get from time
//    T_START = 0 to time T_STOP = 5, we simply take STEP_NUM steps each of size
//    DT = ( T_STOP - T_START ) / STEP_NUM.
//
//    Because finite differences are only an approximation to derivatives, this
//    process only produces estimates of the solution.  And these estimates tend
//    to become more inaccurate for large values of time.  Usually, we can reduce
//    this error by decreasing DT and taking more, smaller time steps.
//
//    In this example, for instance, taking just 100 steps gives nonsensical
//    answers.  Using STEP_NUM = 1000 gives an approximate solution that seems
//    to have the right kind of oscillatory behavior, except that the amplitude
//    of the waves increases with each repetition.  Using 10000 steps, the
//    approximation begins to become accurate enough that we can see that the
//    waves seem to have a fixed period and amplitude.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 May 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Lindfield, John Penny,
//    Numerical Methods Using MATLAB,
//    Second Edition,
//    Prentice Hall, 1999,
//    ISBN: 0-13-012641-1,
//    LC: QA297.P45.
//
//  Parameters:
//
//    Input, int STEP_NUM, the number of steps.
//
{
  double dt;
  int i;
  bool ierror;
  int length;
  string output_filename;
  int step_num;
  double t_start;
  double t_stop;
  double *trf;

  timestamp ( );

  cout << "\n";
  cout << "FD_PREDATOR_PREY\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  A finite difference approximate solution of a pair\n";
  cout << "  of ordinary differential equations for a population\n";
  cout << "  of predators and prey.\n";
  cout << "\n";
  cout << "  The exact solution shows wave behavior, with a fixed\n";
  cout << "  period and amplitude.  The finite difference approximation\n";
  cout << "  can provide a good estimate for this behavior if the stepsize\n";
  cout << "  DT is small enough.\n";
//
//  STEP_NUM is an input argument or else read from the user interactively.
//
  if ( 1 < argc )
  {
    step_num = atoi ( argv[1] );
  }
  else
  {
    cout << "\n";
    cout << "FD_PREDATOR_PREY:\n";
    cout << "  Please enter the number of time steps:\n";

    cin >> step_num;
  }

  t_start = 0.0;
  t_stop =  5.0;
  dt = ( t_stop - t_start ) / ( double ) ( step_num );
//
//  TRF(1:3,1:STEP_NUM+1) contains TIME, PREY, and PREDATOR values for each step.
//
  trf = new double[3*(step_num+1)];

  trf[0+(0)*3] = t_start;
  trf[1+(0)*3] = 5000.0;
  trf[2+(0)*3] =  100.0;

  for ( i = 0; i < step_num; i++ )
  {
    trf[0+(i+1)*3] = trf[0+(i)*3] + dt;

    trf[1+(i+1)*3] = trf[1+(i)*3] + dt 
      * (    2.0 * trf[1+(i)*3] - 0.001 * trf[1+(i)*3] * trf[2+(i)*3] );

    trf[2+(i+1)*3] = trf[2+(i)*3] + dt 
      * ( - 10.0 * trf[2+(i)*3] + 0.002 * trf[1+(i)*3] * trf[2+(i)*3] );
  }
  output_filename = "trf_";
  output_filename = output_filename + argv[1];
  output_filename = output_filename + ".txt";

  r8mat_write ( output_filename, 3, step_num + 1, trf );

  cout << "\n";
  cout << "  Initial time    = " << t_start << "\n";
  cout << "  Final time      = " << t_stop << "\n";
  cout << "  Initial prey    = " << trf[1+0*3] << "\n";
  cout << "  Initial pred    = " << trf[2+0*3] << "\n";
  cout << "  Number of steps = " << step_num << "\n";
  cout << "  Solution data written to \"" << output_filename << "\".\n";

  delete [] trf;
//
//  Terminate.
//
  cout << "\n";
  cout << "FD_PREDATOR_PREY\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
