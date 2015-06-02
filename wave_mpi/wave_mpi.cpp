# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "mpi.h"

int main ( int argc, char *argv[] );
double *update ( int id, int p, int n_global, int n_local, int nsteps, 
  double dt );
void collect ( int id, int p, int n_global, int n_local, int nsteps, 
  double dt, double u_local[] );
double dudt ( double x, double t );
double exact ( double x, double t );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    WAVE_MPI solves the wave equation in parallel using MPI.
//
//  Discussion:
//
//    Discretize the equation for u(x,t):
//      d^2 u/dt^2 - c^2 * d^2 u/dx^2 = 0  for 0 < x < 1, 0 < t
//    with boundary conditions:
//      u(0,t) = u0(t) = sin ( 2 * pi * ( 0 - c * t ) )
//      u(1,t) = u1(t) = sin ( 2 * pi * ( 1 - c * t ) )
//    and initial conditions:
//         u(x,0) = g(x,t=0) =                sin ( 2 * pi * ( x - c * t ) )
//      dudt(x,0) = h(x,t=0) = - 2 * pi * c * cos ( 2 * pi * ( x - c * t ) ) 
//
//    by:
//
//      alpha = c * dt / dx.
//
//      U(x,t+dt) = 2 U(x,t) - U(x,t-dt) 
//        + alpha^2 ( U(x-dx,t) - 2 U(x,t) + U(x+dx,t) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Geoffrey Fox, Mark Johnson, Gregory Lyzenga, Steve Otto, John Salmon, 
//    David Walker,
//    Solving problems on concurrent processors, 
//    Volume 1: General Techniques and Regular Problems,
//    Prentice Hall, 1988,
//    ISBN: 0-13-8230226,
//    LC: QA76.5.F627.
//
//  Local parameters:
//
//    Local, double DT, the time step.
//
//    Local, int ID, the MPI process ID.
//
//    Local, int N_GLOBAL, the total number of points.
//
//    Local, int N_LOCAL, the number of points visible to this process.
//
//    Local, int NSTEPS, the number of time steps.
//
//    Local, int P, the number of MPI processes.
//
{
  double dt = 0.00125;
  int i_global_hi;
  int i_global_lo;
  int id;
  int n_global = 401;
  int n_local;
  int nsteps = 4000;
  int p;
  double *u1_local;
  double wtime;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );

  id = MPI::COMM_WORLD.Get_rank ( );

  p = MPI::COMM_WORLD.Get_size ( );

  if ( id == 0 ) 
  {
    timestamp ( );
    cout << "\n";
    cout << "MPI_WAVE:\n";
    cout << "  C++ version.\n";
    cout << "  Estimate a solution of the wave equation using MPI.\n";
    cout << "\n";
    cout << "  Using " << p << " processes.\n";
    cout << "  Using a total of " << n_global << " points.\n";
    cout << "  Using " << nsteps << " time steps of size " << dt << "\n";
    cout << "  Computing final solution at time " << dt * nsteps << "\n";
  }

  wtime = MPI::Wtime ( );
//
//  Determine N_LOCAL.
//
  i_global_lo = (   id       * ( n_global - 1 ) ) / p;
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
  if ( 0 < id )
  {
    i_global_lo = i_global_lo - 1;
  }

  n_local = i_global_hi + 1 - i_global_lo;
//
//  Update N_LOCAL values.
//
  u1_local = update ( id, p, n_global, n_local, nsteps, dt );
//
//  Collect local values into global array.
//
  collect ( id, p, n_global, n_local, nsteps, dt, u1_local );
//
//  Report elapsed wallclock time.
//
  wtime = MPI::Wtime ( ) - wtime;
  if ( id == 0 )
  {
    cout << "\n";
    cout << "  Elapsed wallclock time was " << wtime << " seconds\n";
  }
//
//  Terminate MPI.
//
  MPI::Finalize ( );
//
//  Free memory.
//
  delete [] u1_local;
//
//  Terminate.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "WAVE_MPI:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }

  return 0;
}
//****************************************************************************80

double *update ( int id, int p, int n_global, int n_local, int nsteps, 
  double dt ) 

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE advances the solution a given number of time steps.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ID, the identifier of this process.
//
//    Input, int P, the number of processes.
//
//    Input, int N_GLOBAL, the total number of points.
//
//    Input, int N_LOCAL, the number of points visible to this process.
//
//    Input, int NSTEPS, the number of time steps.
//
//    Input, double DT, the size of the time step.
//
//    Output, double UPDATE[N_LOCAL], the portion of the solution
//    at the last time, as evaluated by this process.
//
{
  double alpha;
  double c;
  double dx;
  int i;
  int i_global;
  int i_global_hi;
  int i_global_lo;
  int i_local;
  int i_local_hi;
  int i_local_lo;
  int j;
  int ltor = 20;
  int rtol = 10;
  MPI::Status status;
  double t;
  double *u0_local;
  double *u1_local;
  double *u2_local;
  double x;
//
//  Determine the value of ALPHA.
//
  c = 1.0;
  dx = 1.0 / ( double ) ( n_global - 1 );
  alpha = c * dt / dx;

  if ( 1.0 <= fabs ( alpha ) )
  {
    if ( id == 0 )
    {
      cerr << "\n";
      cerr << "UPDATE - Warning!\n";
      cerr << "  1 <= |ALPHA| = | C * dT / dX |.\n";
      cerr << "  C = " << c << "\n";
      cerr << "  dT = " << dt << "\n";
      cerr << "  dX = " << dx << "\n";
      cerr << "  ALPHA = " << alpha << "\n";
      cerr << "  Computation will not be stable!\n";
    }
    MPI_Finalize ( );
    exit ( 1 );
  }

  i_global_lo = (   id       * ( n_global - 1 ) ) / p;
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
  if ( 0 < id )
  {
    i_global_lo = i_global_lo - 1;
  }

  i_local_lo = 0;
  i_local_hi = i_global_hi - i_global_lo;

  u0_local = new double[n_local];
  u1_local = new double[n_local];
  u2_local = new double[n_local];

  t = 0.0;
  for ( i_global = i_global_lo; i_global <= i_global_hi; i_global++ ) 
  {
    x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
    i_local = i_global - i_global_lo;
    u1_local[i_local] = exact ( x, t );
  }

  for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
  {
    u0_local[i_local] = u1_local[i_local];
  }
//
//  Take NSTEPS time steps.
//
  for ( i = 1; i <= nsteps; i++ )
  {
    t = dt * ( double ) i;
//
//  For the first time step, we need to use the initial derivative information.
//
    if ( i == 1 )
    {
      for ( i_local = i_local_lo + 1; i_local < i_local_hi; i_local++ ) 
      {
        i_global = i_global_lo + i_local;
        x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
        u2_local[i_local] = 
          +         0.5 * alpha * alpha   * u1_local[i_local-1]
          + ( 1.0 -       alpha * alpha ) * u1_local[i_local] 
          +         0.5 * alpha * alpha   * u1_local[i_local+1]
          +                            dt * dudt ( x, t );
      }
    }
//
//  After the first time step, we can use the previous two solution estimates.
//
    else
    {
      for ( i_local = i_local_lo + 1; i_local < i_local_hi; i_local++ ) 
      {
        u2_local[i_local] = 
          +               alpha * alpha   * u1_local[i_local-1]
          + 2.0 * ( 1.0 - alpha * alpha ) * u1_local[i_local] 
          +               alpha * alpha   * u1_local[i_local+1]
          -                                 u0_local[i_local];
      }
    }
//
//  Exchange data with "left-hand" neighbor. 
//
    if ( 0 < id ) 
    {
      MPI::COMM_WORLD.Send ( &u2_local[i_local_lo+1], 1, MPI::DOUBLE, id - 1, 
        rtol );
      MPI::COMM_WORLD.Recv ( &u2_local[i_local_lo],   1, MPI::DOUBLE, id - 1, 
        ltor, status );
    }
    else
    {
      x = 0.0;
      u2_local[i_local_lo] = exact ( x, t );
    }
//
//  Exchange data with "right-hand" neighbor.
//
    if ( id < p - 1 ) 
    {
      MPI::COMM_WORLD.Send ( &u2_local[i_local_hi-1], 1, MPI::DOUBLE, id + 1, 
        ltor );
      MPI::COMM_WORLD.Recv ( &u2_local[i_local_hi],   1, MPI::DOUBLE, id + 1, 
        rtol, status );
    }
    else
    {
      x = 1.0;
      u2_local[i_local_hi] = exact ( x, t );
    }
//
//  Shift data for next time step.
//
    for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
    {
      u0_local[i_local] = u1_local[i_local];
      u1_local[i_local] = u2_local[i_local];
    }
  }
//
//  Free memory.
//
  delete [] u0_local;
  delete [] u2_local;

  return u1_local;
}
//****************************************************************************80

void collect ( int id, int p, int n_global, int n_local, int nsteps, 
  double dt, double u_local[] ) 

//****************************************************************************80
//
//  Purpose:
//
//    COLLECT has workers send results to the master, which prints them.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ID, the identifier of this process.
//
//    Input, int P, the number of processes.
//
//    Input, int N_GLOBAL, the total number of points.
//
//    Input, int N_LOCAL, the number of points visible to this process.
//
//    Input, int NSTEPS, the number of time steps.
//
//    Input, double DT, the size of the time step.
//
//    Input, double U_LOCAL[N_LOCAL], the final solution estimate computed
//    by this process.
//
{
  int buffer[2];
  int collect1 = 10;
  int collect2 = 20;
  int i;
  int i_global;
  int i_global_hi;
  int i_global_lo;
  int i_local;
  int i_local_hi;
  int i_local_lo;
  int j;
  int n_local2;
  MPI::Status status;
  double t;
  double *u_global;
  double x;

  i_global_lo = (   id       * ( n_global - 1 ) ) / p;
  i_global_hi = ( ( id + 1 ) * ( n_global - 1 ) ) / p;
  if ( 0 < id )
  {
    i_global_lo = i_global_lo - 1;
  }

  i_local_lo = 0;
  i_local_hi = i_global_hi - i_global_lo;
//
//  Master collects worker results into the U_GLOBAL array.
//
  if ( id == 0 )
  {
//
//  Create the global array.
//
    u_global = new double[n_global];
//
//  Copy the master's results into the global array.
//
    for ( i_local = i_local_lo; i_local <= i_local_hi; i_local++ )
    {
      i_global = i_global_lo + i_local - i_local_lo;
      u_global[i_global] = u_local[i_local];
    }
//
//  Contact each worker.
//
    for ( i = 1; i < p; i++ ) 
    {
//
//  Message "collect1" contains the global index and number of values.
//
      MPI::COMM_WORLD.Recv ( buffer, 2, MPI::INT, i, collect1, status );
      i_global_lo = buffer[0];
      n_local2 = buffer[1];

      if ( i_global_lo < 0 )
      {
        cerr << "  Illegal I_GLOBAL_LO = " << i_global_lo << "\n";
        exit ( 1 );
      }
      else if ( n_global <= i_global_lo + n_local2 - 1 )
      {
        cerr << "  Illegal I_GLOBAL_LO + N_LOCAL2 = " 
             << i_global_lo + n_local2 << "\n";
        exit ( 1 );
      }
//
//  Message "collect2" contains the values.
//
      MPI::COMM_WORLD.Recv ( &u_global[i_global_lo], n_local2, MPI::DOUBLE, i, 
        collect2, status );
    }
//
//  Print the results.
//
    t = dt * ( double ) nsteps;
    cout << "\n";
    cout << "    I      X     F(X)   Exact\n";
    cout << "\n";
    for ( i_global = 0; i_global < n_global; i_global++) 
    {
      x = ( double ) ( i_global ) / ( double ) ( n_global - 1 );
      cout << "  " << setw(3) << i_global
           << "  " << setprecision(3) << setw(6) << x
           << "  " << setprecision(3) << setw(6) << u_global[i_global]
           << "  " << setprecision(3) << setw(6) << exact ( x, t ) << "\n";
    }

    delete [] u_global;
  }
//
//  Workers send results to process 0.
//
  else
  {
//
//  Message "collect1" contains the global index and number of values.
//
    buffer[0] = i_global_lo;
    buffer[1] = n_local;
    MPI::COMM_WORLD.Send ( buffer, 2, MPI::INT, 0, collect1 );
//
//  Message "collect2" contains the values.
//
    MPI::COMM_WORLD.Send ( u_local, n_local, MPI::DOUBLE, 0, collect2 );
  }

  return;
}
//****************************************************************************80

double exact ( double x, double t )

//****************************************************************************80
//
//  Purpose:
//
//    EXACT evaluates the exact solution
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the location.
//
//    Input, double T, the time.
//
//    Output, double EXACT, the value of the exact solution.
//
{
  const double c = 1.0;
  const double pi = 3.141592653589793;
  double value;

  value = sin ( 2.0 * pi * ( x - c * t ) );

  return value;
}
//****************************************************************************80

double dudt ( double x, double t )

//****************************************************************************80
//
//  Purpose:
//
//    DUDT evaluates the partial derivative dudt.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 November 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the location.
//
//    Input, double T, the time.
//
//    Output, double DUDT, the value of the time derivative of the solution.
//
{
  const double c = 1.0;
  const double pi = 3.141592653589793;
  double value;

  value = - 2.0 * pi * c * cos ( 2.0 * pi * ( x - c * t ) );

  return value;
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
