# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "mpi.h"

int main ( int argc, char *argv[] );
double f ( double x );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUADRATURE.
//
//  Discussion:
//
//    QUADRATURE estimates an integral using quadrature.
//
//    The integral of F(X) = 4 / ( 1 + X * X ) from 0 to 1 is PI.
//
//    We break up the interval [0,1] into N subintervals, evaluate
//    F(X) at the midpoint of each subinterval, and multiply the
//    sum of these values by N to get an estimate for the integral.
//
//    If we have M processes available because we are using MPI, then
//    we can ask processes 0, 1, 2, ... M-1 to handle the subintervals
//    in the following order:
//
//          0      1       2            M-1  <-- Process numbers begin at 0
//     ------ ------  ------  -----  ------
//          1      2       3    ...       M
//        M+1    M+2     M+3    ...     2*M
//      2*M+1    2*M+2 2*M+3    ...     3*M
//                              
//    and so on up to subinterval N.  The partial sums collected by 
//    each process are then sent to the master process to be added 
//    together to get the estimated integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: Portable Parallel Programming with the
//    Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323.
//
{
  double h;
  int i;
  int id;
  int ierr;
  int n;
  int n_part;
  int p;
  double q;
  double q_diff;
  double q_exact = 3.141592653589793238462643;
  double q_part;
  double wtime_diff;
  double wtime_end;
  double wtime_start;
  double x;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get the number of processes.
//
  p = MPI::COMM_WORLD.Get_size (  );
//
//  Determine this processes's rank.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Say hello.
//
  if ( id == 0 ) 
  {
    timestamp ( );
    cout << "\n";
    cout << "QUADRATURE - Master process:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  An MPI proram to estimate an integral.\n";
    cout << "\n";
    cout << "  The number of processes available is " << p << "\n";
  }
//
//  Record the starting time.
//
  if ( id == 0 ) 
  {
    wtime_start = MPI::Wtime ( );
  }
//
//  Normally, we will assume, the number of intervals N is read in
//  by process 0.  We'll just simulate this by assigning a value to
//  N, but only on process 0.
//
  if ( id == 0 )
  {
    n = 1000;
    cout << "\n";
    cout << "QUADRATURE - Master process:\n";
    cout << "  Number of intervals used is " << n << "\n";
  }
//
//  Broadcast the number of intervals to the other processes.
//
  MPI::COMM_WORLD.Bcast ( &n, 1, MPI::INT, 0 );
//
//  Integrate F(X) over a subinterval determined by your process ID.
//
  h = 1.0 / ( double ) n;

  q_part = 0.0;
  n_part = 0;

  for ( i = id + 1; i <= n; i = i + p ) 
  {
    x = ( double ) ( 2 * i - 1 )
      / ( double ) ( 2 * n     );

    n_part = n_part + 1;
    q_part = q_part + f ( x );
  }

  q_part = q_part * h;

  cout << "\n";
  cout << "QUADRATURE - Process " << id << "\n";
  cout << "  Points used = " << n_part << "\n";
  cout << "  Estimate " << q_part << ".\n";
//
//  Each process sends its value Q_PART to the MASTER process, to be added to
//  the global result Q.
//
  MPI::COMM_WORLD.Reduce ( &q_part, &q, 1, MPI::DOUBLE, MPI::SUM, 0 );
//
//  The master process scales the sum and prints the result.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "QUADRATURE - Master process:\n";         
    cout << "  The estimate for PI is " << q << "\n";
    cout << "  The exact value is     " << q_exact << "\n";
    q_diff = fabs ( q - q_exact );
    cout << "  The error is           " << q_diff << "\n";

    wtime_end = MPI_Wtime ( );
    wtime_diff = wtime_end - wtime_start;

    cout << "\n";
    cout << "  Wall clock elapsed seconds = " << wtime_diff << "\n";
  }
//
//  Terminate MPI.
//
  MPI::Finalize ( );
//
//  Terminate.
//
  if ( id == 0 ) 
  {
    cout << "\n";
    cout << "QUADRATURE - Master process:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }
  return 0;
}
//****************************************************************************80

double f ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F evaluates the function F(X) which we are integrating.
//
//  Discussion:
//
//    Integral ( 0 <= X <= 1 ) 4/(1+X*X) dX = PI
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double F, the value of the function.
//
{
  double value;

  value = 4.0 / ( 1.0 + x * x );

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
