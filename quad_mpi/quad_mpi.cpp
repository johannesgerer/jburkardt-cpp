# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <iomanip>
# include <cstring>
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
//    MAIN is the main program for QUAD_MPI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int id;
  double my_a;
  double my_b;
  int my_n;
  double my_total;
  int n;
  int p;
  int q;
  int source;
  MPI::Status status;
  int tag;
  int target;
  double total;
  double wtime;
  double x;

  a =  0.0;
  b = 10.0;
  n = 10000000;
  exact = 0.49936338107645674464;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get this processor's ID.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Get the number of processors.
//
  p = MPI::COMM_WORLD.Get_size ( );

  if ( id == 0 )
  {
//
//  We want N to be the total number of evaluations.
//  If necessary, we adjust N to be divisible by the number of processors.
//
    my_n = n / ( p - 1 );
    n = ( p - 1 ) * my_n;

    wtime = MPI::Wtime ( );

    timestamp ( );
    cout << "\n";
    cout << "QUAD_MPI\n";
    cout << "  C++/MPI version\n";
    cout << "  Estimate an integral of f(x) from A to B.\n";
    cout << "  f(x) = 50 / (pi * ( 2500 * x * x + 1 ) )\n";
    cout << "\n";
    cout << "  A = " << a << "\n";
    cout << "  B = " << b << "\n";
    cout << "  N = " << n << "\n";
    cout << "  EXACT = " << setw(24) << setprecision(16) << exact << "\n";
    cout << "\n";
    cout << "  Use MPI to divide the computation among\n";
    cout << "  multiple processes.\n";
  }

  source = 0;
  MPI::COMM_WORLD.Bcast ( &my_n, 1, MPI::INT, source );
//
//  Process 0 assigns each process a subinterval of [A,B].
//
  if ( id == 0 )
  {
    for ( q = 1; q <= p - 1; q++ )
    {
      my_a = ( ( double ) ( p - q     ) * a   
             + ( double ) (     q - 1 ) * b ) 
             / ( double ) ( p     - 1 );

      target = q;
      tag = 1;
      MPI::COMM_WORLD.Send ( &my_a, 1, MPI::DOUBLE, target, tag );

      my_b = ( ( double ) ( p - q - 1 ) * a   
             + ( double ) (     q     ) * b ) 
             / ( double ) ( p     - 1 );

      target = q;
      tag = 2;
      MPI::COMM_WORLD.Send ( &my_b, 1, MPI::DOUBLE, target, tag );
    }
    total = 0.0;
    my_total = 0.0;
  }
//
//  Processes receive MY_A, MY_B, and compute their part of the integral.
//
  else
  {
    source = 0;
    tag = 1;
    MPI::COMM_WORLD.Recv ( &my_a, 1, MPI::DOUBLE, source, tag, status );

    source = 0;
    tag = 2;
    MPI::COMM_WORLD.Recv ( &my_b, 1, MPI::DOUBLE, source, tag, status );

    my_total = 0.0;
    for ( i = 1; i <= my_n; i++ )
    {
      x = ( ( double ) ( my_n - i     ) * my_a 
          + ( double ) (        i - 1 ) * my_b )
          / ( double ) ( my_n     - 1 );
      my_total = my_total + f ( x );
    }

    my_total = ( my_b - my_a ) * my_total / ( double ) ( my_n );

    cout << "  Process " << id << " contributed MY_TOTAL = " 
         << my_total << "\n";
  }
//
//  Each process sends its value to the master process.
//
  MPI::COMM_WORLD.Reduce ( &my_total, &total, 1, MPI::DOUBLE, MPI::SUM, 0 );
//
//  Compute the weighted estimate.
//
  if ( id == 0 )
  {
    error = fabs ( total - exact );
    wtime = MPI::Wtime ( ) - wtime;

    cout << "\n";
    cout << "  Estimate = " << setw(24) << setprecision(16) << total << "\n";
    cout << "  Error = " << error << "\n";
    cout << "  Time = " << wtime << "\n";
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
    cout << "QUAD_MPI:\n";
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
//    F evaluates the function.
//
{
  double pi;
  double value;

  pi = 3.141592653589793;
  value = 50.0 / ( pi * ( 2500.0 * x * x + 1.0 ) );

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
