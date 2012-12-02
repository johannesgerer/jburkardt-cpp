# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "mpi.h"

int main ( int argc, char *argv[] );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUM.
//
//  Discussion:
//
//    SUM is an example of using the MPI message passing interface library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Forrest Hoffman,
//    Message Passing with MPI and PVM,
//    LINUX Magazine,
//    Volume 4, Number 4, April 2002, pages 38-41, 63.
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: Portable Parallel Programming with the
//    Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323.
//
{
# define N 100

  double array[N];
  int i;
  int id;
  int p;
  double PI = 3.141592653589793238462643;
  double seed;
  MPI::Status status;
  double sum;
  double sum_all;
  int tag;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get the number of processes.
//
  p = MPI::COMM_WORLD.Get_size ( );
//
//  Determine the rank of this process.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Say hello.
//
  if ( id == 0 )
  {
    timestamp ( );
    cout << "\n";
    cout << "SUM - Master process:\n";
    cout << "  C++ version\n";
    cout << "\n";
    cout << "  An MPI example program.\n";
    cout << "  The master process computes some coefficients,\n";
    cout << "  sends them to each worker process, which sums them.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  The number of processes available is " << p << "\n";
  }
//
//  The master process initializes the array.
//
  if ( id == 0 ) 
  {
    seed = 1.2345;

    for ( i = 0; i < N; i++ )
    {
      array[i] = ( double ) i * seed * PI;
    }
  }
// 
//  The master process broadcasts the computed initial values 
//  to all the other processes.
//
  MPI::COMM_WORLD.Bcast ( array, N, MPI::DOUBLE, 0 );
//
//  Each process adds up its entries.
//
  sum = 0.0;
  for ( i = 0; i < N; i++ )
  {
    sum = sum + array[i] * ( double ) id;
  }

  cout << "\n";
  cout << "SUM - Process " << id << ":\n";
  cout << "  My contribution to the sum is " << sum << "\n";
// 
//  Each worker process sends its sum back to the master process.
//
  if ( id != 0 ) 
  {
    MPI::COMM_WORLD.Send ( &sum, 1, MPI::DOUBLE, 0, 1 );
  }
  else 
  {
    sum_all = sum;
    for ( i = 1; i < p; i++ ) 
    {
      tag = 1;
      MPI::COMM_WORLD.Recv ( &sum, 1, MPI::DOUBLE, MPI::ANY_SOURCE, tag, 
        status );

      sum_all = sum_all + sum;
    }
  }
  if ( id == 0 )
  {
    cout << "\n";         
    cout << "SUM - Master process:\n";         
    cout << "  The total sum is " << sum_all << "\n";
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
    cout << "SUM - Master process:\n";       
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );      
  }
  return 0;

# undef N
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
