# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "mpi.h"

int main ( int argc, char *argv[] );
void ring_io ( int p, int id );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RING_MPI.
//
//  Discussion:
//
//    RING_MPI sends messages of various size from process 0 to 1 to 2 to
//    ...the last process and then back to 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Pacheco,
//    Parallel Programming with MPI,
//    Morgan Kaufman, 1996,
//    ISBN: 1558603395,
//    LC: QA76.642.P3.
//
{
  int error;
  int id;
  int p;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Get the number of processes.
//
  p = MPI::COMM_WORLD.Get_size ( );
//
//  Get the individual process ID.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Print a message.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "RING_MPI:\n";
    cout << "  C++/MPI version\n";
    cout << "  Measure time required to transmit data around\n";
    cout << "  a ring of processes\n";
    cout << "\n";
    cout << "  The number of processes is " << p << "\n";
  }

  ring_io ( p, id );
//
//  Shut down MPI.
//
  MPI::Finalize ( );
//
//  Terminate.
//
  if ( id == 0 )
  {
    cout << "\n";
    cout << "RING_MPI:\n";
    cout << "  Normal end of execution.\n";
  }
  return 0;
}
//****************************************************************************80

void ring_io ( int p, int id )

//****************************************************************************80
//
//  Purpose:
//
//    RING_IO carries out the tasks of process ID, of a total of P processes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Pacheco,
//    Parallel Programming with MPI,
//    Morgan Kaufman, 1996,
//    ISBN: 1558603395,
//    LC: QA76.642.P3.
//
{
  int dest;
  int i;
  int j;
  int n;
  int n_test[5] = { 100, 1000, 10000, 100000, 1000000 };
  int n_test_num = 5;
  int source;
  MPI::Status status;
  double tave;
  int test;
  int test_num = 10;
  double tmax;
  double tmin;
  double wtime;
  double *x;

  if ( id == 0 )
  {
    cout << "\n";
    cout << "  Timings based on " << test_num << " experiments\n";
    cout << "  N double precision values were sent\n";
    cout << "  in a ring transmission starting and ending at process 0\n";
    cout << "  and using a total of " << p << " processes.\n";
    cout << "\n";
    cout << "         N           T min           T ave           T max\n";
    cout << "\n";
  }
//
//  Choose message size.
//
  for ( i = 0; i < n_test_num; i++ )
  {    
    n = n_test[i];

    x = new double[n];
//
//  Process 0 sends very first message, 
//  then waits to receive the "echo" that has gone around the world.
//
    if ( id == 0 )
    {
      dest = 1;
      source = p - 1;

      tave = 0.0;
      tmin = 1.0E+30;
      tmax = 0.0;

      for ( test = 1; test <= test_num; test++ )
      {
//
//  Just in case, set the entries of X in a way that identifies
//  which iteration of the test is being carried out.
//
        for ( j = 0; j < n; j++ )
        {
          x[j] = ( double ) ( test + j );
        }

        wtime = MPI::Wtime ( );
        MPI::COMM_WORLD.Send ( x, n, MPI::DOUBLE, dest,   0 );
        MPI::COMM_WORLD.Recv ( x, n, MPI::DOUBLE, source, 0, status );
        wtime = MPI::Wtime ( ) - wtime;
//
//  Record the time it took.
//
        tave = tave + wtime;
        if ( wtime < tmin )
        {
          tmin = wtime;
        }
        if ( wtime > tmax )
        {
          tmax = wtime;
        }
      }

      tave = tave / ( double ) ( test_num );

      cout << "  " << setw(8) << n
           << "  " << setw(14) << tmin
           << "  " << setw(14) << tave
           << "  " << setw(14) << tmax << "\n";
    }
//
//  Worker ID must receive first from ID-1, then send to ID+1.
//
    else
    {
      source = id - 1;
      dest = ( ( id + 1 ) % p );
 
      for ( test = 1; test <= test_num; test++ )
      {
        MPI::COMM_WORLD.Recv ( x, n, MPI::DOUBLE, source, 0, status );
        MPI::COMM_WORLD.Send ( x, n, MPI::DOUBLE, dest,   0 );
      }
    }
    delete [] x;
  }

  return;
}
