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
//    MAIN is the main program for DAY1.
//
//  Discussion:
//
//    DAY1 is exercise 3 for first day of the MPI workshop
//
//    The instructions say:
//
//    Process 1 computes the squares of the first 200 integers.
//    It sends this data to process 3.
//
//    Process 3 should divide the integers between 20 and 119 by 53,
//    getting a real result, and passes this data back to process 1.
//
//    * I presume the first 200 integers are the numbers 0 through 199.
//
//    * The instructions literally mean that process 3 should look
//      at integers whose VALUES are between 20 and 119.  I doubt that
//      is what the instructor meant, but it's more interesting than
//      simply picking the entries with index between 20 and 119,
//      so that's what I'll do.
//
//    * It is also not completely clear whether only the selected data
//      should be sent back, or the entire array.  Again, it is more
//      interesting to send back only part of the data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
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
//  Modified:
//
//    26 October 2011
//
//  Author:
//
//    John Burkardt
//
{
# define I_DIM 200
# define R_DIM 200

  int count;
  int count2;
  int dest;
  int i;
  int i_buffer[I_DIM];
  int id;
  int p;
  float r_buffer[R_DIM];
  int source;
  MPI::Status status;
  int tag;
//
//  Initialize MPI.
//
  MPI::Init ( argc, argv );
//
//  Determine this process's rank.
//
  id = MPI::COMM_WORLD.Get_rank ( );
//
//  Get the number of processes.
//
  p = MPI::COMM_WORLD.Get_size ( );
//
//  Have Process 0 say hello.
//
  if ( id == 0 )
  {
    timestamp ( );
    cout << "\n";
    cout << "DAY1:\n";
    cout << "  C++ version\n";
    cout << "  An MPI example program.\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
    cout << "\n";
    cout << "  The number of processes available is " << p << "\n";
  }
//
//  If we don't have at least 4 processes, then bail out now.
//
  if ( p < 4 )
  {
    cout << "\n";
    cout << "DAY1 - Process " << id << ".\n";
    cout << "  Not enough processes for this task!\n";
    cout << "  Bailing out now!\n";
    MPI::Finalize ( );
    return 1;
  }
//
//  Process 1 knows that it will generate 200 integers, and may receive no more
//  than 200 reals.
//
  if ( id == 1 )
  {
    count = 200;

    for ( i = 0; i < count; i++ ) 
    {
      i_buffer[i] = i * i;
    }

    dest = 3;
    tag = 1;

    MPI::COMM_WORLD.Send ( i_buffer, count, MPI::INT, dest, tag );

    cout << "P:" << id << " sent " << count 
         << " integers to process " << dest << ".\n";

    source = 3;
    tag = 2;

    MPI::COMM_WORLD.Recv ( r_buffer, R_DIM, MPI::FLOAT, source, tag, status );

    cout << "P:" << id << " received real values from process 3.\n";

    count = status.Get_count ( MPI::FLOAT );

    cout << "P:" << id << " Number of real values received is "
         << count << ".\n";

    cout << "P:" << id << " First 3 values = "
         << r_buffer[0] << "  "
         << r_buffer[1] << "  "
         << r_buffer[2] << "\n";
  }
//
//  Process 3 receives the integer data from process 1, selects some of the data, does
//  a real computation on it, and sends that part back to process 1.
//
  else if ( id == 3 ) 
  {
    source = 1;
    tag = 1;

    MPI::COMM_WORLD.Recv ( i_buffer, I_DIM, MPI::INT, source, tag, status );

    cout << "\n";
    cout << "P:" << id << " received integer values from process 1.\n";

    count = status.Get_count ( MPI::INT );

    cout << "P:" << id << " - Number of integers received is " 
         << count << ".\n";

    cout << "P:" << id << " First 3 values = "
         << i_buffer[0] << "  "
         << i_buffer[1] << "  "
         << i_buffer[2] << "\n";

    count2 = 0;
     
    for ( i = 0; i < count; i++ ) 
    {
      if ( 20 <= i_buffer[i] && i_buffer[i] <= 119 ) 
      {

        r_buffer[count2] = ( float ) i_buffer[i] / 53.0E+00;
        count2 = count2 + 1;

        if ( count2 <= 3 ) 
        {
          cout << "P:" << id << " Input integer " << i_buffer[i]
               << " becomes " << r_buffer[count2-1] << ".\n";
        }

      }
    }

    dest = 1;
    tag = 2;
  
    MPI::COMM_WORLD.Send ( r_buffer, count2, MPI::FLOAT, dest, tag );

    cout << "P:" << id << " sent " << count2 << " reals to process "
         << dest << ".\n";
  }
  else
  {
    cout << "\n";
    cout << "P:" << id << " - MPI has no work for me!\n";
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
    cout << "DAY1:\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp ( );
  }
  return 0;

# undef I_DIM
# undef R_DIM
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
