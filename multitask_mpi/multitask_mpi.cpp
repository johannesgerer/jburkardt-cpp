# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>

# include "mpi.h"

using namespace std;

int main ( int argc, char *argv[] );
void p0_set_input ( int *input1, int *input2 );
void p0_send_input ( int input1, int input2 );
void p0_receive_output ( int *output1, int *output2 );
int p1_receive_input ( );
int p1_compute_output ( int input1 );
void p1_send_output ( int output1 );
int p2_receive_input ( );
int p2_compute_output ( int input2 );
void p2_send_output ( int output2 );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MPI_MULTITASK.
//
//  Discussion:
//
//    Message tag 1: P0 sends input to P1
//    Message tag 2: P0 sends input to P2
//    Message tag 3: P1 sends output to P0.
//    Message tag 4: P2 sends output to P0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  int id;
  int input1;
  int input2;
  int output1;
  int output2;
  int p;
  double wtime;
//
//  Process 0 is the "monitor".
//  It chooses the inputs, and sends them to the workers.
//  It waits for the outputs.
//  It plots the outputs.
//
  MPI::Init ( argc, argv );

  id = MPI::COMM_WORLD.Get_rank ( );

  p = MPI::COMM_WORLD.Get_size ( );
//
//  Make sure we have enough processes.
//
  if ( p < 3 )
  {
    printf ( "\n" );
    printf ( "MPI_MULTITASK - Fatal error!\n" );
    printf ( "  Number of available processes must be at least 3!\n" );
    MPI::Finalize ( );
    exit ( 1 );
  }
//
//  Run program P0 on process 0, and so on.
//
  if ( id == 0 )
  {
    timestamp ( );

    printf ( "\n" );
    printf ( "MPI_MULTITASK:\n" );
    printf ( "  C++ / MPI version\n" );

    wtime = MPI::Wtime ( );

    p0_set_input ( &input1, &input2 );
    p0_send_input ( input1, input2 );
    p0_receive_output ( &output1, &output2 );

    wtime = MPI::Wtime ( ) - wtime;
    printf ( "  Process 0 time = %g\n", wtime );

    MPI::Finalize ( );

    printf ( "\n" );
    printf ( "MPI_MULTITASK:\n" );
    printf ( "  Normal end of execution.\n" );

    timestamp ( );
  }
//
//  Process 1 works on task 1.
//  It receives input from process 0.
//  It computes the output.
//  It sends the output to process 0.
//
  else if ( id == 1 )
  {
    wtime = MPI::Wtime ( );
    input1 = p1_receive_input ( );
    output1 = p1_compute_output ( input1 );
    p1_send_output ( output1 );
    wtime = MPI::Wtime ( ) - wtime;
    printf ( "  Process 1 time = %g\n", wtime );
    MPI::Finalize ( );
  }
//
//  Process 2 works on task 2.
//  It receives input from process 0.
//  It computes the output.
//  It sends the output to process 0.
//
  else if ( id == 2 )
  {
    wtime = MPI::Wtime ( );
    input2 = p2_receive_input ( );
    output2 = p2_compute_output ( input2 );
    p2_send_output ( output2 );
    wtime = MPI::Wtime ( ) - wtime;
    printf ( "  Process 2 time = %g\n", wtime );
    MPI::Finalize ( );
  }
  return 0;
}
//****************************************************************************80

void p0_set_input ( int *input1, int *input2 )

//****************************************************************************80
//
//  Purpose:
//
//    P0_SET_INPUT sets input.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int *INPUT1, *INPUT2, the values of two
//    inputs used by tasks 1 and 2.
//
{
  *input1 = 10000000;
  *input2 = 100000;

  printf ( "\n" );
  printf ( "P0_SET_PARAMETERS:\n" );
  printf ( "  Set INPUT1 = %d\n", *input1 );
  printf ( "      INPUT2 = %d\n", *input2 );

  return;
}
//****************************************************************************80

void p0_send_input ( int input1, int input2 )

//****************************************************************************80
//
//  Purpose:
//
//    P0_SEND_INPUT sends input to processes 1 and 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int INPUT1, INPUT2, the values of two
//    inputs used by tasks 1 and 2.
//
{
  int id;
  int tag;

  id = 1;
  tag = 1;
  MPI::COMM_WORLD.Send ( &input1, 1, MPI::INT, id, tag );

  id = 2;
  tag = 2;
  MPI::COMM_WORLD.Send ( &input2, 1, MPI::INT, id, tag );

  return;
}
//****************************************************************************80

void p0_receive_output ( int *output1, int *output2 )

//****************************************************************************80
//
//  Purpose:
//
//    P0_RECEIVE_OUTPUT receives output from processes 1 and 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int OUTPUT1, OUTPUT2, the values of the
//    outputs of tasks 1 and 2.
//
{
  int output;
  int output_received;
  int source;
  MPI::Status status;

  output_received = 0;
//
//  Loop until every worker has checked in.
//
  while ( output_received < 2 )
  {
//
//  Receive the next message that arrives.
//
    MPI::COMM_WORLD.Recv ( &output, 1, MPI::INT, MPI::ANY_SOURCE, MPI::ANY_TAG, 
      status );
//
//  The actual source of the message is saved in STATUS.
//
    source = status.Get_source ( );
//
//  Save the value in OUTPUT1 or OUTPUT2.
//
    if ( source == 1 )
    {
      *output1 = output;
    }
    else
    {
      *output2 = output;
    }
    output_received = output_received + 1;
  }

  printf ( "\n" );
  printf ( "  Process 1 returned OUTPUT1 = %d\n", *output1 );
  printf ( "  Process 2 returned OUTPUT2 = %d\n", *output2 );

  return;
}
//****************************************************************************80

int p1_receive_input ( )

//****************************************************************************80
//
//  Purpose:
//
//    P1_RECEIVE_INPUT receives input from process 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P1_RECEIVE_INPUT, the value of the parameter.
//
{
  int id;
  int input1;
  MPI::Status status;
  int tag;

  id = 0;
  tag = 1;
  MPI::COMM_WORLD.Recv ( &input1, 1, MPI::INT, id, tag, status );

  return input1;
}
//****************************************************************************80

int p1_compute_output ( int input1 )

//****************************************************************************80
//
//  Purpose:
//
//    P1_COMPUTE_OUTPUT carries out computation number 1.
//
//  Discussion:
//
//    No MPI calls occur in this function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int INPUT1, the problem input.
//
//    Output, int P1_COMPUTE_OUTPUT1, the problem output.
//
{
  int i;
  int j;
  int k;
  int output1;

  output1 = 0;

  for ( i = 2; i <= input1; i++ )
  {
    j = i;
    k = 0;

    while ( 1 < j )
    {
      if ( ( j % 2 ) == 0 )
      {
        j = j / 2;
      }
      else
      {
        j = 3 * j + 1;
      }
      k = k + 1;
    }
    if ( output1 < k )
    {
      output1 = k;
    }
  }
  return output1;
}
//****************************************************************************80

void p1_send_output ( int output1 )

//****************************************************************************80
//
//  Purpose:
//
//    P1_SEND_OUTPUT sends output to process 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OUTPUT1, the problem output.
//
{
  int id;
  int tag;

  id = 0;
  tag = 3;
  MPI::COMM_WORLD.Send ( &output1, 1, MPI::INT, id, tag );

  return;
}
//****************************************************************************80

int p2_receive_input ( )

//****************************************************************************80
//
//  Purpose:
//
//    P2_RECEIVE_INPUT receives input from process 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int P2_RECEIVE_INPUT, the value of the parameter.
//
{
  int id;
  int input2;
  MPI::Status status;
  int tag;

  id = 0;
  tag = 2;
  MPI::COMM_WORLD.Recv ( &input2, 1, MPI::INT, id, tag, status );

  return input2;
}
//****************************************************************************80

int p2_compute_output ( int input2 )

//****************************************************************************80
//
//  Purpose:
//
//    P2_COMPUTE_OUTPUT carries out computation number 2.
//
//  Discussion:
//
//    No MPI calls occur in this function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int INPUT2, the problem input.
//
//    Output, int P2_COMPUTE_OUTPUT, the problem output.
//
{
  int i;
  int j;
  int output2;
  int prime;

  output2 = 0;

  for ( i = 2; i <= input2; i++ )
  {
    prime = 1;
    for ( j = 2; j < i; j++ )
    {
      if ( ( i % j ) == 0 )
      {
        prime = 0;
        break;
      }
    }
    if ( prime )
    {
      output2 = output2 + 1;
    }
  }
  return output2;
}
//****************************************************************************80

void p2_send_output ( int output2 )

//****************************************************************************80
//
//  Purpose:
//
//    P2_SEND_OUTPUT sends output to process 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int OUTPUT2, the problem output.
//
{
  int id;
  int tag;

  id = 0;
  tag = 4;
  MPI::COMM_WORLD.Send ( &output2, 1, MPI::INT, id, tag );

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
