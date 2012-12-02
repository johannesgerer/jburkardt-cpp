# include <cstdlib>
# include <ctime>
# include <iostream>

using namespace std;

# include "mpi_stubs.hpp"

namespace MPI 
{
  // Define global COMM_WORLD communicator.  Should be initialized in Init
  Comm COMM_WORLD;

//****************************************************************************80

  void Finalize ( )

//****************************************************************************80
//
//  Purpose:
//
//    MPI::FINALIZE shuts down the MPI library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    None
//
  {
    return;
  }
//****************************************************************************80

  void Init ( int &argc, char **&argv )

//****************************************************************************80
//
//  Purpose:
//
//    MPI::INIT initializes the MPI library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    Input, int &ARGC, a reference to the program's ARGC parameter.
//
//    Input, char **&ARGV, a reference to the program's ARGV parameter.
//
  {
    return;
  }
//****************************************************************************80

  void Init ( )

//****************************************************************************80
//
//  Purpose:
//
//    MPI::INIT initializes the MPI library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    None
//
  {
    return;
  }
//****************************************************************************80

  double Wtick ( )

//****************************************************************************80
//
//  Purpose:
//
//    MPI::WTICK returns the number of seconds in one tick of the timer.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    Output, double WTICK, the number of seconds in one tick of the timer.
//
  {
    double value;

    value = 1.0 / ( double ) CLOCKS_PER_SEC;
  
    return value;
  }
//****************************************************************************80

  double Wtime ( )

//****************************************************************************80
//
//  Purpose:
//
//    MPI::WTIME returns the current wallclock time in seconds.
//
//  Discussion:
//
//    This "stub" version always returns a value that is 1.0 greater
//    than the value it returned on the previous call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    Output, double WTIME, the current wallclock time in seconds.
//
  {
    double value;
  
    value = ( double ) clock ( ) 
          / ( double ) CLOCKS_PER_SEC;
  
    return value;
  }
//****************************************************************************80

  int Comm::Get_size ( ) const
  
//****************************************************************************80
//
//  Purpose:
//
//    MPI::COMM::GET_SIZE reports the number of processes in a communicator.
//
//  Discussion:
//
//    This "stub" version always returns a value of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    Output, int GET_SIZE, the number of processors.
//
  {
    return 1;
  }
//****************************************************************************80

  int Comm::Get_rank ( ) const
  
//****************************************************************************80
//
//  Purpose:
//
//    MPI::COMM::GET_RANK reports the rank of the calling process.
//
//  Discussion:
//
//    This "stub" version always returns the value 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    William Gropp, Ewing Lusk, Anthony Skjellum,
//    Using MPI: 
//    Portable Parallel Programming with the Message-Passing Interface,
//    Second Edition,
//    MIT Press, 1999,
//    ISBN: 0262571323,
//    LC: QA76.642.G76.
//
//  Parameters:
//
//    Output, int GET_RANK, the rank of the process.
//
  {
    return 0;
  }
  
}
// 
//  END OF NAMESPACE MPI
//
