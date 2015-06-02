# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    HELLO has each thread print out its ID.
//
//  Discussion:
//
//    HELLO is a "Hello, World" program for OpenMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 May 2010
//
//  Author:
//
//    John Burkardt
//
{
  int id;
  double wtime;

  cout << "\n";
  cout << "HELLO_OPENMP\n";
  cout << "  C++/OpenMP version\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  wtime = omp_get_wtime ( );
//
//  Have each thread say hello
//
# pragma omp parallel \
  private ( id )
  {
    id = omp_get_thread_num ( );
    cout << "  This is process " << id << "\n";
  }
//
//  Finish up by measuring the elapsed time.
//
  wtime = omp_get_wtime ( ) - wtime;
//
//  Terminate.
//
  cout << "\n";
  cout << "HELLO_OPENMP\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  cout << "  Elapsed wall clock time = " << wtime << "\n";

  return 0;
}
