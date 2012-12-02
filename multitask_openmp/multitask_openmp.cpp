# include <cstdlib>
# include <iostream>
# include <cmath>
# include <ctime>
# include <omp.h>

using namespace std;

int main ( );
int *prime_table ( int prime_num );
double *sine_table ( int sine_num );
void timestamp ( );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MULTITASK_OPENMP.
//
//  Discussion:
//
//    This program demonstrates how OpenMP can be used for multitasking, that 
//    is, a simple kind of parallel processing in which a certain number of 
//    perhaps quite unrelated tasks must be done.
//
//    The OpenMP SECTIONS directive identifies the portion of the program where
//    the code for these tasks is given.
//
//    The OpenMP SECTION directive is used repeatedly to divide this area of
//    the program into independent tasks.
//
//    The code will get the benefit of parallel processing up to the point where
//    there are as many threads as there are tasks.
//
//    The code will get a substantial speedup if the tasks take roughly the
//    same amount of time.  However, if one task takes substantially more time
//    than the others, this results in a limit to the parallel speedup that is
//    possible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2011
//
//  Author:
//
//    John Burkardt
//
//
{
  int prime_num;
  int *primes;
  int sine_num;
  double *sines;
  double wtime;
  double wtime1;
  double wtime2;

  timestamp ( );
  cout << "\n";
  cout << "MULTITASK_OPENMP:\n";
  cout << "  C++/OpenMP version\n";
  cout << "  Demonstrate how OpenMP can \"multitask\" by using the\n";
  cout << "  SECTIONS directive to carry out several tasks in parallel.\n";

  prime_num = 20000;
  sine_num = 20000;

  wtime = omp_get_wtime ( );

# pragma omp parallel shared ( prime_num, primes, sine_num, sines )
{
  # pragma omp sections
  {
    # pragma omp section
    {
      wtime1 = omp_get_wtime ( );
      primes = prime_table ( prime_num );
      wtime1 = omp_get_wtime ( ) - wtime1;
    }
    # pragma omp section
    {
      wtime2 = omp_get_wtime ( );
      sines = sine_table ( sine_num );
      wtime2 = omp_get_wtime ( ) - wtime2;
    }
  }
}
  wtime = omp_get_wtime ( ) - wtime;

  cout << "\n";
  cout << "  Number of primes computed was " << prime_num << "\n";
  cout << "  Last prime was " << primes[prime_num-1] << "\n";
  cout << "  Number of sines computed was " << sine_num << "\n";
  cout << "  Last sine computed was " << sines[sine_num-1] << "\n";
  cout << "\n";
  cout << "  Elapsed time = " << wtime << "\n";
  cout << "  Task 1 time = " << wtime1 << "\n";
  cout << "  Task 2 time = " << wtime2 << "\n";

  free ( primes );
  free ( sines );
//
//  Terminate.
//
  cout << "\n";
  cout << "MULTITASK_OPENMP:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int *prime_table ( int prime_num )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_TABLE computes a table of the first PRIME_NUM prime numbers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PRIME_NUM, the number of primes to compute.
//
//    Output, int PRIME_TABLE[PRIME_NUM], the computed primes.
//
{
  int i;
  int j;
  int p;
  int prime;
  int *primes;

  primes = ( int * ) malloc ( prime_num * sizeof ( int ) );

  i = 2;
  p = 0;

  while ( p < prime_num )
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
      primes[p] = i;
      p = p + 1;
    }
    i = i + 1;
  }

  return primes;
}
//****************************************************************************80

double *sine_table ( int sine_num )

//****************************************************************************80
//
//  Purpose:
//
//    SINE_TABLE computes a table of sines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SINE_NUM, the number of sines to compute.
//
//    Output, double SINE_TABLE[SINE_NUM], the sines.
//
{
  double a;
  int i;
  int j;
  double pi = 3.141592653589793;
  double *sines;

  sines = ( double * ) malloc ( sine_num * sizeof ( double ) );

  for ( i = 0; i < sine_num; i++ )
  {
    sines[i] = 0.0;
    for ( j = 0; j <= i; j++ )
    {
      a = ( double ) ( j ) * pi / ( double ) ( sine_num - 1 );
      sines[i] = sines[i] + sin ( a );
    }
  }

  return sines;
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
