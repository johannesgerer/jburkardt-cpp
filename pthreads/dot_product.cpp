# include <cstdlib>
# include <iostream>
# include <pthread.h>

using namespace std;

# define NUMTHRDS 4
double sum;
double a[256], b[256];
int status;
int n = 256;
pthread_t thds[NUMTHRDS];
pthread_mutex_t mutex_sum;

int main ( int argc, char *argv[] );
void *dotprod ( void *arg );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DOT_PRODUCT.
//
//  Reference:
//
//    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
//    Using OpenMP: Portable Shared Memory Parallel Processing,
//    MIT Press, 2007,
//    ISBN13: 978-0262533027,
//    LC: QA76.642.C49.
//
{
  int i;

  pthread_attr_t attr;

  cout << "\n";
  cout << "DOT_PRODUCT:\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Demonstrate the use of Posix Threads by computing the\n";
  cout << "  dot product of two vectors.\n";
  cout << "\n";
  cout << "  The length of the vectors N = " << n << "\n";
  cout << "  The number of threads used is " << NUMTHRDS << "\n";

  for ( i = 0; i < n; i++ )
  {
    a[i] = i * 0.5;
    b[i] = i * 2.0;
  }

  pthread_mutex_init ( &mutex_sum, NULL );
  pthread_attr_init ( &attr );
  pthread_attr_setdetachstate ( &attr, PTHREAD_CREATE_JOINABLE );

  for ( i = 0; i < NUMTHRDS; i++ )
  {
    pthread_create ( &thds[i], &attr, dotprod, ( void * ) i );
  }

  pthread_attr_destroy ( &attr );

  for ( i = 0; i < NUMTHRDS; i++ )
  {
    pthread_join ( thds[i], ( void ** ) &status );
  }

  cout << "\n";
  cout << "  Sum = " << sum << "\n";

  pthread_mutex_destroy ( &mutex_sum );
  pthread_exit ( NULL );

  cout << "\n";
  cout << "DOT_PRODUCT:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void *dotprod ( void *arg )

//****************************************************************************80
//
//  Purpose:
//
//    DOT_PRODUCT is executed by each thread.
//
//  Reference:
//
//    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
//    Using OpenMP: Portable Shared Memory Parallel Processing,
//    MIT Press, 2007,
//    ISBN13: 978-0262533027,
//    LC: QA76.642.C49.
//
{
  int i;
  int my_first;
  int my_last;
  int myid;
  double sum_local;

  myid = ( int ) arg;
/*
  Determine the portion of the dot product to be computed by this thread.
*/
  my_first = myid * n / NUMTHRDS;
  my_last = ( myid + 1 ) * n / NUMTHRDS;
/*
  Compute a part of the dot product.
*/
  sum_local = 0;
  for ( i = my_first; i <= my_last; i++ )
  {
    sum_local = sum_local + a[i] * b[i];
  }
/*
  Lock the variable MUTEX_SUM, update it, and then unlock it.
*/
  pthread_mutex_lock ( &mutex_sum );
  sum = sum + sum_local;
  pthread_mutex_unlock ( &mutex_sum );

  pthread_exit ( ( void * ) 0 );
}
