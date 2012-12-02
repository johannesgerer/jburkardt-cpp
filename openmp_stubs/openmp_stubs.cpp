# include <cstdlib>
# include <iostream>
# include <ctime>

using namespace std;

# include "openmp_stubs.hpp"

//****************************************************************************80

void omp_destroy_lock ( omp_lock_t *lock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_DESTROY_LOCK destroys a simple lock.
//
//  Discussion:
//
//    The routine is intended to return the state of the lock to the 
//    uninitialized state.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, omp_lock_t *LOCK, the simple lock.
//
{
  *lock = INIT;

  return;
}
//****************************************************************************80

void omp_destroy_nest_lock ( omp_nest_lock_t *nlock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_DESTROY_NEST_LOCK destroys a nestable lock.
//
//  Discussion:
//
//    The routine is intended to return the state of the lock to the 
//    uninitialized state.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, omp_nest_lock_t *NLOCK, the nestable lock.
//
{
  nlock->owner = NOOWNER;
  nlock->count = UNLOCKED;

  return;
}
//****************************************************************************80

int omp_get_dynamic ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_DYNAMIC reports if dynamic adjustment of thread number is allowed.
//
//  Discussion:
//
//    The user can request dynamic thread adjustment by calling OMP_SET_DYNAMIC.
//
//    For this stub library, the value FALSE is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GET_DYNAMIC, is TRUE (1) if dynamic adjustment of thread
//    number has been enable, by default or by a user call, and FALSE (0)
//    otherwise.
//
{
  int value;

  value = 0;

  return value;
}
//****************************************************************************80

int omp_get_max_threads ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_MAX_THREADS returns the default number of threads.
//
//  Discussion:
//
//    If a parallel region is reached, and no number of threads has been
//    specified explicitly, there is a default number of threads that will
//    be used to form the new team.  That value is returned by this function.
//
//    For this stub library, the value 1 is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GET_MAX_THREADS, the default number of threads.
//
{
  int value;

  value = 1;

  return value;
}
//****************************************************************************80

int omp_get_nested ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_NESTED reports if nested parallelism has been enabled.
//
//  Discussion:
//
//    The user can request nested parallelism by calling OMP_SET_NESTED.
//
//    For this stub library, the value FALSE is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GET_NESTED, is TRUE (1) if nested parallelism has been
//    enable by default or by a user call, and FALSE (0) otherwise.
//
{
  int value;

  value = 0;

  return value;
}
//****************************************************************************80

int omp_get_num_procs ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_NUM_PROCS returns the number of processors available to the program.
//
//  Discussion:
//
//    For this stub library, the value 1 is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GEN_NUM_PROCS, the number of processors available.
//
{
  int value;

  value = 1;

  return value;
}
//****************************************************************************80

int omp_get_num_threads ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_NUM_THREADS returns the number of threads in the current team.
//
//  Discussion:
//
//    For this stub library, the value 1 is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GET_NUM_THREADS, the number of threads in the 
//    current team.
//
{
  int value;

  value = 1;

  return value;
}
//****************************************************************************80

int omp_get_thread_num ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_THREAD_NUM returns the thread number of a thread in a team.
//
//  Discussion:
//
//    Thread numbers start at 0.
//
//    If this function is not called from a parallel region, then only one
//    thread is executing, so the value returned is 0.
//
//    For this stub library, the value 0 is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_GET_THREAD_NUM, the thread number.
//
{
  int value;

  value = 0;

  return value;
}
//****************************************************************************80

double omp_get_wtick ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_WTICK returns the precision of the timer used by OMP_GET_WTIME.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, double OMP_GET_WTICK, the number of seconds between
//    successive "ticks" of the wall clock timer.
//
{
  double value;

  value = 1.0 / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double omp_get_wtime ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_GET_WTIME returns elapsed wall clock time in seconds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, double OMP_GET_WTIME, the current reading of the
//    wall clock timer.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

int omp_in_parallel ( void )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_IN_PARALLEL returns TRUE if the call is made from a parallel region.
//
//  Discussion:
//
//    For this stub library, the value FALSE is always returned.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, int OMP_IN_PARALLEL, is "TRUE" (1) if the routine was called
//    from a parallel region and "FALSE" (0) otherwise.
//
{
  int value;

  value = 0;

  return value;
}
//****************************************************************************80

void omp_init_lock ( omp_lock_t *lock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_INIT_LOCK initializes a simple lock.
//
//  Discussion:
//
//    This routine is intended to initialize the lock to the unlocked state.
//
//    For this stub library, the lock points to -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, omp_lock_t *LOCK, the simple lock.
//
{
  *lock = UNLOCKED;

  return;
}
//****************************************************************************80

void omp_init_nest_lock ( omp_nest_lock_t *nlock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_INIT_NEST_LOCK initializes a nestable lock.
//
//  Discussion:
//
//    This routine is intended to initialize the lock to the unlocked state.
//
//    For this stub library, the lock points to -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Output, omp_nest_lock_t *NLOCK, the lock.
//
{
  nlock->owner = NOOWNER;
  nlock->count = 0;

  return;
}
//****************************************************************************80

void omp_set_dynamic ( int dynamic_threads )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_SET_DYNAMIC enables dynamic adjustment of the number of threads.
//
//  Discussion:
//
//    If DYNAMIC_THREADS is TRUE, then the number of threads available for
//    execution in a parallel region may be dynamically adjusted.
//
//    For this stub library, the input value is ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input, int DYNAMIC_THREADS, is TRUE (1) if the user wishes to allow
//    dynamic adjustment of the number of threads available for execution
//    in any parallel region.
//
{
  return;
}
//****************************************************************************80

void omp_set_lock ( omp_lock_t *lock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_SET_LOCK sets a simple lock.
//
//  Discussion:
//
//    The lock must already have been initialized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_lock_t *LOCK, the simple lock.
//
{
  if ( *lock == UNLOCKED )
  {
    *lock = LOCKED;
  }
  else if ( *lock == LOCKED )
  {
    cout << "\n";
    cout << "OMP_SET_LOCK - Fatal error!\n";
    cout << "  Deadlock using lock variable.\n";
    exit ( 1 );
  }
  else
  {
    cout << "\n";
    cout << "OMP_SET_LOCK - Fatal error!\n";
    cout << "  Lock not initialized.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void omp_set_nest_lock ( omp_nest_lock_t *nlock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_SET_NEST_LOCK sets a nestable lock.
//
//  Discussion:
//
//    The lock must already have been initialized.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_nest_lock_t *NLOCK, the nestable lock.
//
{
  if ( nlock->owner == MASTER && 1 <= nlock->count )
  {
    nlock->count = nlock->count + 1;
  }
  else if ( nlock->owner == NOOWNER && nlock->count == 0 )
  {
    nlock->owner = MASTER;
    nlock->count = 1;
  }
  else
  {
    cout << "\n";
    cout << "OMP_SET_NEST_LOCK - Fatal error!\n";
    cout << "  Lock corrupted or not initialized.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void omp_set_nested ( int nested )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_SET_NESTED controls the use of nested parallelism.
//
//  Discussion:
//
//    For this stub library, the input value is ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input, int NESTED, is TRUE (1) if nested parallelism is to be enabled.
//
{
  return;
}
//****************************************************************************80

void omp_set_num_threads ( int num_threads )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_SET_NUM_THREADS sets the number of threads.
//
//  Discussion:
//
//    This routine sets the number of threads to be used in all subsequent
//    parallel regions for which an explicit number of threads is not
//    specified.
//
//    For this stub library, the input value is ignored.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input, int NUM_THREADS, the number of threads to be used in all
//    subsequent parallel regions for which an explicit number of threads
//    is not specified.  0 < NUM_THREADS.
//
{
  if ( num_threads <= 0 )
  {
    cout << "\n";
    cout << "OMP_SET_NUM_THREADS - Fatal error!\n";
    cout << "  Number of threads must be at least 1.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

int omp_test_lock ( omp_lock_t *lock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_TEST_LOCK tests a simple lock.
//
//  Discussion:
//
//    Calling this routine with an uninitialized lock causes an error.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_lock_t *LOCK, the simple lock.
//    If the lock was initialized but not set, it is set by this call.
//
//    Output, int OMP_TEST_LOCK:
//    TRUE (1), on input, the lock was initialized and not set;
//    FALSE (0), on input the lock was initialized, and set.
//
{
  int value;

  if ( *lock == UNLOCKED )
  {
    *lock = LOCKED;
    value = 1;
  }
  else if ( *lock == LOCKED )
  {
    value = 0;
  }
  else
  {
    cout << "\n";
    cout << "OMP_TEST_LOCK - Fatal error!\n";
    cout << "  Lock not initialized.\n";
    exit ( 1 );
  }

  return value;
}
//****************************************************************************80

int omp_test_nest_lock ( omp_nest_lock_t *nlock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_TEST_NEST_LOCK tests a nestable lock.
//
//  Discussion:
//
//    Calling this routine with an uninitialized lock causes an error.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_nest_lock_t NLOCK, the nestable lock.
//    If the lock was initialized but not set, it is set by this call.
//
//    Output, int OMP_TEST_NEST_LOCK, returns the new nesting count,
//    if the call was successful.  Otherwise, the value 0 is returned.
//
{
  int value;

  omp_set_nest_lock ( nlock );

  value = nlock->count;

  return value;
}
//****************************************************************************80

void omp_unset_lock ( omp_lock_t *lock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_UNSET_LOCK unsets a simple lock.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_lock_t *LOCK, the simple lock.
//
{
  if ( *lock == LOCKED )
  {
    *lock = UNLOCKED;
  }
  else if ( *lock == UNLOCKED )
  {
    cout << "\n";
    cout << "OMP_UNSET_LOCK - Fatal error!\n";
    cout << "  Attempt to unset lock, but lock was not set first.\n";
    exit ( 1 );
  }
  else
  {
    cout << "\n";
    cout << "OMP_UNSET_LOCK - Fatal error!\n";
    cout << "  Lock not initialized.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

void omp_unset_nest_lock ( omp_nest_lock_t *nlock )

//****************************************************************************80
//
//  Purpose:
//
//    OMP_UNSET_NEST_LOCK unsets a nestable lock.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    OpenMP Application Program Interface,
//    Version 2.5,
//    May 2005.
//
//  Parameters:
//
//    Input/output, omp_nest_lock_t *NLOCK, the nestable lock.
//
{
  if ( nlock->owner == NOOWNER && 1 <= nlock->count )
  {
    nlock->count = nlock->count -1;
    if ( nlock->count == 0 )
    {
      nlock->owner = NOOWNER;
    }
  }
  else if ( nlock->owner == NOOWNER && nlock->count == 0 )
  {
    cout << "\n";
    cout << "OMP_UNSET_NEST_LOCK - Fatal error!\n";
    cout << "  Nested lock not set.\n";
    exit ( 1 );
  }
  else
  {
    cout << "\n";
    cout << "OMP_UNSET_NEST_LOCK - Fatal error!\n";
    cout << "  Nested lock corrupted or not initialized.\n";
    exit ( 1 );
  }

  return;
}
