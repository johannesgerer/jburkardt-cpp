# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
int i4_power ( int i, int j );
void i4mat_memory_test ( int n_log );
void i4vec_memory_test ( int n_log );
float r4_cpu_time ( );
float r4_real_time ( );
void r4mat_memory_test ( int n_log );
void r4vec_memory_test ( int n_log );
double r8_cpu_time ( );
double r8_real_time ( );
void r8mat_memory_test ( int n_log );
void r8vec_memory_test ( int n_log );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MEMORY_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int n_log;
  int n_log_max;
  int n_log_min;

  timestamp ( );

  cout << "\n";
  cout << "MEMORY_TEST\n";
  cout << "  C++ version\n";
  cout << "  Try to see how big vectors and matrices can be.\n";

  if ( argc <= 1 )
  {
    n_log_min = 0;
    cout << "\n";
    cout << "  Using default value of N_LOG_MIN = " << n_log_min << "\n";
  }
  else
  {
    n_log_min = atoi ( argv[1] );
    cout << "\n";
    cout << "  User value of N_LOG_MIN = " << n_log_min << "\n";
  }
  if ( argc <= 2 )
  {
    n_log_max = 27;
    cout << "\n";
    cout << "  Using default value of N_LOG_MAX = " << n_log_max << "\n";
  }
  else
  {
    n_log_max = atoi ( argv[2] );
    cout << "\n";
    cout << "  User value of N_LOG_MAX = " << n_log_max << "\n";
  }
//
//  I4VEC test.
//
  cout << "\n";
  cout << "I4VEC Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    i4vec_memory_test ( n_log );
  }
//
//  R4VEC test.
//
  cout << "\n";
  cout << "R4VEC Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r4vec_memory_test ( n_log );
  }
//
//  R8VEC test.
//
  cout << "\n";
  cout << "R8VEC Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r8vec_memory_test ( n_log );
  }
//
//  I4MAT test.
//
  cout << "\n";
  cout << "I4MAT Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N            N1            N2     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    i4mat_memory_test ( n_log );
  }
//
//  R4MAT test.
//
  cout << "\n";
  cout << "R4MAT Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N            N1            N2     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r4mat_memory_test ( n_log );
  }
//
//  R8MAT test.
//
  cout << "\n";
  cout << "R8MAT Memory Test\n";
  cout << "\n";
  cout << "Log2(N)            N            N1            N2     Ave       CPU        Real\n";
  cout << "\n";

  for ( n_log = n_log_min; n_log <= n_log_max; n_log++ )
  {
    r8mat_memory_test ( n_log );
  }

  cout << "\n";
  cout << "MEMORY_TEST\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void i4mat_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_MEMORY_TEST declares and uses an I4MAT of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int *i4mat;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  i4mat = new int[n1*n2];

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( float ) random ( ) / ( float ) RAND_MAX;
      i4mat[i+j*n1] = ( int ) ( 3.0 * x );
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + ( float ) i4mat[i+j*n1];
    }
  }
  average = average / ( float ) n1 / ( float ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(12) << n1
       << "  " << setw(12) << n2
       << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] i4mat;

  return;
}
//****************************************************************************80

void i4vec_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_MEMORY_TEST declares and uses an I4VEC of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int *i4vec;
  int n;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  i4vec = new int[n];

  for ( i = 0; i < n; i++ )
  {
    x = ( float ) random ( ) / ( float ) RAND_MAX;
    i4vec[i] = ( int ) ( 3.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + ( float ) i4vec[i];
  }
  average = average / ( float ) n;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] i4vec;

  return;
}
//****************************************************************************80

float r4_cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_CPU_TIME, the current total elapsed CPU time in second.
//
{
  float value;

  value = ( float ) clock ( ) / ( float ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

float r4_real_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_REAL_TIME returns the current real time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, float R4_REAL_TIME, the real time in seconds.
//
{
  time_t now;
  static time_t then = 0;
  float value;

  if ( then == 0 )
  {
    then = time ( NULL );
  }

  now = time ( NULL );

  value = ( float ) difftime ( now, then );

  return value;
}
//****************************************************************************80

void r4mat_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    R4MAT_MEMORY_TEST declares and uses an R4MAT of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  float *r4mat;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r4mat = new float[n1*n2];

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( float ) random ( ) / ( float ) RAND_MAX;
      r4mat[i+j*n1] = 2.0 * x;
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + r4mat[i+j*n1];
    }
  }
  average = average / ( float ) n1 / ( float ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(12) << n1
       << "  " << setw(12) << n2
       << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] r4mat;

  return;
}
//****************************************************************************80

void r4vec_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    R4VEC_MEMORY_TEST declares and uses an R4VEC of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  float average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int n;
  float *r4vec;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  float x;

  n = i4_power ( 2, n_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r4vec = new float[n];

  for ( i = 0; i < n; i++ )
  {
    x = ( float ) random ( ) / ( float ) RAND_MAX;
    r4vec[i] = ( 2.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + r4vec[i];
  }
  average = average / ( float ) n;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] r4vec;

  return;
}
//****************************************************************************80

double r8_cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CPU_TIME reports the elapsed CPU time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

double r8_real_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_REAL_TIME returns the current real time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_REAL_TIME, the real time in seconds.
//
{
  time_t now;
  static time_t then = 0;
  double value;

  if ( then == 0 )
  {
    then = time ( NULL );
  }

  now = time ( NULL );

  value = difftime ( now, then );

  return value;
}
//****************************************************************************80

void r8mat_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MEMORY_TEST declares and uses an R8MAT of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  double average;
  float cpu_diff;
  float cpu_time1;
  float cpu_time2;
  int i;
  int j;
  int n;
  int n1;
  int n1_log;
  int n2;
  int n2_log;
  double *r8mat;
  float real_diff;
  float real_time1;
  float real_time2;
  unsigned int seed = 123456789;
  double x;

  n = i4_power ( 2, n_log );

  n1_log = n_log / 2;
  n1 = i4_power ( 2, n1_log );
  n2_log = n_log - n1_log;
  n2 = i4_power ( 2, n2_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r4_real_time ( );
  cpu_time1 = r4_cpu_time ( );

  r8mat = new double[n1*n2];

  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      x = ( double ) random ( ) / ( double ) RAND_MAX;
      r8mat[i+j*n1] = 2.0 * x;
    }
  }

  average = 0.0;
  for ( j = 0; j < n2; j++ )
  {
    for ( i = 0; i < n1; i++ )
    {
      average = average + r8mat[i+j*n1];
    }
  }
  average = average / ( double ) n1 / ( double ) n2;

  cpu_time2 = r4_cpu_time ( );
  real_time2 = r4_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(12) << n1
       << "  " << setw(12) << n2
       << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] r8mat;

  return;
}
//****************************************************************************80

void r8vec_memory_test ( int n_log )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MEMORY_TEST declares and uses an R8VEC of size N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N_LOG, the logarithm base 2 of N.
//
{
  double average;
  double cpu_diff;
  double cpu_time1;
  double cpu_time2;
  int i;
  int n;
  double *r8vec;
  double real_diff;
  double real_time1;
  double real_time2;
  unsigned int seed = 123456789;
  double x;

  n = i4_power ( 2, n_log );

  cout << "  " << setw(4) << n_log
       << "  " << setw(12) << n;

  srandom ( seed );

  real_time1 = r8_real_time ( );
  cpu_time1 = r8_cpu_time ( );

  r8vec = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x = ( double ) random ( ) / ( double ) RAND_MAX;
    r8vec[i] = ( 2.0 * x );
  }

  average = 0.0;
  for ( i = 0; i < n; i++ )
  {
    average = average + r8vec[i];
  }
  average = average / ( double ) n;

  cpu_time2 = r8_cpu_time ( );
  real_time2 = r8_real_time ( );

  cpu_diff = cpu_time2 - cpu_time1;
  real_diff = real_time2 - real_time1;

  cout << "  " << setw(4) << setprecision(2) << average
       << "  " << setw(10) << setprecision(2) << cpu_diff
       << "  " << setw(10) << setprecision(2) << real_diff << "\n";

  delete [] r8vec;

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2003
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
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
