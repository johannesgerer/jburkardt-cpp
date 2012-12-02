# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

int main ( int argc, char *argv[] );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double r8_sum ( double x, double y );
double r8_uniform_01 ( int *seed );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for the FUNCTION POINTER demonstration.
//
//  Discussion:
//
//    This is an example of how to declare and use a function pointer.
//    Essentially, this is a symbolic name that can be associated with
//    a function; the association can be changed to a different function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double ( *func ) ( double x, double y );
  int i;
  int seed;
  
  timestamp ( );
  cout << "\n";
  cout << "FUNCTION_POINTER\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Examples of function pointers.\n";
  cout << "  We define a variable FUNC which can point to a function.\n";
  cout << "  Our declaration indicates only that FUNC has two inputs\n";
  cout << "  of type double, and returns a double as its value.\n";

  cout << "\n";
  cout << "  We will now set FUNC to a MIN, a MAX and a SUM function\n";
  cout << "  successively, and invoke it with 5 random sets of input.\n";

  seed = 12345689;
  func = r8_min;

  cout << "\n";
  cout << "  FUNC = R8_MIN\n";
  cout << "  -----A----  -----B----  C = FUNC ( A, B )\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    a = r8_uniform_01 ( &seed );
    b = r8_uniform_01 ( &seed );
    c = func ( a, b );
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(10) << c << "\n";
  }
  
  seed = 12345689;
  func = r8_max;

  cout << "\n";
  cout << "  FUNC = R8_MAX\n";
  cout << "  -----A----  -----B----  C = FUNC ( A, B )\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    a = r8_uniform_01 ( &seed );
    b = r8_uniform_01 ( &seed );
    c = func ( a, b );
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(10) << c << "\n";
  }

  seed = 123456789;
  func = r8_sum;

  cout << "\n";
  cout << "  FUNC = R8_SUM\n";
  cout << "  -----A----  -----B----  C = FUNC ( A, B )\n";
  cout << "\n";

  for ( i = 0; i < 5; i++ )
  {
    a = r8_uniform_01 ( &seed );
    b = r8_uniform_01 ( &seed );
    c = func ( a, b );
    cout << "  " << setw(10) << a
         << "  " << setw(10) << b
         << "  " << setw(10) << c << "\n";
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "FUNCTION_PONTER:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  } 
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = y;
  } 
  else
  {
    value = x;
  }
  return value;
}
//****************************************************************************80

double r8_sum ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_SUM returns the sum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_SUM, the sum of X and Y.
//
{
  double value;

  value = x + y;

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number.
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
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
