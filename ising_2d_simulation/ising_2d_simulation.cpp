# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <string>

using namespace std;

int main ( int argc, char *argv[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
void ising_2d_agree ( int m, int n, int c1[], int c5[] );
int *ising_2d_initialize ( int m, int n, double thresh, int *seed );
void ising_2d_stats ( int step, int m, int n, int c1[] );
void neighbor_2d_stats ( int step, int m, int n, int c1[], int c5[] );
void r8mat_uniform_01 ( int m, int n, int *seed, double r[] );
void timestamp ( );
void transition ( int m, int n, int iterations, double prob[], 
  double thresh, int *seed );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ISING_2D_SIMULATION.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2011
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int iterations;
  int m;
  int n;
  double prob[5] = { 0.98, 0.85, 0.50, 0.15, 0.02 };
  int seed;
  int step;
  double thresh;

  timestamp ( );
  cout << "\n";
  cout << "ISING_2D_SIMULATION\n";
  cout << "  C++ version\n";
  cout << "  Monte Carlo simulation of a 2D Ising model.\n";
//
//  Get input.
//
  if ( 1 < argc )
  {
    m = atoi ( argv[1] );
  }
  else
  {
    m = 10;
  }
  if ( 2 < argc )
  {
    n = atoi ( argv[2] );
  }
  else
  {
    n = 10;
  }
  if ( 3 < argc )
  {
    iterations = atoi ( argv[3] );
  }
  else
  {
    iterations = 15;
  }
  if ( 4 < argc )
  {
    thresh = atof ( argv[4] );
  }
  else
  {
    thresh = 0.50;
  }
  if ( 5 < argc )
  {
    seed = atoi ( argv[5] );
  }
  else
  {
    seed = 123456789;
  }
  cout << "\n";
  cout << "  The number of rows is M = " << m << "\n";
  cout << "  The number of columns is N = " << n << "\n";
  cout << "  The number of iterations taken is ITERATIONS = " << iterations << "\n";
  cout << "  The threshhold THRESH = " << thresh << "\n";
  cout << "  The seed SEED = " << seed << "\n";
  cout << "\n";
  cout << "  The transition probability table, based on the number of\n";
  cout << "  neighbors with the same spin.\n";
  cout << "\n";
  cout << "      1         2         3         4         5\n";
  cout << "\n";
  for ( i = 0; i < 5; i++ )
  {
    cout << setw[10] <<  prob[i];
  }
  cout << "\n";
//
//  Do the simulation.
//
  transition ( m, n, iterations, prob, thresh, &seed );
//
//  Terminate.
//
  cout << "\n";
  cout << "ISING_2D_SIMULATION\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
{
  int value;

  if ( i1 < i2 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I   Value
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

void ising_2d_agree ( int m, int n, int c1[], int c5[] )

//****************************************************************************80
//
//  Purpose:
//
//    ISING_2D_AGREE returns the number of neighbors agreeing with each cell.
//
//  Discussion:
//
//    The count includes the cell itself, so it is between 1 and 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 Noveber 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of cells in each 
//    spatial dimension.
//
//    Input, int C1[M*N], an array of 1's and -1's.
//
//    Output, int C5[M*N], the number of neighbors 
//    that agree.  1, 2, 3, 4, or 5.
//
{
  int i;
  int im;
  int ip;
  int j;
  int jm;
  int jp;

  for ( j = 0; j < n; j++ )
  {
    jp = i4_wrap ( j + 1, 0, n - 1 );
    jm = i4_wrap ( j - 1, 0, n - 1 );
    for ( i = 0; i < m; i++ )
    {
      ip = i4_wrap ( i + 1, 0, m - 1 );
      im = i4_wrap ( i - 1, 0, m - 1 );
      c5[i+j*m] = c1[i+j*m] + c1[ip+j*m] + c1[im+j*m] + c1[i+jm*m] + c1[i+jp*m];
      if ( 0 < c1[i+j*m] )
      {
        c5[i+j*m] = ( 5 + c5[i+j*m] ) / 2;
      }
      else
      {
        c5[i+j*m] = ( 5 - c5[i+j*m] ) / 2;
      }
    }
  }
  return;
}
//****************************************************************************80

int *ising_2d_initialize ( int m, int n, double thresh, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    ISING_2D_INITIALIZE initializes the Ising array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double THRESH, the threshhold.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
//    Output, in ISING_2D_INITIALIZE[M*N], the initial Ising array.
//
{
  int *c1;
  int i;
  int j;
  double *r;

  r = new double[m*n];

  r8mat_uniform_01 ( m, n, seed, r );

  c1 = new int[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( r[i+j*m] <= thresh )
      {
        c1[i+j*m] = -1;
      }
      else
      {
        c1[i+j*m] = +1;
      }
    }
  }
  delete [] r;

  return c1;
}
//****************************************************************************80

void ising_2d_stats ( int step, int m, int n, int c1[] )

//****************************************************************************80
//
//  Purpose:
//
//    ISING_2D_STATS prints information about the current step.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the step number.
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int C1[M*N], the current state of the system.
//
{
  int i;
  int j;
  int pos_count;
  double pos_percent;
  int neg_count;
  double neg_percent;

  if ( step == 0 )
  {
    cout << "\n";
    cout << "  Step     Positives       Negatives\n";
    cout << "             #    %%          #    %%\n";
    cout << "\n";
  }

  pos_count = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( 0 < c1[i+j*m] )
      {
        pos_count = pos_count + 1;
      }
    }
  }
  neg_count = m * n - pos_count;
  pos_percent = ( double ) ( 100 * pos_count ) / ( double ) ( m * n );
  neg_percent = ( double ) ( 100 * neg_count ) / ( double ) ( m * n );

  cout << "  " << setw(4) << step
       << "  " << setw(6) << pos_count
       << "  " << setw(6) << pos_percent
       << "  " << setw(6) << neg_count
       << "  " << setw(6) << neg_percent << "\n";

  return;
}
//****************************************************************************80

void neighbor_2d_stats ( int step, int m, int n, int c1[], int c5[] )

//****************************************************************************80
/*
  Purpose:

    NEIGHBOR_2D_STATS prints neighbor statistics about the current step.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    23 November 2011

  Author:

    John Burkardt

  Parameters:

    Input, int STEP, the step number.

    Input, int M, N, the number of rows and columns.

    Input, int C1[M*N], the current state of the system.

    Input, int C5[M*N], the number of agreeable neighbors.
*/
{
  int i;
  int j;
  int k;
  int stats[11];

  if ( step == 0 )
  {
    cout << "\n";
    cout << "  Step     Neighborhood Charge:\n";
    cout << "           -5    -4    -3    -2    -1     +1    +2    +3    +4    +5\n";
    cout << "\n";
  }

  for ( i = - 5; i <= 5; i++ )
  {
    stats[i+5] = 0;
  }

  for (j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      stats[c5[i+j*m]-1+5] = stats[c5[i+j*m]-1+5] + 1;
    }
  }
  cout << "  " << setw(4) << step;
  for ( i = - 5; i <= 5; i++ )
  {
    if ( i != 0 )
    {
      cout << "  " << setw(4) << stats[i+5];
    }
  }
  cout << "\n";
      
  return;
}
//****************************************************************************80

void r8mat_uniform_01 ( int m, int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is an array of R8's.
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
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
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has
//    been updated.
//
//    Output, double R[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int j;
  int k;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8MAT_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + i4_huge;
      }

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }
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
//****************************************************************************80

void transition ( int m, int n, int iterations, double prob[], 
  double thresh, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, int ITERATIONS, the number of iterations.
//
//    Input, double PROB[5].  PROB[I-1] represents the probability 
//    that the spin of a given cell will be reversed, given that it has I 
//    immediate neighbors (including itself) with spin the same as its own.
//
//    Input, double THRESH, the threshhold.
//
//    Input/output, int *SEED, a seed for the random number 
//    generator.
//
{
  int *c1;
  int *c5;
  int i;
  int j;
  double *r;
  int step;

  c5 = new int[m*n];

  r = new double[m*n];

  c1 = ising_2d_initialize ( m, n, thresh, seed );

  step = 0;
  ising_2d_stats ( step, m, n, c1 );

  for ( step = 1; step <= iterations; step++ )
  {
//
//  C5 contains 1 through 5, the number of cells that agree with the center cell.
//
    ising_2d_agree ( m, n, c1, c5 );

    if ( 0 )
    {
      neighbor_2d_stats ( step, m, n, c1, c5 );
    }
//
//  Determine the chances of flipping cell (I,J).
//
    r8mat_uniform_01 ( m, n, seed, r );

    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        if ( r[i+j*m] < prob[c5[i+j*m]-1] )
        {
          c1[i+j*m] = - c1[i+j*m];
        }
      }
    }
    ising_2d_stats ( step, m, n, c1 );
  }
  delete [] c1;
  delete [] c5;
  delete [] r;

  return;
}
