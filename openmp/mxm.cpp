# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
void r8_mxm ( int l, int m, int n );
double r8_uniform_01 ( int *seed );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MXM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int id;
  int l;
  int m;
  int n;

  cout << "\n";
  cout << "MXM\n";
  cout << "  C++/OpenMP version.\n";
  cout << "\n";
  cout << "  Matrix multiplication tests.\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";

  l = 500;
  m = 500;
  n = 500;

  r8_mxm ( l, m, n );
//
//  Terminate.
//
  cout << "\n";
  cout << "MXM:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

void r8_mxm ( int l, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MXM carries out a matrix-matrix multiplication in R8 arithmetic.
//
//  Discussion:
//
//    A(LxN) = B(LxM) * C(MxN).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int L, M, N, the dimensions that specify the sizes of the
//    A, B, and C matrices.
//
{
  double *a;
  double *b;
  double *c;
  int i;
  int j;
  int k;
  int ops;
  double rate;
  int seed;
  double time_begin;
  double time_elapsed;
  double time_stop;
//
//  Allocate the matrices.
//
  a = new double [ l * n ];
  b = new double [ l * m ];
  c = new double [ m * n ];
//
//  Assign values to the B and C matrices.
//
  seed = 123456789;

  for ( k = 0; k < l * m; k++ )
  {
    b[k] = r8_uniform_01 ( &seed );
  }

  for ( k = 0; k < m * n; k++ )
  {
    c[k] = r8_uniform_01 ( &seed );
  }
//
//  Compute A = B * C.
//
  time_begin = omp_get_wtime ( );

# pragma omp parallel \
  shared ( a, b, c, l, m, n ) \
  private ( i, j, k )

# pragma omp for

  for ( j = 0; j < n; j++)
  {
    for ( i = 0; i < l; i++ )
    {
      a[i+j*l] = 0.0;
      for ( k = 0; k < m; k++ )
      {
        a[i+j*l] = a[i+j*l] + b[i+k*l] * c[k+j*m];
      }
    }
  }
  time_stop = omp_get_wtime ( );
//
//  Report.
//
  ops = l * n * ( 2 * m );
  time_elapsed = time_stop - time_begin;
  rate = ( double ) ( ops ) / time_elapsed / 1000000.0;

  cout << "\n";
  cout << "R8_MXM matrix multiplication timing.\n";
  cout << "  A(LxN) = B(LxM) * C(MxN).\n";
  cout << "  L = " << l  << "\n";
  cout << "  M = " << m  << "\n";
  cout << "  N = " << n  << "\n";
  cout << "  Floating point OPS roughly " << ops << "\n";
  cout << "  Elapsed time dT = " << time_elapsed  << "\n";
  cout << "  Rate = MegaOPS/dT = " << rate  << "\n";

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}

