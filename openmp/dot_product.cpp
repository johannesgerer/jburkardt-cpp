# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
double test01 ( int n, double x[], double y[] );
double test02 ( int n, double x[], double y[] );
double test03 ( int n, double x[], double y[] );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for DOT_PRODUCT.
//
//  Discussion:
//
//    This program illustrates how a vector dot product could be set up
//    in a FORTRAN90 program using OpenMP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2009
//
//  Author:
//
//    John Burkardt
//
{
  double factor;
  int i;
  int n;
  double wtime;
  double *x;
  double xdoty;
  double *y;

  cout << "\n";
  cout << "DOT_PRODUCT\n";
  cout << "  C++/OpenMP version\n";
  cout << "\n";
  cout << "  A program which computes a vector dot product.\n";

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";
//
//  Set up the vector data.
//  N may be increased to get better timing data.
//
//  The value FACTOR is chosen so that the correct value of the dot product 
//  of X and Y is N.
//
  n = 100;

  while ( n < 10000000 )
  {
    n = n * 10;

    x = new double[n];
    y = new double[n];

    factor = ( double ) ( n );
    factor = 1.0 / sqrt ( 2.0 * factor * factor + 3 * factor + 1.0 );

    for ( i = 0; i < n; i++ )
    {
      x[i] = ( i + 1 ) * factor;
    }

    for ( i = 0; i < n; i++ )
    {
      y[i] = ( i + 1 ) * 6 * factor;
    }

    cout << "\n";
//
//  Test #1
//
    wtime = omp_get_wtime ( );

    xdoty = test01 ( n, x, y );

    wtime = omp_get_wtime ( ) - wtime;

    cout << "  Sequential"
         << "  " << setw(8) << n
         << "  " << setw(14) << xdoty
         << "  " << setw(15) << wtime << "\n";
//
//  Test #2
//
    wtime = omp_get_wtime ( );

    xdoty = test02 ( n, x, y );

    wtime = omp_get_wtime ( ) - wtime;

    cout << "  Parallel  "
         << "  " << setw(8) << n
         << "  " << setw(14) << xdoty
         << "  " << setw(15) << wtime << "\n";

    delete [] x;
    delete [] y;
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "DOT_PRODUCT\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

double test01 ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 computes the dot product with no parallel processing directives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the vectors.
//
//    Input, double X[N], Y[N], the vectors.
//
//    Output, double TEST01, the dot product of X and Y.
//
{
  int i;
  double xdoty;

  xdoty = 0.0;

  for ( i = 0; i < n; i++ )
  {
    xdoty = xdoty + x[i] * y[i];
  }

  return xdoty;
}
//****************************************************************************80

double test02 ( int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 computes the dot product with parallel processing directives.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 April 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the vectors.
//
//    Input, double X[N], Y[N], the vectors.
//
//    Output, double TEST02, the dot product of X and Y.
//
{
  int i;
  double xdoty;

  xdoty = 0.0;

# pragma omp parallel \
  shared ( n, x, y ) \
  private ( i )

# pragma omp for reduction ( + : xdoty )

  for ( i = 0; i < n; i++ )
  {
    xdoty = xdoty + x[i] * y[i];
  }

  return xdoty;
}

