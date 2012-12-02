# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <omp.h>

using namespace std;

int main ( int argc, char *argv[] );
double f ( double x, double y );
double cpu_time ( );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUAD2D_OPENMP.
//
{
  double a;
  double b;
  double error;
  double exact;
  int i;
  int j;
  int n;
  int nx;
  int ny;
  double pi;
  double total;
  double wtime;
  double x;
  double y;

  a = 0.0;
  b = 1.0;
  nx = 32768;
  ny = 32768;
  n = nx * ny;
  pi = 3.141592653589793;
  exact = pi * pi / 6.0;

  cout << "\n";
  cout << "QUAD2D_OPENMP:\n";
  cout << "  C++/OpenMP version\n";
  cout << "  Estimate the integral of f(x,y) over [0,1]x[0,1].\n";
  cout << "  f(x,y) = 1 / ( 1 - x * y ).\n";
  cout << "\n";
  cout << "  A        = " << a << "\n";
  cout << "  B        = " << b << "\n";
  cout << "  NX       = " << nx << "\n";
  cout << "  NY       = " << ny << "\n";
  cout << "  N        = " << n << "\n";
  cout << "  Exact    = " << setprecision(16) << setw(24) << exact << "\n";

  wtime = omp_get_wtime ( );

  total = 0.0;
# pragma omp parallel shared ( a, b, n ) \
                      private ( i, j, x, y )

# pragma omp for reduction ( + : total )

  for ( i = 1; i <= nx; i++ )
  {
    x = ( ( 2 * nx - 2 * i + 1 ) * a + ( 2 * i - 1 ) * b ) / ( 2 * nx );
    for ( j = 1; j <= ny; j++ )
    {
      y = ( ( 2 * ny - 2 * j + 1 ) * a + ( 2 * j - 1 ) * b ) / ( 2 * ny );
      total = total + f ( x, y );
    }
  }

  wtime = omp_get_wtime ( ) - wtime;

  total = ( b - a ) * ( b - a ) * total / ( double ) ( nx ) / ( double ) ( ny );
  error = fabs ( total - exact );
 
  cout << "\n";
  cout << "  Estimate = " << setprecision(16) << setw(24) << total << "\n";
  cout << "  Error    = " << error << "\n";
  cout << "  Time     = " << wtime << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "QUAD_OPENMP:\n";
  cout << "  Normal end of execution.\n";

  return 0;
}
//****************************************************************************80

double f ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    F evaluates the function.
//
{
  double value;

  value = 1.0 / ( 1.0 - x * y );

  return value;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
// 
//    CPU_TIME reports the elapsed CPU time.
//
//  Modified:
//
//    27 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current total elapsed CPU time in second.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
