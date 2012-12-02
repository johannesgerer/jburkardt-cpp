# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "test_int_2d.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TEST_INT_2D_PRB.
//
//  Discussion:
//
//    TEST_INT_2D_PRB demonstrates the TEST_INT_2D integration test functions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TEST_INT_2D_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the TEST_INT_2D library.\n";
 
  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TEST_INT_2D_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 applies a Monte Carlo rule.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double *fx;
  int i;
  int j;
  int n;
  int problem;
  int problem_num;
  double quad;
  int seed;
  double volume;
  double *x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use a Monte Carlo rule.\n";
  cout << "\n";
  cout << "  Repeatedly multiply the number of points by 4.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "   Problem      Points         Approx            Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    cout << "\n";
    n = 1;
    for ( i = 1; i <= 12; i++ )
    {
      seed = 123456789;

      x = new double[2*n];
      fx = new double[n];

      x = r8mat_uniform_01 ( 2, n, &seed );

      p00_lim ( problem, a, b );

      for ( dim = 0; dim < 2; dim++ )
      {
        for ( j = 0; j < n; j++ )
        {
          x[dim+j*2] = ( 1.0 - x[dim+j*2] ) * a[dim]
                     +         x[dim+j*2]   * b[dim];
        }
      }
      volume = ( b[1] - a[1] ) * ( b[0] - a[0] );

      p00_fun ( problem, n, x, fx );

      quad = volume * r8vec_sum ( n, fx ) / ( double ) ( n );
   

      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      cout << "  " << setw(8)  << problem
           << "  " << setw(10) << n
           << "  " << setw(14) << quad
           << "  " << setw(14) << error << "\n";

      delete [] fx;
      delete [] x;

      n = n * 4;
    }
    cout << "  " << setw(8)  << problem
         << "       Exact"
         << "  " << setw(14) << exact << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 applies a product of composite midpoint rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double *fx;
  int i;
  int ix;
  int iy;
  int k;
  int n;
  int nx;
  int ny;
  int problem;
  int problem_num;
  double quad;
  double volume;
  double *x;
  double xval;
  double yval;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use a product of composite midpoint rules..\n";
  cout << "  Repeatedly multiply the number of points by 4.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "   Problem      Points         Approx            Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    cout << "\n";
    nx = 1;
    ny = 1;

    for ( i = 1; i <= 12; i++ )
    {
      n = nx * ny;

      x = new double[2*n];
      fx = new double[n];

      p00_lim ( problem, a, b );

      k = 0;

      for ( ix = 1; ix <= nx; ix++ )
      {
        xval = ( ( double ) ( 2 * nx - 2 * ix + 1 ) * a[0]   
               + ( double ) (          2 * ix - 1 ) * b[0] ) 
               / ( double ) ( 2 * nx              );

        for ( iy = 1; iy <= ny; iy++ )
        {
          yval = ( ( double ) ( 2 * ny - 2 * iy + 1 ) * a[1]
                 + ( double ) (          2 * iy - 1 ) * b[1] ) 
                 / ( double ) ( 2 * ny              );

          x[0+k*2] = xval;
          x[1+k*2] = yval;
          k = k + 1;
        }
      }

      volume = ( b[1] - a[1] ) * ( b[0] - a[0] );

      p00_fun ( problem, n, x, fx );

      quad = volume * r8vec_sum ( n, fx ) / ( double ) ( n );
   
      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      cout << "  " << setw(8)  << problem
           << "  " << setw(10) << n
           << "  " << setw(14) << quad
           << "  " << setw(14) << error << "\n";

      delete [] fx;
      delete [] x;

      nx = nx * 2;
      ny = ny * 2;
    }
    cout << "  " << setw(8)  << problem
         << "       Exact"
         << "  " << setw(14) << exact << "\n";
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 applies a product of Gauss-Legendre rules.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2011
//
//  Author:
//
//    John Burkardt
//
{
  double a[2];
  double b[2];
  int dim;
  double error;
  double exact;
  double  *fxy;
  int i;
  int ix;
  int iy;
  int j;
  int k;
  int nx;
  int nxy;
  int ny;
  int problem;
  int problem_num;
  double quad;
  double volume;
  double *w;
  double *wxy;
  double *x;
  double *xy;
  double xval;
  double yval;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Use a product of Gauss-Legendre rules.\n";
  cout << "  The 1D rules essentially double in order.\n";

  problem_num = p00_problem_num ( );

  cout << "\n";
  cout << "   Problem      Points       Approx         Error\n";

  for ( problem = 1; problem <= problem_num; problem++ )
  {
    cout << "\n";

    nx = 1;
    ny = 1;

    for ( i = 1; i <= 8; i++ )
    {
      x = new double[nx];
      w = new double[nx];

      legendre_dr_compute ( nx, x, w );

      nxy = nx * ny;

      wxy = new double[nxy];
      xy = new double[2*nxy];
      fxy = new double[nxy];

      p00_lim ( problem, a, b );

      k = 0;

      for ( ix = 0; ix < nx; ix++ )
      {
        xval = ( ( 1.0 + x[ix] ) * a[0]   
               + ( 1.0 - x[ix] ) * b[0] ) 
               /   2.0;
        for ( iy = 0; iy < ny; iy++ )
        {
          yval = ( ( 1.0 + x[iy] ) * a[1]   
                 + ( 1.0 - x[iy] ) * b[1] ) 
                 /   2.0;
          xy[0+k*2] = xval;
          xy[1+k*2] = yval;
          wxy[k] = w[ix] * w[iy];
          k = k + 1;
        }
      }
      volume = ( b[0] - a[0] ) * ( b[1] - a[1] );

      p00_fun ( problem, nxy, xy, fxy );

      quad = 0.0;
      for ( j = 0; j < nxy; j++ )
      {
        quad = quad + wxy[j] * fxy[j];
      }
      quad = quad * volume / 4.0;

      exact = p00_exact ( problem );

      error = r8_abs ( quad - exact );

      cout << "  " << setw(8) << problem
           << "  " << setw(10) << nxy
           << "  " << setw(14) << quad
           << "  " << setw(14) << error << "\n";

      delete [] fxy;
      delete [] w;
      delete [] wxy;
      delete [] x;
      delete [] xy;

      nx = 2 * nx + 1;
      ny = nx;
    }
    cout << "  " << setw(8) << problem
         << "  " << "     Exact"
         << "  " << setw(14) << exact << "\n";
  }
  return;
}
