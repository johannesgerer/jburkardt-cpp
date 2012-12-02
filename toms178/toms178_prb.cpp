# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "toms178.hpp"

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
double rosenbrock ( double x[], int n );
double woods ( double x[], int n );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TOMS178_PRB.
//
//  Discussion:
//
//    TOMS178_PRB calls the TOMS178 routines.
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "TOMS178_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the TOMS178 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TOMS178_PRB:\n";
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
//    TEST01 tests HOOKE with the Rosenbrock function.
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *endpt;
  double eps;
  int i;
  int it;
  int itermax;
  int nvars = 2;
  double rho;
  double *startpt;
  double value;

  endpt = new double[nvars];
  startpt = new double[nvars];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  HOOKE seeks a minimizer of F(X).\n";
  cout << "  Here we use the Rosenbrock function.\n";
//
//  Starting guess for Rosenbrock.
//
  startpt[0] = -1.2;
  startpt[1] = 1.0;

  cout << "\n";
  cout << "  Initial estimate X =\n";
  cout << "\n";
  for ( i = 0; i < nvars; i++ )
  {
    cout << "  " << setw(8) << i + 1
         << "  " << setw(14) << startpt[i] << "\n";
  }

  value = rosenbrock ( startpt, nvars );

  cout << "\n";
  cout << "  F(X) = " << value << "\n";
//
//  Call HOOKE.
//
  itermax = 5000;
  rho = 0.5;
  eps = 1.0E-06;

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, &rosenbrock );
//
//  Results.
//
  cout << "\n";
  cout << "  Number of iterations taken = " << it << "\n";
  cout << "\n";
  cout << "  X* = \n";
  cout << "\n";
  for ( i = 0; i < nvars; i++ )
  {
    cout << "  " << setw(8) << i + 1
         << "  " << setw(14) << endpt[i] << "\n";
  }

  value = rosenbrock ( endpt, nvars );

  cout << "\n";
  cout << "  F(X*) = " << value << "\n";

  delete [] endpt;
  delete [] startpt;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests HOOKE with the WOODS function.
//
//  Discussion:
//
//    The Hooke and Jeeves algorithm works well when RHO = 0.5, but
//    does poorly when RHO = 0.6, and better when RHO = 0.8
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *endpt;
  double eps;
  int i;
  int it;
  int itermax;
  int nvars = 4;
  double rho;
  double *startpt;
  double value;

  endpt = new double[nvars];
  startpt = new double[nvars];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  HOOKE seeks a minimizer of F(X).\n";
  cout << "  Here we use the Rosenbrock function.\n";
  cout << "  Here we use the Woods function.\n";
//
//  Starting guess.
//
  startpt[0] = -3.0;
  startpt[1] = -1.0;
  startpt[2] = -3.0;
  startpt[3] = -1.0;

  cout << "\n";
  cout << "  Initial estimate X =\n";
  cout << "\n";
  for ( i = 0; i < nvars; i++ )
  {
    cout << "  " << setw(8) << i + 1
         << "  " << setw(14) << startpt[i] << "\n";
  }

  value = woods ( startpt, nvars );

  cout << "\n";
  cout << "  F(X) = " << value << "\n";
//
//  Call HOOKE.
//
  itermax = 5000;
  rho = 0.5;
  eps = 1.0E-06;

  it = hooke ( nvars, startpt, endpt, rho, eps, itermax, &woods );
//
//  Results.
//
  cout << "\n";
  cout << "  Number of iterations taken = " << it << "\n";
  cout << "\n";
  cout << "  X* = \n";
  cout << "\n";
  for ( i = 0; i < nvars; i++ )
  {
    cout << "  " << setw(8) << i + 1
         << "  " << setw(14) << endpt[i] << "\n";
  }

  value = woods ( endpt, nvars );

  cout << "\n";
  cout << "  F(X*) = " << value << "\n";

  delete [] endpt;
  delete [] startpt;
  return;
}
//****************************************************************************80

double rosenbrock ( double x[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK evaluates the Rosenbrock function.
//
//  Discussion:
//
//    The Hooke and Jeeves algorithm works reasonably well on
//    Rosenbrock's test function, depending on the value of RHO chosen.
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real X(N), the argument of the function.
//
//    Input, integer N, the spatial dimension.
//
//    Output, real ROSENBROCK, the value of the function.
//
{
  double value;

  value = 100.0 * pow ( x[1] - x[0] * x[0], 2 )
             +    pow ( 1.0 - x[0], 2 );

  return value;
}
//****************************************************************************80

double woods ( double x[], int n )

//****************************************************************************80
//
//  Purpose:
//
//    WOODS evaluates the Woods function.
//
//  Modified:
//
//    12 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, real X(N), the argument of the function.
//
//    Input, int N, the spatial dimension.
//
//    Output, real WOODS, the value of the function.
//
{
  double s1;
  double s2;
  double s3;
  double t1;
  double t2;
  double t3;
  double t4;
  double t5;
  double value;

  s1 = x[1] - x[0] * x[0];
  s2 = 1.0 - x[0];
  s3 = x[1] - 1.0;
  t1 = x[3] - x[2] * x[2];
  t2 = 1.0 - x[2];
  t3 = x[3] - 1.0;
  t4 = s3 + t3;
  t5 = s3 - t3;

  value = 100.0 * s1 * s1
        +         s2 * s2
        +  90.0 * t1 * t1
        +         t2 * t2
        +  10.0 * t4 * t4
        +   0.1 * t5 * t5;

  return value;
}

