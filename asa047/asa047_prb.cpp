# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

# include "asa047.hpp"

using namespace std;

int main ( );
void test01 ( );
double rosenbrock ( double x[2] );
void test02 ( );
double powell ( double x[4] );
void test03 ( );
double helical ( double x[3] );
void test04 ( );
double quartic ( double x[10] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA047_PRB.
//
//  Discussion:
//
//    ASA047_PRB tests the ASA047 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA047_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA047 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA047_PRB:\n";
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
//    TEST01 demonstrates the use of NELMIN on ROSENBROCK.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = 2;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Apply NELMIN to ROSENBROCK function.\n";

  start[0] = -1.2;
  start[1] =  1.0;

  reqmin = 1.0E-08;

  step[0] = 1.0;
  step[1] = 1.0;

  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = rosenbrock ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( rosenbrock, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************80

double rosenbrock ( double x[2] )

//****************************************************************************80
//
//  Purpose:
//
//    ROSENBROCK evaluates the Rosenbrock parabolic value function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[2], the argument.
//
//    Output, double ROSENBROCK, the value of the function.
//
{
  double fx;
  double fx1;
  double fx2;

  fx1 = x[1] - x[0] * x[0];
  fx2 = 1.0 - x[0];

  fx = 100.0 * fx1 * fx1
     +         fx2 * fx2;

  return fx;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 demonstrates the use of NELMIN on POWELL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = 4;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Apply NELMIN to POWELL quartic function.\n";

  start[0] =   3.0;
  start[1] = - 1.0;
  start[2] =   0.0;
  start[3] =   1.0;

  reqmin = 1.0E-08;

  step[0] = 1.0;
  step[1] = 1.0;
  step[2] = 1.0;
  step[3] = 1.0;

  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = powell ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( powell, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************80

double powell ( double x[4] )

//****************************************************************************80
//
//  Purpose:
//
//    POWELL evaluates the Powell quartic function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[4], the argument.
//
//    Output, double POWELL, the value of the function.
//
{
  double fx;
  double fx1;
  double fx2;
  double fx3;
  double fx4;

  fx1 = x[0] + 10.0 * x[1];
  fx2 = x[2] - x[3];
  fx3 = x[1] - 2.0 * x[2];
  fx4 = x[0] - x[3];

  fx =            fx1 * fx1
     +  5.0 * fx2 * fx2
     +            fx3 * fx3 * fx3 * fx3
     + 10.0 * fx4 * fx4 * fx4 * fx4;;

  return fx;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates the use of NELMIN on HELICAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = 3;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Apply NELMIN to the HELICAL function.\n";

  start[0] = - 1.0;
  start[1] =   0.0;
  start[2] =   0.0;

  reqmin = 1.0E-08;

  step[0] = 1.0;
  step[1] = 1.0;
  step[2] = 1.0;

  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = helical ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( helical, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************80

double helical ( double x[3] )

//****************************************************************************80
//
//  Purpose:
//
//    HELICAL evaluates the Fletcher-Powell helical valley function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[3], the argument.
//
//    Output, double HELICAL, the value of the function.
//
{
  double fx;
  double fx1;
  double fx2;
  double fx3;
  double pi = 3.141592653589793;
  double theta;

  if ( 0.0 < x[0] )
  {
    theta = atan2 ( x[1], x[0] ) / 2.0 / pi;
  }
  else if ( x[0] < 0.0 )
  {
    theta = 0.5 + atan2 ( x[1], x[0] ) / 2.0 / pi;
  }
  else if ( x[0] == 0.0 )
  {
    theta = 0.25;
  }

  fx1 = x[2] - 10.0 * theta;
  fx2 = sqrt ( x[0] * x[0] + x[1] * x[1] );
  fx3 = x[2];

  fx = 100.0 * fx1 * fx1
     +         fx2 * fx2
     +         fx3 * fx3;

  return fx;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrates the use of NELMIN on QUARTIC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 February 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int icount;
  int ifault;
  int kcount;
  int konvge;
  int n;
  int numres;
  double reqmin;
  double *start;
  double *step;
  double *xmin;
  double ynewlo;

  n = 10;

  start = new double[n];
  step = new double[n];
  xmin = new double[n];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Apply NELMIN to the QUARTIC function.\n";

  for ( i = 0; i < n; i++ )
  {
    start[i] = 1.0;
  }

  reqmin = 1.0E-08;

  for ( i = 0; i < n; i++ )
  {
    step[i] = 1.0;
  }

  konvge = 10;
  kcount = 500;

  cout << "\n";
  cout << "  Starting point X:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << start[i] << "\n";
  }

  ynewlo = quartic ( start );

  cout << "\n";
  cout << "  F(X) = " << ynewlo << "\n";

  nelmin ( quartic, n, start, xmin, &ynewlo, reqmin, step,
    konvge, kcount, &icount, &numres, &ifault );

  cout << "\n";
  cout << "  Return code IFAULT = " << ifault << "\n";
  cout << "\n";
  cout << "  Estimate of minimizing value X*:\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(14) << xmin[i] << "\n";
  }

  cout << "\n";
  cout << "  F(X*) = " << ynewlo << "\n";

  cout << "\n";
  cout << "  Number of iterations = " << icount << "\n";
  cout << "  Number of restarts =   " << numres << "\n";

  delete [] start;
  delete [] step;
  delete [] xmin;

  return;
}
//****************************************************************************80

double quartic ( double x[10] )

//****************************************************************************80
//
//  Purpose:
//
//    QUARTIC evaluates a function defined by a sum of fourth powers.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double X[10], the argument.
//
//    Output, double QUARTIC, the value of the function.
//
{
  double fx;
  int i;

  fx = 0.0;

  for ( i = 0; i < 10; i++ )
  {
    fx = fx + x[i] * x[i] * x[i] * x[i];
  }

  return fx;
}
