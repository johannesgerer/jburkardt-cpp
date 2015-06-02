# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "nintlib.hpp"

int main ( );

void testnd ( int dim_num, double func ( int dim_num, double x[] ) );
void test01 ( int dim_num, double func ( int dim_num, double x[] ) );
void test02 ( int dim_num, double func ( int dim_num, double x[] ) );
void test03 ( int dim_num, double func ( int dim_num, double x[] ) );
void test04 ( int dim_num, double func ( int dim_num, double x[] ) );
void test05 ( int dim_num, double func ( int dim_num, double x[] ) );
void test06 ( int dim_num, double func ( int dim_num, double x[] ) );
double f1dn ( int dim_num, double x[] );
double fbdn ( int dim_num, double x[] );
double fedn ( int dim_num, double x[] );
double fxdn ( int dim_num, double x[] );
double fx2dn ( int dim_num, double x[] );
double fx3dn ( int dim_num, double x[] );

//****************************************************************************8080

int main ( )

//****************************************************************************8080
//
//  Purpose:
//
//    MAIN is the main program for NINTLIB_PRB.
//
//  Discussion:
//
//    NINTLIB_PRB tests the NINTLIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 3

  double a;
  double b;
  int dim_num;
  int dim_num_test[TEST_NUM] = { 2, 3, 4 };
  int test;

  timestamp ( );
  cout << "\n";
  cout << "NINTLIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the NINTLIB library.\n";

  a = 0.0;
  b = 1.0;

  cout << "\n";
  cout << "TESTND\n";
  cout << "  Test routines for estimating the integral of\n";
  cout << "  of F(X) in the hypercube [A,B]**DIM_NUM.\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    dim_num = dim_num_test[test];

    cout << "\n";
    cout << "\n";
    cout << "  DIM_NUM = " << dim_num << "\n";
    cout << "\n";
    cout << "\n";
    cout << "  A(1:DIM_NUM) = " << a << "\n";
    cout << "  B(1:DIM_NUM) = " << b << "\n";

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = 1\n";
    cout << "\n";

    testnd ( dim_num, &f1dn );

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )\n";
    cout << "\n";

    testnd ( dim_num, &fxdn );

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = sum( X(1:DIM_NUM)^2 )\n";
    cout << "\n";

    testnd ( dim_num, &fx2dn );

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)^3 )\n";
    cout << "\n";

    testnd ( dim_num, &fx3dn );

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = exp(sum(X(1:DIM_NUM)))\n";
    cout << "\n";

    testnd ( dim_num, &fedn );

    cout << "\n";
    cout << "\n";
    cout << "  F(X(1:DIM_NUM)) = 1/(1+sum(X(1:DIM_NUM)^2))\n";
    cout << "\n";

    testnd ( dim_num, &fbdn );
  }
//
//  Terminate.
//
  cout << "\n";
  cout << "NINTLIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
# undef TEST_NUM
}
//****************************************************************************8080

void testnd ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************8080
//
//  Purpose:
//
//    TESTND tests the integrators on a particular function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
{
  test01 ( dim_num, func );
  test02 ( dim_num, func );
  test03 ( dim_num, func );
  test04 ( dim_num, func );
  if ( dim_num == 2 )
  {
    test05 ( dim_num, func );
  }
  test06 ( dim_num, func );

  return;
}
//****************************************************************************80

void test01 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BOX_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
# define ORDER 5

  int eval_num;
  int i;
  double result;
  double wtab[ORDER] = {
    0.236926885056189087514264040720,
    0.478628670499366468041291514836,
    0.568888888888888888888888888889,
    0.478628670499366468041291514836,
    0.236926885056189087514264040720 };
  double wtab2[ORDER];
  double xtab[ORDER] = {
    -0.906179845938663992797626878299,
    -0.538469310105683091036314420700,
     0.0,
     0.538469310105683091036314420700,
     0.906179845938663992797626878299 };
  double xtab2[ORDER];
//
//  Adjust the quadrature rule from [-1,1] to [0,1]:
//
  for ( i = 0; i < ORDER; i++ )
  {
    xtab2[i] = ( xtab[i] + 1.0 ) / 2.0;
  }
  for ( i = 0; i < ORDER; i++ )
  {
    wtab2[i] = 0.5 * wtab[i];
  }

  result = box_nd ( func, dim_num, ORDER, xtab2, wtab2, &eval_num );

  cout << "  BOX_ND:         "
       << setprecision(12) << setw(20) << result
       << setw(8) << "  " << eval_num << "\n";

  return;
# undef ORDER
}
//****************************************************************************80

void test02 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests P5_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
  double *a;
  double *b;
  int dim;
  int eval_num;
  double result;
//
//  Set the integration limits.
//
  a = new double[dim_num];
  b = new double[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    a[dim] = 0.0;
  }
  for ( dim = 0; dim < dim_num; dim++ )
  {
    b[dim] = 1.0;
  }

  result = p5_nd ( func, dim_num, a, b, &eval_num );

  cout << "  P5_ND:          "
         << setprecision(12) << setw(20) << result
         << setw(8) << "  " << eval_num << "\n";

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test03 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests ROMBERG_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
  double *a;
  double *b;
  int dim;
  int eval_num;
  int ind;
  int it_max = 3;
  double result;
  int *sub_num;
  double tol;
//
//  Set the integration limits.
//
  a = new double[dim_num];
  b = new double[dim_num];
  sub_num = new int[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    a[dim] = 0.0;
  }
  for ( dim = 0; dim < dim_num; dim++ )
  {
    b[dim] = 1.0;
  }
  tol = 0.001;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    sub_num[dim] = 10;
  }

  result = romberg_nd ( func, a, b, dim_num, sub_num, it_max, tol,
    &ind, &eval_num );

  cout << "  ROMBERG_ND:     "
         << setprecision(12) << setw(20) << result
         << setw(8) << "  " << eval_num << "\n";

  delete [] a;
  delete [] b;
  delete [] sub_num;

  return;
}
//****************************************************************************80

void test04 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests SAMPLE_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
# define K2 4

  double dev1[K2];
  double dev2[K2];
  double err1[K2];
  double est1[K2];
  double est2[K2];
  double err2[K2];
  int eval_num;
  int k1;

  k1 = 1;

  sample_nd ( func, k1, K2, dim_num, est1, err1, dev1, est2, err2,
    dev2, &eval_num );

  cout << "  SAMPLE_ND:      "
         << setprecision(12) << setw(20) << est2[K2-1]
         << setw(8) << "  " << eval_num << "\n";

  return;
# undef K2
}
//****************************************************************************80

void test05 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 demonstrates how to refine multi-dimensional integration results.
//
//  Discussion:
//
//    This routine is only set up for DIM_NUM = 2 for now.
//
//    We are given a routine, NDP5, which will integrate over a
//    DIM_NUM dimensional hypercube using a fixed method.  In order to
//    improve the approximation to an integral, we can subdivide
//    the hypercube and call NDP5 to integrate again over each of
//    these regions.
//
//    The information that we gather can be used to tell us when
//    to expect that we have achieved a certain degree of accuracy.
//
//    With a little more work, we could make this code adaptive.
//    That is, it would only refine SOME of the subregions, where
//    the approximation to the integral was still not good enough.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
  double *a;
  double *b;
  int dim;
  int eval_num;
  int eval_total;
  int i;
  int igrid;
  int j;
  int ngrid;
  double result;
  double result_total;
  double *xlo;
  double *xhi;

  a = new double[dim_num];
  b = new double[dim_num];
  xlo = new double[dim_num];
  xhi = new double[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    xlo[dim] = 0.0;
  }
  for ( dim = 0; dim < dim_num; dim++ )
  {
    xhi[dim] = 1.0;
  }

  for ( igrid = 1; igrid <= 6; igrid++ )
  {
    ngrid = i4_power ( 2, igrid - 1 );

    result_total = 0.0;
    eval_total = 0;

    for ( i = 1; i <= ngrid; i++ )
    {
      a[0] = ( ( double ) ( ngrid - i + 1 ) * xlo[0]
             + ( double ) (         i - 1 ) * xhi[0] )
             / ( double ) ( ngrid         );

      b[0] = ( ( double ) ( ngrid - i ) * xlo[0]
             + ( double ) (         i ) * xhi[0] )
             / ( double ) ( ngrid     );

      for ( j = 1; j <= ngrid; j++ )
      {
        a[1] = ( ( double ) ( ngrid - j + 1 ) * xlo[1]
               + ( double ) (         j - 1 ) * xhi[1] )
               / ( double ) ( ngrid         );

        b[1] = ( ( double ) ( ngrid - j ) * xlo[1]
               + ( double ) (         j ) * xhi[1] )
               / ( double ) ( ngrid     );

        result = p5_nd ( func, dim_num, a, b, &eval_num );

        result_total = result_total + result;
        eval_total = eval_total + eval_num;
      }
    }
    cout << "  P5_ND+:         "
         << setprecision(12) << setw(20) << result_total
         << setw(8) << "  " << eval_total << "\n";

  }
  delete [] a;
  delete [] b;
  delete [] xhi;
  delete [] xlo;

  return;
}
//****************************************************************************80

void test06 ( int dim_num, double func ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests MONTE_CARLO_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates the function
//    to be integrated.
//
{
  double *a;
  double *b;
  int dim;
  int eval_num;
  double result;
  int seed;
  int test;
  int test_num = 3;

  seed = 123456789;
//
//  Set the integration limits.
//
  a = new double[dim_num];
  b = new double[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    a[dim] = 0.0;
  }
  for ( dim = 0; dim < dim_num; dim++ )
  {
    b[dim] = 1.0;
  }

  for ( test = 1; test <= test_num; test++ )
  {
    eval_num = i4_power ( 8, test ) * 10000;

    result = monte_carlo_nd ( func, dim_num, a, b, eval_num, &seed );

    cout << "  MONTE_CARLO_ND: "
         << setprecision(12) << setw(20) << result
         << setw(8) << "  " << eval_num << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

double fbdn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FBDN(X(1:DIM_NUM)) = 1 / ( 1 + sum ( X(1:DIM_NUM)**2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double FBDN, the value of the function at X.
//
{
  double arg;
  int dim;
  double value;

  arg = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    arg = arg + x[dim] * x[dim];
  }

  value = 1.0 / ( 1.0 + arg );

  return value;
}
//****************************************************************************80

double fedn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FEDN(X(1:DIM_NUM)) = EXP ( sum ( X(1:DIM_NUM) ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double FEDN, the value of the function at X.
//
{
  double arg;
  int dim;
  double value;

  arg = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    arg = arg + x[dim];
  }

  value = exp ( arg );

  return value;
}
//****************************************************************************80

double f1dn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F1DN(X(1:DIM_NUM)) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double F1DN, the value of the function at X.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double fxdn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FXDN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double FXDN, the value of the function at X.
//
{
  double arg;
  int dim;
  double value;

  arg = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    arg = arg + x[dim];
  }

  value = arg;

  return value;
}
//****************************************************************************80

double fx2dn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX2DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double FX2DN, the value of the function at X.
//
{
  double arg;
  int dim;
  double value;

  arg = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    arg = arg + x[dim] * x[dim];
  }

  value = arg;

  return value;
}
//****************************************************************************80

double fx3dn ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    FX3DN(X(1:DIM_NUM)) = sum ( X(1:DIM_NUM)**3 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the argument.
//
//    Output, double FX3DN, the value of the function at X.
//
{
  double arg;
  int dim;
  double value;

  arg = 0.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    arg = arg + x[dim] * x[dim] * x[dim];
  }

  value = arg;

  return value;
}
