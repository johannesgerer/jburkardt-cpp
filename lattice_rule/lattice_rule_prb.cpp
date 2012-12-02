# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "lattice_rule.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test085 ( );
void test09 ( );

void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LATTICE_RULE_PRB.
//
//  Discussion:
//
//    LATTICE_RULE_PRB calls the LATTICE_RULE test routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "LATTICE_RULE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LATTICE_RULE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test085 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LATTICE_RULE_PRB\n";
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
//    TEST01 tests FIBONACCI_LATTICE_Q.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  FIBONACCI_LATTICE_Q applies a Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_q ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests FIBONACCI_LATTICE_T.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  FIBONACCI_LATTICE_T applies a symmetric Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_t ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests FIBONACCI_LATTICE_B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  FIBONACCI_LATTICE_B applies an optimal Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_b ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests FIBONACCI_LATTICE_Q1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  FIBONACCI_LATTICE_Q1 applies a Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  A nonlinear coordinate transformation is applied.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_q1 ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests FIBONACCI_LATTICE_Q2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  FIBONACCI_LATTICE_Q2 applies a Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  A nonlinear coordinate transformation is applied.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_q2 ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests FIBONACCI_LATTICE_Q3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  FIBONACCI_LATTICE_Q3 applies a Fibonacci lattice rule\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "  A nonlinear coordinate transformation is applied.\n";
  cout << "  These Fibonacci rules are only available in 2D.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = fibonacci_lattice_q3 ( k, f_01_2d );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests LATTICE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int i;
  int m;
  double quad;
  int *z;;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  LATTICE applies a lattice rule to integrate\n";
  cout << "  a function over the unit hypercube.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  The lattice rule order M will vary.\n";

  z = new int[dim_num];

  z[0] = 1;
  z[1] = 2;

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

  i4vec_print ( dim_num, z, "  The lattice generator vector:" );

  cout << "\n";
  cout << "         I         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    m = prime ( 3 * i );

    quad = lattice ( dim_num, m, z, f_01_2d );

    exact = e_01_2d ( dim_num, a, b );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] z;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests LATTICE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int i;
  int m = 53;
  double quad;
  int *z;;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  LATTICE applies a lattice rule to integrate\n";
  cout << "  a function over the unit hypercube.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  The lattice rule order M will vary.\n";
  cout << "  The lattice generator vector Z will vary.\n";

  z = new int[dim_num];

  z[0] = 1;

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

  i4vec_print ( dim_num, z, "  The lattice generator vector:" );

  cout << "\n";
  cout << "         M      Z[0]      Z[1]      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( i = 1; i <= m - 1; i++ )
  {
    z[1] = i;

    quad = lattice ( dim_num, m, z, f_01_2d );

    exact = e_01_2d ( dim_num, a, b );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << m
         << "  " << setw(8) << z[0]
         << "  " << setw(8) << z[1]
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] z;

  return;
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests LATTICE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;
  int *z;;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  LATTICE is a lattice rule for periodic functions.\n";
  cout << "  However, we apply it to a nonperiodic function\n";
  cout << "  just to see how it does.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

  z = new int[dim_num];

  z[0] = 1;
  z[1] = 2;

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

  i4vec_print ( dim_num, z, "  The lattice generator vector:" );

  cout << "\n";
  cout << "         I         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = lattice ( dim_num, m, z, f_01_2d );

    exact = e_01_2d ( dim_num, a, b );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] z;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests LATTICE_NP0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;
  int *z;;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  LATTICE_NP0 applies a lattice rule to a\n";
  cout << "  nonperiodic function by reflecting the function\n";
  cout << "  about the midpoint and averaging.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

  z = new int[dim_num];

  z[0] = 1;
  z[1] = 2;

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

  i4vec_print ( dim_num, z, "  The lattice generator vector:" );

  cout << "\n";
  cout << "         I         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = lattice_np0 ( dim_num, m, z, f_01_2d );

    exact = e_01_2d ( dim_num, a, b );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] z;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests LATTICE_NP1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;
  int *z;;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  LATTICE_NP1 applies a lattice rule to a\n";
  cout << "  nonperiodic function using a nonlinear transformation\n";
  cout << "  to integrate a function over the unit square.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

  z = new int[dim_num];

  z[0] = 1;
  z[1] = 2;

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

  i4vec_print ( dim_num, z, "  The lattice generator vector:" );

  cout << "\n";
  cout << "         I         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 3; k <= 18; k++ )
  {
    m = fibonacci ( k );

    quad = lattice_np1 ( dim_num, m, z, f_01_2d );

    exact = e_01_2d ( dim_num, a, b );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;
  delete [] z;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests MONTE_CARLO.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 November 2008
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *b;
  int dim;
  int dim_num = 2;
  double error;
  double exact;
  int k;
  int m;
  double quad;
  int seed;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  MONTE_CARLO applies a Monte Carlo scheme\n";
  cout << "  to estimate the integral of a function\n";
  cout << "  over the unit hypercube.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

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
  seed = 123456789;

  exact = e_01_2d ( dim_num, a, b );

  cout << "\n";
  cout << "         K         M      EXACT     ESTIMATE  ERROR\n";
  cout << "\n";

  for ( k = 2; k <= 5; k++ )
  {
    m = i4_power ( 10, k );

    quad = monte_carlo ( dim_num, m, f_01_2d, &seed );

    error = r8_abs ( exact - quad );

    cout << "  " << setw(8) << k
         << "  " << setw(8) << m
         << "  " << setw(10) << exact
         << "  " << setw(10) << quad
         << "  " << setw(10) << error << "\n";
  }

  delete [] a;
  delete [] b;

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests LATTICE_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  int dim_num = 2;
  int m = 8;
  int *z;

  z = new int[dim_num];

  z[0] = 1;
  z[1] = 3;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  LATTICE_PRINT prints out the lattice generated\n";
  cout << "  by a single generator vector.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";

  i4vec_print ( dim_num, z, "  The generator vector:" );

  lattice_print ( dim_num, m, z, "  The total lattice:" );

  delete [] z;

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests FIND_Z20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    23 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  int dim;
  int dim_num = 2;
  int i;
  int m;
  int *z;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  FIND_Z20 finds the optimal lattice generator Z\n";
  cout << "  with Fourier coefficient smoothness ALPHA = 2,\n";
  cout << "'  and copy exponent 0,\n";
  cout << "  for a rank 1 \"method of good lattice points\" rule.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "\n";
  cout << "     M      Z(1)  Z(2)\n";
  cout << "\n";
  cout << "  (M = Fibonacci)\n";
  cout << "\n";

  for ( i = 3; i <= 10; i++ )
  {
    m = fibonacci(i);

    z = find_z20 ( dim_num, m );

    cout << "  " << setw(8) << m;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << z[dim];
    }
    cout << "\n";

    delete [] z;
  }

  cout << "\n";
  cout << "  (M = 2**K)\n";
  cout << "\n";

  for ( i = 2; i <= 10; i++ )
  {
    m = i4_power ( 2, i );

    z = find_z20 ( dim_num, m );

    cout << "  " << setw(8) << m;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << z[dim];
    }
    cout << "\n";

    delete [] z;
  }

  cout << "\n";
  cout << "  (M = 3*2**K)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    m = 3 * i4_power ( 2, i );

    z = find_z20 ( dim_num, m );

    cout << "  " << setw(8) << m;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << z[dim];
    }
    cout << "\n";

    delete [] z;
  }

  cout << "\n";
  cout << "  (M = Prime)\n";
  cout << "\n";

  for ( i = 3; i <= 10; i++ )
  {
    m = prime ( 10 * i );

    z = find_z20 ( dim_num, m );

    cout << "  " << setw(8) << m;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(8) << z[dim];
    }
    cout << "\n";

    delete [] z;
  }

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests FIBONACCI_LATTICE_Q_NODES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 January 2005.
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994, page 18.
//
{
  int dim_num = 2;
  int k;
  int m;
  double *x;

  k = 12;
  m = fibonacci ( k );

  cout << "\n";
  cout << "TEST14\n";
  cout << "  FIBONACCI_LATTICE_Q_NODES...\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  The Fibonacci index K =   " << k << "\n";
  cout << "  The Fibonacci value M =   " << m << "\n";

  x = fibonacci_lattice_q_nodes ( k );

  r8mat_transpose_print ( dim_num, m, x, "  The Fibonacci lattice nodes:" );

  delete [] x;

  return;
}
