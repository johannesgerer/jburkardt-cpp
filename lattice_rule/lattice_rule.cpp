# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "lattice_rule.hpp"

//****************************************************************************80

double e_01_2d ( int dim_num, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    E_01_2D is the exact integral of 2d test function #1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
//
//    Output, double E_01_2D, the integral of the function 
//    over the limits.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

double f_01_2d ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F_01_2D is the 2D test function #1.
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the point where the function 
//    is to be evaluated.
//
//    Output, double F_01_2D, the value of the function at X.
//
{
  double e = 2.718281828459045;
  double value;

  value = x[1] * exp ( x[0] * x[1] ) / ( e - 2.0 );

  return value;
}
//****************************************************************************80

double f2 ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    F2 evaluates a function of a scalar used in defining P2(Q).
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, double X, the value of the argument.
//
//    Output, double F2, the value of F2(X).
//
{
  double pi = 3.141592653589793;
  double value;

  value = 1.0 + 2.0 * pi * pi * ( x * x - x + 1.0 / 6.0 );

  return value;
}
//****************************************************************************80

double f20_s ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F20_S evaluates a function of a vector used in defining P2(Q).
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the value of the argument.
//
//    Output, double F20_S, the value of F20_S(X).
//
{
  int i;
  double value;

  value = 1.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value * ( 1.0 + ( f2 ( x[i] ) - 1.0 ) );
  }
  value = value - 1.0;

  return value;
}
//****************************************************************************80

int fibonacci ( int k )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI returns the Fibonacci number of given index.
//
//  Example:
//
//    K   Fibonacci
//
//    0   0
//    1   1
//    2   1
//    3   2
//    4   3
//    5   5
//    6   8
//    7  13
//    8  21
//    9  34
//   10  55
//   11  89
//   12 144
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 1.
//
//    Output, int FIBONACCI, the value of the K-th Fibonacci number.
//
{
  int a;
  int b;
  int c;
  int kk;

  if ( k < 0 )
  {
    a = - i4_huge ( );
    return a;
  }
  else if ( k == 0 )
  {
    a = 0;
    return a;
  }
  else if ( k == 1 )
  {
    a = 1;
    return a;
  }

  c = 0;
  b = 0;
  a = 1;

  for ( kk = 2; kk <= k; kk++ )
  {
    c = b;
    b = a;
    a = c + b;
  }
  return a;
}
//****************************************************************************80

double fibonacci_lattice_b ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_B applies an optimal Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    This routine may be applied to integrands which are not periodic.
//
//    When K is odd, this is the same as the symmetric Fibonacci lattice
//    integration rule.  But when K is even, a correction is made to the
//    corner weights which is expected to improve the results.
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_B, the estimated integral.
//
{
  double delta;
  int dim;
  int dim_num = 2;
  int j;
  int m;
  int n;
  double quad;
  double quad1;
  double quad2;
  int rank;
  double w[2*2];
  double *x;
  int *z;

  x = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );
  n = fibonacci ( k - 1 );
//
//  Get the corner weights.
//
  if ( ( k % 2 ) == 1 )
  {
    w[0+0*2] = 1.0 / ( double ) ( 4 * m );
    w[1+0*2] = 1.0 / ( double ) ( 4 * m );
    w[0+1*2] = 1.0 / ( double ) ( 4 * m );
    w[1+1*2] = 1.0 / ( double ) ( 4 * m );
  }
  else
  {
     delta = 0.0;
     for ( j = 1; j <= m - 1; j++ )
     {
       delta = delta + ( double ) ( j * ( ( j * n ) % m ) ) 
                     / ( double ) ( m * m );
     }
     w[0+0*2] = 0.25 - delta / ( double ) ( m );

     delta = 0.0;
     for ( j = 1; j <= m - 1; j++ )
     {
       delta = delta + ( double ) ( j * ( m - ( ( j * n ) % m ) ) ) 
         / ( double ) ( m * m );
     }
     w[0+1*2] = 0.25 - delta / ( double ) ( m );

     w[1+0*2] = w[0+1*2];
     w[1+1*2] = w[0+0*2];
  }
//
//  Get all the corner values.
//
  rank = 0;
  quad1 = 0.0;

  for ( ; ; )
  {
    tuple_next ( 0, 1, dim_num, &rank, z );

    if ( rank == 0 )
    {
      break;
    }
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = ( double ) z[dim];
    }
    quad1 = quad1 + w[z[0]+z[1]*2] * f ( dim_num, x );
  }
//
//  Get the interior values.
//
  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  quad2 = 0.0;
  for ( j = 1; j <= m - 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = fmod ( ( double ) ( j * z[dim] ) / ( double ) ( m ), 1.0 );
    }
    quad2 = quad2 + f ( dim_num, x );
  }

  quad = quad1 + quad2 / ( double ) ( m );

  delete [] x;
  delete [] z;

  return quad;
}
//****************************************************************************80

double fibonacci_lattice_q ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_Q applies a Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    Because this is a standard lattice rule, it is really only suited
//    for functions which are periodic, of period 1, in both X and Y.
//
//    The related routines FIBONACCI_LATTICE_S and FIBONACCI_LATTICE_B
//    may be helpful in cases where the integrand does not satisfy this
//    requirement.
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_Q, the estimated integral.
//
{
  int dim;
  int dim_num = 2;
  int j;
  int m;
  double quad;
  double *x;
  int *z;

  x = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );

  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = fmod ( ( double ) ( j * z[dim] ) / ( double ) ( m ), 1.0 );
    }
    quad = quad + f ( dim_num, x );
  }
  quad = quad / ( double ) ( m );

  delete [] z;
  delete [] x;

  return quad;
}
//****************************************************************************80

double *fibonacci_lattice_q_nodes ( int k )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_Q_NODES returns Fibonacci lattice nodes in 2D.
//
//  Discussion:
//
//    Because this is a standard lattice rule, it is really only suited
//    for functions which are periodic, of period 1, in both X and Y.
//
//    The number of nodes returned is 
//
//      M = fibonacci ( k ).
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Output, double X[2*M], the nodes.
//
{
  int dim;
  int dim_num = 2;
  int j;
  int m;
  double *x;
  int *z;

  m = fibonacci ( k );

  x = new double[2*m];
  z = new int[dim_num];

  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim+j*dim_num] = fmod ( ( double ) ( j * z[dim] ) 
        / ( double ) ( m ), 1.0 );
    }
  }

  delete [] z;

  return x;
}
//****************************************************************************80

double fibonacci_lattice_q1 ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_Q1 applies a Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
//    It uses a nonlinear transformation on the integrand, which makes
//    the lattice rule more suitable for nonperiodic integrands.
//
//    The transformation replaces the integration variable X by
//
//      PHI(X) = 3*X^2 - 2*X**3
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_Q1, the estimated integral.
//
{
  int dim_num = 2;
  double dphi;
  int i;
  int j;
  int m;
  double quad;
  double *x;
  double *y;
  int *z;

  x = new double[dim_num];
  y = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );

  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = fmod ( ( double ) ( j * z[i] ) / ( double ) ( m ), 1.0 );
    }
    dphi = 1.0;
    for ( i = 0; i < dim_num; i++ )
    {
      y[i] = ( 3.0 - 2.0 * x[i] ) * x[i] * x[i];
      dphi = dphi * 6.0 * ( 1.0 - x[i] ) * x[i];
    }
    quad = quad + f ( dim_num, y ) * dphi;
  }

  quad = quad / ( double ) ( m );

  delete [] x;
  delete [] y;
  delete [] z;

  return quad;
}
//****************************************************************************80

double fibonacci_lattice_q2 ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_Q2 applies a Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
//    It uses a nonlinear transformation on the integrand, which makes
//    the lattice rule more suitable for nonperiodic integrands.
//
//    The transformation replaces the integration variable X by
//
//      PHI(X) = 3*X^3 *( 10 - 15 * X + 6 * X^2 )
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_Q2, the estimated integral.
//
{
  int dim_num = 2;
  double dphi;
  int i;
  int j;
  int m;
  double quad;
  double *x;
  double *y;
  int *z;

  x = new double[dim_num];
  y = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );

  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = fmod ( ( double ) ( j * z[i] ) / ( double ) ( m ), 1.0 );
    }
    dphi = 1.0;
    for ( i = 0; i < dim_num; i++ )
    {
      y[i] = ( 10.0 - 15.0 * x[i] + 6.0 * pow ( x[i], 2 ) ) * pow ( x[i], 3 );
      dphi = dphi * 30.0 * pow ( 1.0 - x[i], 2 ) * pow ( x[i], 2 ); 
    }
    quad = quad + f ( dim_num, y ) * dphi;
  }

  quad = quad / ( double ) ( m );

  delete [] x;
  delete [] y;
  delete [] z;

  return quad;
}
//****************************************************************************80

double fibonacci_lattice_q3 ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_Q3 applies a Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    This is a modification of the algorithm in FIBONACCI_LATTICE_Q.
//    It uses a nonlinear transformation on the integrand, which makes
//    the lattice rule more suitable for nonperiodic integrands.
//
//    The transformation replaces the integration variable X by
//
//      PHI(X) = X - sin ( 2 * PI * X ) / ( 2 * PI )
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_Q3, the estimated integral.
//
{
  int dim_num = 2;
  double dphi;
  int i;
  int j;
  int m;
  double quad;
  double pi = 3.141592653589793;
  double two_pi;
  double *x;
  double *y;
  int *z;

  x = new double[dim_num];
  y = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );

  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  two_pi = 2.0 * pi;

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = fmod ( ( double ) ( j * z[i] ) / ( double ) ( m ), 1.0 );
    }
    dphi = 1.0;
    for ( i = 0; i < dim_num; i++ )
    {
      y[i] = x[i] - sin ( two_pi * x[i] ) / two_pi;
      dphi = dphi * ( 1.0 - cos ( two_pi * x[i] ) );
    }
    quad = quad + f ( dim_num, y ) * dphi;
  }

  quad = quad / ( double ) ( m );

  delete [] x;
  delete [] y;
  delete [] z;

  return quad;
}
//****************************************************************************80

double fibonacci_lattice_t ( int k, double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    FIBONACCI_LATTICE_T applies a symmetric Fibonacci lattice integration rule in 2D.
//
//  Discussion:
//
//    This routine may be applied to integrands which are not periodic.
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int K, the index of the Fibonacci number to be used.
//    K must be at least 3.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double FIBONACCI_LATTICE_T, the estimated integral.
//
{
  int dim_num = 2;
  int i;
  int j;
  int m;
  double quad;
  double quad1;
  double quad2;
  int rank;
  double w;
  double *x;
  int *z;

  x = new double[dim_num];
  z = new int[dim_num];

  quad = 0.0;

  m = fibonacci ( k );
//
//  Get all the corner values.
//
  rank = 0;
  quad1 = 0.0;
  w = 1.0 / ( double ) i4_power ( 2, dim_num );

  for ( ; ; )
  {
    tuple_next ( 0, 1, dim_num, &rank, z );

    if ( rank == 0 )
    {
      break;
    }

    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = ( double ) z[i];
    }
    quad1 = quad1 + w * f ( dim_num, x );
  }
//
//  Get the interior values.
//
  z[0] = 1;
  z[1] = fibonacci ( k - 1 );

  quad2 = 0.0;
  for ( j = 1; j <= m - 1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = fmod ( ( double ) ( j * z[i] ) / ( double ) ( m ), 1.0 );
    }
    quad2 = quad2 + f ( dim_num, x );
  }

  quad = ( quad1 + quad2 ) / ( double ) ( m );

  delete [] x;
  delete [] z;

  return quad;
}
//****************************************************************************80

int *find_z20 ( int dim_num, int m )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_Z20 finds the appropriate Z vector to minimize P2(QS).
//
//  Discussion:
//
//    For the method of good lattice points, a number of points M, and
//    a single generator vector Z is chosen.  The integrand is assumed
//    to be periodic of period 1 in each argument, and is evaluated at
//    each of the points X(I)(1:DIM_NUM) = I * Z(1:DIM_NUM) / M, 
//    for I = 0 to M-1.  The integral is then approximated by the average
//    of these values.
//
//    Assuming that DIM_NUM and M are known, and that the integrand is not
//    known beforehand, the accuracy of the method depends entirely
//    on the choice of Z.  One method of choosing Z is to search for
//    the Z among all candidates which minimizes a particular quantity
//    P_ALPHA(QS).
//
//    We only need to look at vectors Z of the form
//    (1, L, L^2, ..., L^(DIM_NUM-1)),
//    for L = 1 to M/2.
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the number of points to be used.
//
//    Output, int FIND_Z20[DIM_NUM], the optimal vector.
//
{
  int dim;
  int i;
  double q0;
  double q0_min;
  int value;
  int *z;
  int *z_min;

  z = new int[dim_num];
  z_min = new int[dim_num];

  q0_min = r8_huge ( );

  for ( i = 1; i <= m / 2; i++ )
  {
    value = 1;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      z[dim] = value;
      value = ( value * i ) % m;
    }
//
//  Use this Z and the lattice integral method Q0 of order M,
//  to approximate the integral of P2.
//
    q0 = lattice ( dim_num, m, z, f20_s );
//
//  If this result is the smallest so far, save the corresponding Z.
//
    if ( q0 < q0_min )
    {
      q0_min = q0;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        z_min[dim] = z[dim];
      }
    }
  }
  delete [] z;
//
//  Return the best Z.
//
  return z_min;
}
//****************************************************************************80

void gray_next ( int n, int *change )

//****************************************************************************80
//
//  Purpose:
//
//    GRAY_NEXT generates the next Gray code by switching one item at a time.
//
//  Discussion:
//
//    On the first call only, the user must set CHANGE = -N.
//    This initializes the routine to the Gray code for N zeroes.
//
//    Each time it is called thereafter, it returns in CHANGE the index
//    of the item to be switched in the Gray code.  The sign of CHANGE
//    indicates whether the item is to be added or subtracted (or
//    whether the corresponding bit should become 1 or 0).  When
//    CHANGE is equal to N+1 on output, all the Gray codes have been
//    generated.
//
//    The routine has internal memory that is set up on call with
//    CHANGE = -N, and released on final return.
//
//  Example:
//
//    N  CHANGE         Subset in/out   Binary Number
//                      Interpretation  Interpretation
//                       1 2 4 8
//   --  ---------      --------------  --------------
//
//    4   -4 / 0         0 0 0 0         0
//
//        +1             1 0 0 0         1
//          +2           1 1 0 0         3
//        -1             0 1 0 0         2
//            +3         0 1 1 0         6
//        +1             1 1 1 0         7
//          -2           1 0 1 0         5
//        -1             0 0 1 0         4
//              +4       0 0 1 1        12
//        +1             1 0 1 1        13
//          +2           1 1 1 1        15
//        -1             0 1 1 1        14
//            -3         0 1 0 1        10
//        +1             1 1 0 1        11
//          -2           1 0 0 1         9
//        -1             0 0 0 1         8
//              -4       0 0 0 0         0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Albert Nijenhuis, Herbert Wilf,
//    Combinatorial Algorithms for Computers and Calculators,
//    Second Edition,
//    Academic Press, 1978,
//    ISBN: 0-12-519260-6,
//    LC: QA164.N54.
//
//  Parameters:
//
//    Input, int N, the order of the total set from which
//    subsets will be drawn.
//
//    Input/output, int *CHANGE.  This item is used for input only
//    on the first call for a particular sequence of Gray codes,
//    at which time it must be set to -N.  This corresponds to
//    all items being excluded, or all bits being 0, in the Gray code.
//    On output, CHANGE indicates which of the N items must be "changed", 
//    and the sign indicates whether the item is to be added or removed
//    (or the bit is to become 1 or 0).  Note that on return from the 
//    first call, CHANGE is set to 0, indicating that we begin with
//    the empty set.
//
{
  static int *a = NULL;
  int i;
  static int k = 0;
  static int n_save = -1;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "GRAY_NEXT - Fatal error!\n";
    cout << "  Input value of N <= 0.\n";
    exit ( 1 );
  }

  if ( *change == -n )
  {
    if ( a )
    {
      delete [] a;
    }

    a = new int[n];
    for ( i = 0; i < n; i++ )
    {
      a[i] = 0;
    }

    n_save = n;
    k = 1;
    *change = 0;

    return;
  }

  if ( n != n_save )
  {
    cout << "\n";
    cout << "GRAY_NEXT - Fatal error!\n";
    cout << "  Input value of N has changed from definition value.\n";
    exit ( 1 );
  }
//
//  First determine WHICH item is to be changed.
//
  if ( ( k % 2 ) == 1 )
  {
    *change = 1;
  }
  else
  {
    for ( i = 1; i <= n_save; i++ )
    {
      if ( a[i-1] == 1 )
      {
        *change = i + 1;
        break;
      }
    }
  }
//
//  Take care of the terminal case CHANGE = N + 1.
//
  if ( *change == n + 1 )
  {
    *change = n;
  }
//
//  Now determine HOW the item is to be changed.
//
  if ( a[*change-1] == 0 )
  {
    a[*change-1] = 1;
  }
  else if ( a[*change-1] == 1 )
  {
    a[*change-1] = 0;
    *change = -( *change );
  }
//
//  Update the counter.
//
  k = k + 1;
//
//  If the output CHANGE = -N_SAVE, then we're done.
//
  if ( *change == -n_save )
  {
    delete [] a;
    a = NULL;
    n_save = 0;
    k = 0;
  }

  return;
}
//****************************************************************************80

int i4_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_HUGE returns a "huge" I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int I4_HUGE, a "huge" I4.
//
{
  return 2147483647;
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

int i4_power ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_POWER returns the value of I^J.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    01 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, J, the base and the power.  J should be nonnegative.
//
//    Output, int I4_POWER, the value of I^J.
//
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J negative.\n";
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      cout << "\n";
      cout << "I4_POWER - Fatal error!\n";
      cout << "  I^J requested, with I = 0 and J = 0.\n";
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, int A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i < n; i++ ) 
  {
    cout << "  " << setw(8) << i 
         << ": " << setw(8) << a[i]  << "\n";
  }
  return;
}
//****************************************************************************80

double lattice ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    LATTICE applies a lattice integration rule.
//
//  Discussion:
//
//    Because this is a standard lattice rule, it is really only suited
//    for functions which are periodic, of period 1, in both X and Y.
//
//    For a suitable F, and a given value of M (the number of lattice points),
//    the performance of the routine is affected by the choice of the
//    generator vector Z.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 November 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the order (number of points) of the rule.
//
//    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
//    of Z satisfy 1 <= Z[*] < M, and are relatively prime to M.
//    This is easy to guarantee if M is itself a prime number.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double LATTICE, the estimated integral.
//
{
  int i;
  int j;
  double quad;
  double *x;

  x = new double[dim_num];

  quad = 0.0;

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i] = fmod ( ( double ) ( j * z[i] ) / ( double ) ( m ), 1.0 );
    }
    quad = quad + f ( dim_num, x );
  }

  quad = quad / ( double ) ( m );

  delete [] x;

  return quad;
}
//****************************************************************************80

double lattice_np0 ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    LATTICE_NP0 applies a lattice integration rule to a nonperiodic function.
//
//  Discussion:
//
//    This routine attempts to modify a lattice rule, suitable for use
//    with a periodic function, for use with a nonperiodic function F(X),
//    essentially by applying the lattice rule to the function
//
//      G(X) = ( F(X) + F(1-X) ) / 2
//
//    This is the rule in 1 dimension.  In two dimensions, we have
//
//      G(X,Y) = ( F(X,Y) + F(X,1-Y) + F(1-X,Y) + F(1-X,1-Y) ) / 4
//
//    with the obvious generalizations to higher dimensions.  
//
//    Drawbacks of this approach include:
//
//    * in dimension DIM_NUM, we must evaluate the function F at 
//      2**DIM_NUM points for every single evaluation of G;
//
//    * the function G, regarded as a periodic function, is continuous,
//      but not generally differentiable, at 0 and 1.
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
//    Seymour Haber,
//    Parameters for Integrating Periodic Functions of Several Variables,
//    Mathematics of Computation,
//    Volume 41, Number 163, July 1983, pages 115-129.
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the order (number of points) of the rule.
//
//    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
//    of Z satisfy 1 <= Z[*] < M, and are relatively prime to M.
//    This is easy to guarantee if M is itself a prime number.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double LATTICE_NP0, the estimated integral.
//
{
  int dim;
  int gray_bit;
  int j;
  double quad;
  double *x;
  double *y;

  x = new double[dim_num];
  y = new double[dim_num];

  quad = 0.0;

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = fmod ( ( double ) ( j * z[dim] ) / ( double ) ( m ), 1.0 );
    }
//
//  Generate all DIM_NUM-tuples for which the I-th element is X(I) or 1-X(I).
//
    gray_bit = - dim_num;

    for ( ; ; )
    {
      gray_next ( dim_num, &gray_bit );

      if ( gray_bit == - dim_num )
      {
        break;
      }

      if ( gray_bit == 0 )
      {
        for ( dim = 0; dim < dim_num; dim++ )
        {
          y[dim] = x[dim];
        }
      }
      else
      {
        dim = abs ( gray_bit ) - 1;
        y[dim] = 1.0 - y[dim];
      }
      quad = quad + f ( dim_num, y );
    }
  }

  quad = quad / ( double ) ( i4_power ( 2, dim_num ) * m );

  delete [] x;
  delete [] y;

  return quad;
}
//****************************************************************************80

double lattice_np1 ( int dim_num, int m, int z[], 
  double f ( int dim_num, double x[] ) )

//****************************************************************************80
//
//  Purpose:
//
//    LATTICE_NP1 applies a lattice integration rule to a nonperiodic function.
//
//  Discussion:
//
//    This routine applies the transformation function
//
//      PHI(T) = 3*T^2 - 2*T^3
//
//    to a nonperiodic integrand to make it suitable for treatment
//    by a lattice rule.
//
//    For a suitable F, and a given value of M (the number of lattice points),
//    the performance of the routine is affected by the choice of the
//    generator vector Z.
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the order (number of points) of the rule.
//
//    Input, int Z[DIM_NUM], the generator vector.  Typically, the elements
//    of Z satisfy 1 <= Z(1:DIM_NUM) < M, and are relatively prime to M.
//    This is easy to guarantee if M is itself a prime number.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Output, double LATTICE_NP1, the estimated integral.
//
{
  int dim;
  double dphi;
  int i;
  int j;
  double quad;
  double *x;
  double *y;

  x = new double[dim_num];
  y = new double[dim_num];

  quad = 0.0;

  for ( j = 0; j <= m - 1; j++ )
  {
    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = fmod ( ( double ) ( j * z[dim] ) / ( double ) ( m ), 1.0 );
    }
    dphi = 1.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      y[dim] = ( 3.0 - 2.0 * x[dim] ) * x[dim] * x[dim];
      dphi = dphi * 6.0 * ( 1.0 - x[dim] ) * x[dim];
    }
    quad = quad + f ( dim_num, y ) * dphi;
  }

  quad = quad / ( double ) ( m );

  delete [] x;
  delete [] y;

  return quad;
}
//****************************************************************************80

void lattice_print ( int dim_num, int m, int z[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    LATTICE_PRINT prints the points in a lattice rule.
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
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the number of points to use.
//
//    Input, int Z[DIM_NUM], the generator vector.
//
//    Input, string TITLE, a title.
//
{
  int dim;
  int i;
  int *y;

  y = new int[dim_num];

  cout << "\n";
  cout << title << "\n";
  cout << "\n";

  for ( i = 0; i <= m - 1; i++ )
  {
    cout << setw(4) << i + 1 << "    ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      y[dim] = ( i * z[dim] ) % m;
      cout << setw(4) << y[dim];
    }
    cout << "\n";
  }

  delete [] y;

  return;
}
//****************************************************************************80

double monte_carlo ( int dim_num, int m, double f ( int dim_num, double x[] ),
  int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MONTE_CARLO applies a Monte Carlo integration rule.
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
//  Reference:
//
//    Ian Sloan, Stephen Joe,
//    Lattice Methods for Multiple Integration,
//    Oxford, 1994,
//    ISBN: 0198534728,
//    LC: QA311.S56
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int M, the number of points to use.
//
//    Input, double F ( int DIM_NUM, double X[] ), the name of the 
//    user-supplied routine which evaluates the function.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double MONTE_CARLO, the estimated integral.
//
{
  int j;
  double quad;
  double *x;

  quad = 0.0;

  for ( j = 1; j <= m; j++ )
  {
    x = r8vec_uniform_01 ( dim_num, seed );
    quad = quad + f ( dim_num, x );
    delete [] x;
  }

  quad = quad / ( double ) m;

  return quad;
}
//****************************************************************************80

int prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, 
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, 
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, 
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, 
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, 
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, 
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, 
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, 
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, 
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, 
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, 
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, 
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, 
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, 
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, 
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, 
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, 
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, 
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, 
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, 
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, 
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, 
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, 
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, 
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, 
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, 
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, 
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, 
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, 
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, 
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, 
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, 
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, 
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, 
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, 
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, 
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, 
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, 
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, 
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, 
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, 
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, 
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, 
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, 
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, 
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, 
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641, 
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739, 
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829, 
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923, 
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007, 
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109, 
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187, 
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309, 
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411, 
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cout << "\n";
    cout << "PRIME - Fatal error//\n";
    cout << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

double r8_abs ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ABS returns the absolute value of an R8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the quantity whose absolute value is desired.
//
//    Output, double R8_ABS, the absolute value of X.
//
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  } 
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, string TITLE, a title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo, 
  int ihi, int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input, double A[M*N], an M by N matrix to be printed.
//
//    Input, int ILO, JLO, the first row and column to print.
//
//    Input, int IHI, JHI, the last row and column to print.
//
//    Input, string TITLE, a title.
//
{
# define INCX 5

  int i;
  int i2;
  int i2hi;
  int i2lo;
  int inc;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  for ( i2lo = i4_max ( ilo, 1 ); i2lo <= i4_min ( ihi, m ); i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    inc = i2hi + 1 - i2lo;

    cout << "\n";
    cout << "  Row: ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i - 1 << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j - 1 << ":";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}
//****************************************************************************80

double *r8vec_uniform_01 ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
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
//    19 August 2004
//
//  Author:
//
//    John Burkardt
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
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
  int i;
  int k;
  double *r;

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + 2147483647;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

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
//****************************************************************************80

void tuple_next ( int m1, int m2, int n, int *rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT computes the next element of a tuple space.
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between M1 and M2.  The elements are produced one at a time.
//    The first element is
//      (M1,M1,...,M1),
//    the second element is
//      (M1,M1,...,M1+1),
//    and the last element is
//      (M2,M2,...,M2)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2, M1 = 1, M2 = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      Rank   X
//    ----  ---    -----  ---
//    0     * *    1      1 1
//    1     1 1    2      1 2
//    2     1 2    3      1 3
//    3     1 3    4      2 1
//    4     2 1    5      2 2
//    5     2 2    6      2 3
//    6     2 3    7      3 1
//    7     3 1    8      3 2
//    8     3 2    9      3 3
//    9     3 3    0      0 0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M1, M2, the minimum and maximum entries.
//
//    Input, int N, the number of components.
//
//    Input/output, int *RANK, counts the elements.
//    On first call, set RANK to 0.  Thereafter, the output value of RANK
//    will indicate the order of the element returned.  When there are no
//    more elements, RANK will be returned as 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  int i;
  int j;

  if ( m2 < m1 )
  {
    *rank = 0;
    return;
  }

  if ( *rank <= 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = m1;
    }
    *rank = 1;
  }
  else
  {
    *rank = *rank + 1;
    i = n - 1;

    for ( ; ; )
    {

      if ( x[i] < m2 )
      {
        x[i] = x[i] + 1;
        break;
      }

      x[i] = m1;

      if ( i == 0 )
      {
        *rank = 0;
        for ( j = 0; j < n; j++ )
        {
          x[j] = m1;
        }
        break;
      }
      i = i - 1;
    }
  }
  return;
}
