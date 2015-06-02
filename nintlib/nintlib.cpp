# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "nintlib.hpp"

//****************************************************************************80

double box_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  int order, double xtab[], double weight[], int *eval_num )

//****************************************************************************80
//
//  Purpose:
//
//    BOX_ND estimates a multidimensional integral using a product rule.
//
//  Discussion:
//
//    The routine creates a DIM_NUM-dimensional product rule from a 1D rule
//    supplied by the user.  The routine is fairly inflexible.  If
//    you supply a rule for integration from -1 to 1, then your product
//    box must be a product of DIM_NUM copies of the interval [-1,1].
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
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//      
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int ORDER, the number of points used in the 1D rule.
//
//    Input, double XTAB[ORDER], the abscissas of the 1D rule.
//
//    Input, double WEIGHT[ORDER], the weights of the 1D rule.
//
//    Output, int *EVAL_NUM, the number of function evaluations.
//
//    Output, double BOX_ND, the approximate value of the integral.
//
{
  int dim;
  int *indx;
  int k;
  double result;
  double w;
  double *x;

  *eval_num = 0;

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "BOX_ND - Fatal error!\n";
    cout << "  DIM_NUM < 1.\n";
    cout << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  if ( order < 1 )
  {
    cout << "\n";
    cout << "BOX_ND - Fatal error!\n";
    cout << "  ORDER < 1.\n";
    cout << "  ORDER = " << order << "\n";
    exit ( 1 );
  }

  k = 0;
  result = 0.0;

  indx = new int[dim_num];
  x = new double[dim_num];

  for ( ; ; )
  {
    tuple_next ( 1, order, dim_num, &k, indx );

    if ( k == 0  )
    {
      break;
    }

    w = 1.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      w = w * weight[indx[dim]-1];
    }

    for ( dim = 0; dim < dim_num; dim++ )
    {
      x[dim] = xtab[indx[dim]-1];
    }

    result = result + w * func ( dim_num, x );
    *eval_num = *eval_num + 1;
  }

  delete [] indx;
  delete [] x;

  return result;
}
//****************************************************************************80

int i4_huge ( void )

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

double monte_carlo_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  double a[], double b[], int eval_num, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    MONTE_CARLO_ND estimates a multidimensional integral using Monte Carlo.
//
//  Discussion:
//
//    Unlike the other routines, this routine requires the user to specify
//    the number of function evaluations as an INPUT quantity.  
//    No attempt at error estimation is made.
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
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
//
//    Input, int EVAL_NUM, the number of function evaluations.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double MONTE_CARLO_ND, the approximate value of the integral.
//
{
  int dim;
  int i;
  double result;
  double volume;
  double *x;

  result = 0.0;

  for ( i = 0; i < eval_num; i++ )
  {
    x = r8vec_uniform_01_new ( dim_num, seed );

    result = result + func ( dim_num, x );

    delete [] x;
  }

  volume = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    volume = volume * ( b[dim] - a[dim] );
  }

  result = result * volume / ( double ) ( eval_num );

  return result;
}
//****************************************************************************80

double p5_nd ( double func ( int dim_num, double x[] ), int dim_num, 
  double a[], double b[], int *eval_num )

//****************************************************************************80
//
//  Purpose:
//
//    P5_ND estimates a multidimensional integral with a formula of exactness 5.
//
//  Discussion:
//
//    The routine uses a method which is exact for polynomials of total 
//    degree 5 or less.
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
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
//
//    Output, int *EVAL_NUM, the number of function evaluations.
//
//    Output, double P5_ND, the approximate value of the integral.
//
{
  double a0;
  double a1;
  double a2;
  double a3;
  double a4;
  double a5;
  int dim;
  double en;
  int i;
  int j;
  double result;
  double sum1;
  double sum2;
  double sum3;
  double volume;
  double *work;

  *eval_num = 0;

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "P5_ND - Fatal error!\n";
    cout << "  DIM_NUM < 1, DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }
 
  a2 = 25.0 / 324.0;
  a3 = sqrt ( 0.6 );
  en = ( double ) ( dim_num );
  a0 = ( 25.0 * en * en - 115.0 * en + 162.0 ) / 162.0;
  a1 = ( 70.0 - 25.0 * en ) / 162.0;
 
  volume = 1.0;
  for ( dim = 0; dim < dim_num; dim++ )
  {
    volume = volume * ( b[dim] - a[dim] );
  }

  work = new double[dim_num];
  for ( dim = 0; dim < dim_num; dim++ )
  {
    work[dim] = 0.5 * ( a[dim] + b[dim] );
  }
 
  result = 0.0;
  if ( volume == 0.0 )
  {
    cout << "\n";
    cout << "P5_ND - Warning!\n";
    cout << "  Volume = 0, integral = 0.\n";
    delete [] work;
    return result;
  }
 
  sum1 = a0 * func ( dim_num, work );
  *eval_num = *eval_num + 1;

  sum2 = 0.0;
  sum3 = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    work[i] = 0.5 * ( ( a[i] + b[i] ) + a3 * ( b[i] - a[i] ) );
    sum2 = sum2 + func ( dim_num, work );
    *eval_num = *eval_num + 1;

    work[i] = 0.5 * ( ( a[i] + b[i] ) - a3 * ( b[i] - a[i] ) );
    sum2 = sum2 + func ( dim_num, work );
    *eval_num = *eval_num + 1;

    work[i] = 0.5 * ( a[i] + b[i] );
  }

  if ( 1 < dim_num )
  {
    a4 = a3;
 
    for ( ; ; )
    {
      for ( i = 0; i < dim_num - 1; i++ )
      {
        work[i] = 0.5 * ( ( a[i] + b[i] ) + a4 * ( b[i] - a[i] ) );
        a5 = a3;
 
        for ( ; ; )
        { 
          for ( j = i + 1; j < dim_num; j++ )
          {
            work[j] = 0.5 * ( ( a[j] + b[j] ) + a5 * ( b[j] - a[j] ) );
            sum3 = sum3 + func ( dim_num, work );
            *eval_num = *eval_num + 1;
            work[j] = 0.5 * ( a[j]+ b[j] );
          }
 
          a5 = -a5;

          if ( 0.0 <= a5 )
          {
            break;
          }
        }
        work[i] = 0.5 * ( a[i] + b[i] );
      }
 
      a4 = -a4;

      if ( 0.0 <= a4 )
      {
        break;
      }
    }
  }
 
  result = volume * ( sum1 + a1 * sum2 + a2 * sum3 );

  delete [] work;
 
  return result;
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
    value = -x;
  }
  return value;
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
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
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int k;
  double *r;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = *seed / 127773;

    *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

    if ( *seed < 0 )
    {
      *seed = *seed + i4_huge ( );
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double romberg_nd ( double func ( int dim_num, double x[] ), double a[], 
  double b[], int dim_num, int sub_num[], int it_max, double tol, int *ind, 
  int *eval_num )

//****************************************************************************80
//
//  Purpose:
//
//    ROMBERG_ND estimates a multidimensional integral using Romberg integration.
//
//  Discussion:
//
//    The routine uses a Romberg method based on the midpoint rule.
//
//    In the reference, this routine is called "NDIMRI".
//
//    Thanks to Barak Bringoltz for pointing out problems in a previous
//    FORTRAN90 implementation of this routine.
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
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
//    Input, double A[DIM_NUM], B[DIM_NUM], the integration limits.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int SUB_NUM[DIM_NUM], the number of subintervals into
//    which the I-th integration interval (A(I), B(I)) is
//    initially subdivided.  SUB_NUM(I) must be greater than 0.
//
//    Input, int IT_MAX, the maximum number of iterations to
//    be performed.  The number of function evaluations on
//    iteration J is at least J**DIM_NUM, which grows very rapidly.
//    IT_MAX should be small!
//
//    Input, double TOL, an error tolerance for the approximation
//    of the integral.
//
//    Output, int *IND, error return flag.
//    IND = -1 if the error tolerance could not be achieved.
//    IND = 1 if the error tolerance was achieved.
//
//    Output, int *EVAL_NUM, the number of function evaluations.
//
//    Output, double ROMBERG_ND, the approximate value of the integral.
//
//  Local Parameters:
//
//    Local, int IWORK[DIM_NUM], a pointer used to generate all the
//    points X in the product region.
//
//    Local, int IWORK2[IT_MAX], a counter of the number of points
//    used at each step of the Romberg iteration.
//
//    Local, int SUB_NUM2[DIM_NUM], the number of subintervals used
//    in each direction, a refinement of the user's input SUB_NUM.
//
//    Local, double TABLE[IT_MAX], the difference table.
//
//    Local, double X[DIM_NUM], an evaluation point.
//
{
  int dim;
  double en;
  double factor;
  int i;
  int it;
  int *iwork;
  int *iwork2;
  int kdim;
  int ll;
  double result;
  double result_old;
  double rnderr;
  double submid;
  int *sub_num2;
  double sum1;
  double weight;
  double *table;
  double *x;

  *eval_num = 0;

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "ROMBERG_ND - Fatal error!\n";
    cout << "  DIM_NUM is less than 1.  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }
 
  if ( it_max < 1 )
  {
    cout << "\n";
    cout << "ROMBERG_ND - Fatal error!\n";
    cout << "  IT_MAX is less than 1.  IT_MAX = " << it_max << "\n";
    exit ( 1 );
  }
 
  for ( i = 0; i < dim_num; i++ )
  {
    if ( sub_num[i] <= 0 )
    {
      cout << "\n";
      cout << "ROMBERG_ND - Fatal error!\n";
      cout << "  SUB_NUM(I) is less than 1.\n";
      cout << "  for I = " << i << "\n";
      cout << "  SUB_NUM(I) = " << sub_num[i] << "\n";
      exit ( 1 );
    }
  }
 
  iwork = new int[dim_num];
  iwork2 = new int[it_max];
  sub_num2 = new int[dim_num];
  table = new double[it_max];
  x = new double[dim_num];

  *ind = 0;
  rnderr = r8_epsilon ( );
  iwork2[0] = 1;

  for ( dim = 0; dim < dim_num; dim++ )
  {
    sub_num2[dim] = sub_num[dim];
  } 
  if ( 1 < it_max )
  {
    iwork2[1] = 2;
  }
 
  it = 1;
 
  for ( ; ; )
  {
    sum1 = 0.0;
 
    weight = 1.0;
    for ( dim = 0; dim < dim_num; dim++ )
    {
      weight = weight * ( b[dim] - a[dim] ) / ( double ) sub_num2[dim];
    }
//
//  Generate every point X in the product region, and evaluate F(X).
//
    for ( dim = 0; dim < dim_num; dim++ )
    {
      iwork[dim] = 1;
    }
    
    for ( ; ; )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        x[dim] = 
         ( ( double ) ( 2 * sub_num2[dim] - 2 * iwork[dim] + 1 ) * a[dim] 
         + ( double ) (                   + 2 * iwork[dim] - 1 ) * b[dim] )
         / ( double ) ( 2 * sub_num2[dim]                      );
      }

      sum1 = sum1 + func ( dim_num, x );
      *eval_num = *eval_num + 1;

      kdim = dim_num;

      while ( 0 < kdim ) 
      {
        if ( iwork[kdim-1] < sub_num2[kdim-1] )
        {
          iwork[kdim-1] = iwork[kdim-1] + 1;
          break;
        }
        iwork[kdim-1] = 1;
        kdim = kdim - 1;
      }

      if ( kdim == 0 )
      {
        break;
      }
    }
//
//  Done with summing.
//
    table[it-1] = weight * sum1;

    if ( it <= 1 )
    {
      result = table[0];
      result_old = result;

      if ( it_max <= it )
      {
        *ind = 1;
        break;
      }

      it = it + 1;
      for ( dim = 0; dim < dim_num; dim++ )
      {
        sub_num2[dim] = iwork2[it-1] * sub_num2[dim];
      }
      continue;
    }
//
//  Compute the difference table for Richardson extrapolation.
// 
    for ( ll = 2; ll <= it; ll++ )
    {
      i = it + 1 - ll;
      factor = ( double ) ( iwork2[i-1] * iwork2[i-1] ) 
        / ( double ) ( iwork2[it-1] * iwork2[it-1] - iwork2[i-1] * iwork2[i-1] );
      table[i] = table[i] + ( table[i] - table[i-1] ) * factor;
    }
 
    result = table[0];
//
//  Terminate successfully if the estimated error is acceptable.
//
    if ( r8_abs ( result - result_old ) <= r8_abs ( result * ( tol + rnderr ) ) )
    {
      *ind = 1;
      break;
    }
//
//  Terminate unsuccessfully if the iteration limit has been reached.
//
    if ( it_max <= it )
    {
      *ind = -1;
      break;
    }
//
//  Prepare for another step.
//
    result_old = result;

    it = it + 1;

    iwork2[it-1] = ( int ) ( 1.5 * ( double ) ( iwork2[it-2] ) );

    for ( dim = 0; dim < dim_num; dim++ )
    {
      sub_num2[dim] = ( int ) ( 1.5 * ( double ) ( sub_num2[dim] ) );
    }
  }

  delete [] iwork;
  delete [] iwork2;
  delete [] sub_num2;
  delete [] table;
  delete [] x;

  return result;
}
//****************************************************************************80

void sample_nd ( double func ( int dim_num, double x[] ), int k1, int k2, 
  int dim_num, double est1[], double err1[], double dev1[], double est2[], 
  double err2[], double dev2[], int *eval_num )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_ND estimates a multidimensional integral using sampling.
//
//  Discussion:
//
//    This routine computes two sequences of integral estimates, EST1 
//    and EST2, for indices K going from K1 to K2.  These estimates are 
//    produced by the generation of 'random' abscissas in the region.  
//    The process can become very expensive if high accuracy is needed.
//
//    The total number of function evaluations is
//    4*(K1^DIM_NUM+(K1+1)^DIM_NUM+...+(K2-1)^DIM_NUM+K2^DIM_NUM), and K2
//    should be chosen so as to make this quantity reasonable.
//    In most situations, EST2(K) are much better estimates than EST1(K).
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
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
//    Input, int K1, the beginning index for the iteration.
//    1 <= K1 <= K2.
//
//    Input, int K2, the final index for the iteration.  K1 <= K2.
//    Increasing K2 increases the accuracy of the calculation,
//    but vastly increases the work and running time of the code.
//
//    Input, int DIM_NUM, the spatial dimension.  1 <= DIM_NUM <= 10.
//
//    Output, double EST1[K2].  Entries K1 through K2 contain
//    successively better estimates of the integral.
//
//    Output, double ERR1[K2].  Entries K1 through K2 contain
//    the corresponding estimates of the integration errors.
//
//    Output, double DEV1[K2].  Entries K1 through K2 contain
//    estimates of the reliability of the the integration.
//    If consecutive values DEV1(K) and DEV1(K+1) do not differ
//    by more than 10 percent, then ERR1(K) can be taken as
//    a reliable upper bound on the difference between EST1(K)
//    and the true value of the integral.
//
//    Output, double EST2[K2].  Entries K2 through K2 contain
//    successively better estimates of the integral.
//
//    Output, double ERR2[K2].  Entries K2 through K2 contain
//    the corresponding estimates of the integration errors.
//
//    Output, double DEV2[K2].  Entries K2 through K2 contain
//    estimates of the reliability of the the integration.
//    If consecutive values DEV2(K) and DEV2(K+2) do not differ
//    by more than 10 percent, then ERR2(K) can be taken as
//    a reliable upper bound on the difference between EST2(K)
//    and the true value of the integral.
//
//    Output, int *EVAL_NUM, the number of function evaluations.
//
{
# define DIM_MAX 10

  double ak;
  double ak1;
  double akn;
  double al[DIM_MAX] = {
       0.4142135623730950, 
       0.7320508075688773, 
       0.2360679774997897, 
       0.6457513110645906, 
       0.3166247903553998, 
       0.6055512754639893, 
       0.1231056256176605, 
       0.3589989435406736, 
       0.7958315233127195, 
       0.3851648071345040 };
  double b;
  double *be;
  double bk;
  double d1;
  double d2;
  double *dex;
  int dim;
  double g;
  double *ga;
  int i;
  int j;
  int k;
  int key;
  bool more;
  double *p1;
  double *p2;
  double *p3;
  double *p4;
  double s1;
  double s2;
  double t;
  double y1;
  double y2;
  double y3;
  double y4;

  *eval_num = 0;
//
//  Check input.
//
  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "SAMPLE_ND - Fatal error!\n";
    cout << "  DIM_NUM must be at least 1,\n";
    cout << "  but DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }
 
  if ( DIM_MAX < dim_num )
  {
    cout << "\n";
    cout << "SAMPLE_ND - Fatal error!\n";
    cout << "  DIM_NUM must be no more than DIM_MAX = " << DIM_MAX << "\n";
    cout << "  but DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }
 
  if ( k1 < 1 )
  {
    cout << "\n";
    cout << "SAMPLE_ND - Fatal error!\n";
    cout << "  K1 must be at least 1, but K1 = " << k1 << "\n";
    exit ( 1 );
  }
 
  if ( k2 < k1 )
  {
    cout << "\n";
    cout << "SAMPLE_ND - Fatal error!\n";
    cout << "  K1 may not be greater than K2, but \n";
    cout << "  K1 = " << k1 << "\n";
    cout << "  K2 = " << k2 << "\n";
    exit ( 1 );
  }
 
  be = new double[dim_num];
  dex = new double[dim_num];
  ga = new double[dim_num];
  p1 = new double[dim_num];
  p2 = new double[dim_num];
  p3 = new double[dim_num];
  p4 = new double[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    be[dim] = al[dim];
  }
 
  for ( dim = 0; dim < dim_num; dim++ )
  {
    ga[dim] = al[dim];
  }
 
  for ( dim = 0; dim < dim_num; dim++ )
  {
    dex[dim] = 0.0;
  }
 
  for ( k = k1; k <= k2; k++ )
  {
    ak = ( double ) ( k );
    key = 0;
    ak1 = ak - 1.1;
    s1 = 0.0;
    d1 = 0.0;
    s2 = 0.0;
    d2 = 0.0;
    akn = pow ( ak, dim_num );
    t = sqrt ( pow ( ak, dim_num ) ) * ak;
    bk = 1.0 / ak;
 
    for ( ; ; )
    { 
      key = key + 1;

      if ( key != 1 )
      {
        key = key - 1;
        more = false;
        for ( j = 0; j < dim_num; j++ )
        {
          if ( dex[j] <= ak1 )
          {
            dex[j] = dex[j] + 1.0;
            more = true;
            break;
          }
 
          dex[j] = 0.0;
        }

        if ( !more )
        {
          break;
        }
      }
 
      for ( i = 0; i < dim_num; i++ )
      {
        b = be[i] + al[i];
        if ( 1.0 < b )
        {
          b = b - 1.0;
        }

        g = ga[i] + b;
        if ( 1.0 < g )
        {
          g = g - 1.0;
        }

        be[i] = b + al[i];
        if ( 1.0 < be[i] )
        {
          be[i] = be[i] - 1.0;
        }

        ga[i] = be[i] + g;
        if ( 1.0 < ga[i] )
        {
          ga[i] = ga[i] - 1.0;
        }

        p1[i] = ( dex[i] + g ) * bk;
        p2[i] = ( dex[i] + 1.0 - g ) * bk;
        p3[i] = ( dex[i] + ga[i] ) * bk;
        p4[i] = ( dex[i] + 1.0 - ga[i] ) * bk;
      }
 
      y1 = func ( dim_num, p1 );
      *eval_num = *eval_num + 1;
//
//  There may be an error in the next two lines,
//  but oddly enough, that is how the original reads
//
      y3 = func ( dim_num, p2 );
      *eval_num = *eval_num + 1;
      y2 = func ( dim_num, p3 );
      *eval_num = *eval_num + 1;
      y4 = func ( dim_num, p4 );
      *eval_num = *eval_num + 1;

      s1 = s1 + y1 + y2;
      d1 = d1 + ( y1 - y2 ) * ( y1 - y2 );
      s2 = s2 + y3 + y4;
      d2 = d2 + ( y1 + y3 - y2 - y4 ) * ( y1 + y3 - y2 - y4 );
    }
  
    est1[k-1] = 0.5 * s1 / akn;
    err1[k-1] = 1.5 * sqrt ( d1 ) / akn;
    dev1[k-1] = err1[k-1] * t;
    est2[k-1] = 0.25 * ( s1 + s2 ) / akn;
    err2[k-1] = 0.75 * sqrt ( d2 ) / akn;
    dev2[k-1] = err2[k-1] * t * ak; 
  }
 
  delete [] be;
  delete [] dex;
  delete [] ga;
  delete [] p1;
  delete [] p2;
  delete [] p3;
  delete [] p4;

  return;
# undef DIM_MAX
}
//****************************************************************************80

double sum2_nd ( double func ( int dim_num, double x[] ), double xtab[], 
  double weight[], int order[], int dim_num, int *eval_num )

//****************************************************************************80
//
//  Purpose:
//
//    SUM2_ND estimates a multidimensional integral using a product rule.
//
//  Discussion:
//
//    The routine uses a product rule supplied by the user.
//
//    The region may be a product of any combination of finite,
//    semi-infinite, or infinite intervals.
//
//    For each factor in the region, it is assumed that an integration
//    rule is given, and hence, the region is defined implicitly by
//    the integration rule chosen.
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
//    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
//    C++ version by John Burkardt
//
//  Reference:
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//  Parameters:
//
//    Input, double FUNC ( int dim_num, double x[] ), evaluates
//    the function to be integrated.
//
//    Input, double XTAB[DIM_NUM*ORDER_MAX].  XTAB(I,J) is the 
//    I-th abscissa of the J-th rule.
//
//    Input, double WEIGHT[DIM_NUM*ORDER_MAX].  WEIGHT(I,J) is the 
//    I-th weight for the J-th rule.
//
//    Input, int ORDER[DIM_NUM].  ORDER(I) is the number of
//    abscissas to be used in the J-th rule.  ORDER(I) must be
//    greater than 0 and less than or equal to ORDER_MAX.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Output, int EVAL_NUM, the number of function evaluations.
//
//    Output, double SUM2_ND, the approximate value of the integral.
//
{
  int dim;
  int i;
  int *iwork;
  int k;
  int m1;
  double result;
  double w1;
  double *work;
//
//  Default values.
//
  result = 0.0;
  *eval_num = 0;

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "SUM2_ND - Fatal error!\n";
    cout << "  DIM_NUM < 1\n";
    cout << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }
 
  for ( i = 0; i < dim_num; i++ )
  {
    if ( order[i] < 1 )
    {
      cout << "\n";
      cout << "SUM2_ND - Fatal error!\n";
      cout << "  ORDER(I) < 1.\n";
      cout << "  For I = " << i << "\n";
      cout << "  ORDER(I) = " << order[i] << "\n";
      exit ( 1 );
    }
  }

  iwork = new int[dim_num];
  work = new double[dim_num];

  for ( dim = 0; dim < dim_num; dim++ )
  {
    iwork[dim] = 1;
  }
 
  for ( ; ; )
  {
    k = 1;
 
    w1 = 1.0;
    for ( i = 0; i < dim_num; i++ )
    {
      m1 = iwork[i];
      work[i] = xtab[i+(m1-1)*dim_num];
      w1 = w1 * weight[i+(m1-1)*dim_num];
    }

    result = result + w1 * func ( dim_num, work );
    *eval_num = *eval_num + 1;

    while ( iwork[k-1] == order[k-1] )
    {
      iwork[k-1] = 1;
      k = k + 1;

      if ( dim_num < k )
      {
        return result;
      }
    }
    iwork[k-1] = iwork[k-1] + 1;
  }

  delete [] iwork;
  delete [] work;

  return result;
}
//****************************************************************************80

void timestamp ( void )

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
