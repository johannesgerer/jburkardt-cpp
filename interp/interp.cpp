# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "interp.hpp"

//****************************************************************************80

double *cc_abscissas ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_ABSCISSAS computes the Clenshaw Curtis abscissas.
//
//  Discussion:
//
//    The interval is [ -1, 1 ].
//
//    The abscissas are the cosines of equally spaced angles between
//    180 and 0 degrees, including the endpoints.
//
//      X(I) = cos ( ( ORDER - I ) * PI / ( ORDER - 1 ) )
//
//    except for the basic case ORDER = 1, when
//
//      X(1) = 0.
//
//    If the value of ORDER is increased in a sensible way, then
//    the new set of abscissas will include the old ones.  One such
//    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
//
//    When doing interpolation with Lagrange polynomials, the Clenshaw Curtis
//    abscissas can be a better choice than regularly spaced abscissas.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Charles Clenshaw, Alan Curtis,
//    A Method for Numerical Integration on an Automatic Computer,
//    Numerische Mathematik,
//    Volume 2, Number 1, December 1960, pages 197-205.
//
//    Philip Davis, Philip Rabinowitz,
//    Methods of Numerical Integration,
//    Second Edition,
//    Dover, 2007,
//    ISBN: 0486453391,
//    LC: QA299.3.D28.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double CC_ABSCISSAS[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CC_ABSCISSA - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( n - 1 - i )
                  / ( double ) ( n - 1 );
    x[i] = cos ( theta );
  }

  return x;
}
//****************************************************************************80

double *cc_abscissas_ab ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    CC_ABSCISSAS_AB computes Clenshaw Curtis abscissas for the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the order of the rule.
//
//    Output, double CC_ABSCISSAS_AB[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "CC_ABSCISSAS_AB - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.5 * ( b + a );
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    theta = r8_pi * ( double ) ( n - 1 - i )
                  / ( double ) ( n - 1 );
    x[i] = 0.5 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
  }

  return x;
}
//****************************************************************************80

double *f1_abscissas ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F1_ABSCISSAS computes Fejer type 1 abscissas.
//
//  Discussion:
//
//    The interval is [ -1, +1 ].
//
//    The abscissas are the cosines of equally spaced angles, which
//    are the midpoints of N equal intervals between 0 and PI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
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
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double F1_ABSCISSAS[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "F1_ABSCISSAS - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( 2 * n - 2 * i - 1 ) * r8_pi 
          / ( double ) ( 2 * n             );
    x[i] = cos ( theta );
  }

  return x;
}
//****************************************************************************80

double *f1_abscissas_ab ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F1_ABSCISSAS_AB computes Fejer type 1 abscissas for the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
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
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the order of the rule.
//
//    Output, double F1_ABSCISSAS_AB[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "F1_ABSCISSAS_AB - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.5 * ( b + a );
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( 2 * n - 2 * i - 1 ) * r8_pi 
          / ( double ) ( 2 * n             );
    x[i] = 0.5 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
  }

  return x;
}
//****************************************************************************80

double *f2_abscissas ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    F2_ABSCISSAS computes Fejer Type 2 abscissas.
//
//  Discussion:
//
//    The interval is [-1,+1].
//
//    The abscissas are the cosines of equally spaced angles.
//    The angles are computed as N+2 equally spaced values between 0 and PI,
//    but with the first and last angle omitted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
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
//    Walter Gautschi,
//    Numerical Quadrature in the Presence of a Singularity,
//    SIAM Journal on Numerical Analysis,
//    Volume 4, Number 3, 1967, pages 357-362.
//
//    Joerg Waldvogel,
//    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
//    BIT Numerical Mathematics,
//    Volume 43, Number 1, 2003, pages 1-18.
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double F2_ABSCISSAS[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "F2_ABSCISSAS - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
    return x;
  }
  else if ( n == 2 )
  {
    x[0] = -0.5;
    x[1] =  0.5;
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( n - i ) * r8_pi 
          / ( double ) ( n + 1 );
    x[i] = cos ( theta );
  }

  return x;
}
//****************************************************************************80

double *f2_abscissas_ab ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    F2_ABSCISSAS_AB computes Fejer Type 2 abscissas for the interval [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the order of the rule.
//
//    Output, double F2_ABSCISSAS_AB[N], the abscissas.
//
{
  int i;
  const double r8_pi = 3.141592653589793;
  double theta;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "F2_ABSCISSAS_AB - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    theta = ( double ) ( n - i ) * r8_pi 
          / ( double ) ( n + 1 );
    x[i] = 0.5 * ( ( b + a ) + ( b - a ) * cos ( theta ) );
  }

  return x;
}
//****************************************************************************80

double *interp_lagrange ( int m, int data_num, double t_data[], 
  double p_data[], int interp_num, double t_interp[] )

//****************************************************************************80
//
//  Purpose:
//
//    INTERP_LAGRANGE: Lagrange polynomial interpolant to a curve in M dimensions.
//
//  Discussion:
//
//    From a space of M dimensions, we are given a sequence of
//    DATA_NUM points, which are presumed to be successive samples
//    from a curve of points P.
//
//    We are also given a parameterization of this data, that is,
//    an associated sequence of DATA_NUM values of a variable T.
//
//    Thus, we have a sequence of values P(T), where T is a scalar,
//    and each value of P is of dimension M.
//
//    We are then given INTERP_NUM values of T, for which values P
//    are to be produced, by linear interpolation of the data we are given.
//
//    The user may request extrapolation.  This occurs whenever
//    a T_INTERP value is less than the minimum T_DATA or greater than the
//    maximum T_DATA.  In that case, extrapolation is used.
//
//    For each spatial component, a polynomial of degree
//    ( DATA_NUM - 1 ) is generated for the interpolation.  In most cases,
//    such a polynomial interpolant begins to oscillate as DATA_NUM
//    increases, even if the original data seems well behaved.  Typically,
//    values of DATA_NUM should be no greater than 10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Input, double T_DATA[DATA_NUM], the value of the
//    independent variable at the sample points.
//
//    Input, double P_DATA[M*DATA_NUM], the value of the
//    dependent variables at the sample points.
//
//    Input, int INTERP_NUM, the number of points
//    at which interpolation is to be done.
//
//    Input, double T_INTERP[INTERP_NUM], the value of the
//    independent variable at the interpolation points.
//
//    Output, double INTERP_LAGRANGE[M*DATA_NUM], the interpolated
//    values of the dependent variables at the interpolation points.
//
{
  double *l_interp;
  double *p_interp;
//
//  Evaluate the DATA_NUM Lagrange polynomials associated with T_DATA(1:DATA_NUM)
//  for the interpolation points T_INTERP(1:INTERP_NUM).
//
  l_interp = lagrange_value ( data_num, t_data, interp_num, t_interp );
//
//  Multiply P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
//  to get P_INTERP(1:M,1:INTERP_NUM).
//
  p_interp = r8mat_mm_new ( m, data_num, interp_num, p_data, l_interp );

  delete [] l_interp;

  return p_interp;
}
//****************************************************************************80

double *interp_linear ( int m, int data_num, double t_data[], double p_data[], 
  int interp_num, double t_interp[] )

//****************************************************************************80
//
//  Purpose:
//
//    INTERP_LINEAR: piecewise linear interpolation to a curve in M dimensions.
//
//  Discussion:
//
//    From a space of M dimensions, we are given a sequence of
//    DATA_NUM points, which are presumed to be successive samples
//    from a curve of points P.
//
//    We are also given a parameterization of this data, that is,
//    an associated sequence of DATA_NUM values of a variable T.
//    The values of T are assumed to be strictly increasing.
//
//    Thus, we have a sequence of values P(T), where T is a scalar,
//    and each value of P is of dimension M.
//
//    We are then given INTERP_NUM values of T, for which values P
//    are to be produced, by linear interpolation of the data we are given.
//
//    Note that the user may request extrapolation.  This occurs whenever
//    a T_INTERP value is less than the minimum T_DATA or greater than the
//    maximum T_DATA.  In that case, linear extrapolation is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Input, double T_DATA[DATA_NUM], the value of the
//    independent variable at the sample points.  The values of T_DATA
//    must be strictly increasing.
//
//    Input, double P_DATA[M*DATA_NUM], the value of the
//    dependent variables at the sample points.
//
//    Input, int INTERP_NUM, the number of points
//    at which interpolation is to be done.
//
//    Input, double T_INTERP[INTERP_NUM], the value of the
//    independent variable at the interpolation points.
//
//    Output, double INTERP_LINEAR[M*DATA_NUM], the interpolated
//    values of the dependent variables at the interpolation points.
//
{
  int i;
  int interp;
  int left;
  double *p_interp;
  int right;
  double t;

  if ( ! r8vec_ascends_strictly ( data_num, t_data ) )
  {
    cerr << "\n";
    cerr << "INTERP_LINEAR - Fatal error\n";
    cerr << 
      "  Independent variable array T_DATA is not strictly increasing.\n";
    exit ( 1 );
  }

  p_interp = new double[m*interp_num];

  for ( interp = 0; interp < interp_num; interp++ )
  {
    t = t_interp[interp];
//
//  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
//  nearest to, TVAL.
//
    r8vec_bracket0 ( data_num, t_data, t, left, right );

    for ( i = 0; i < m; i++ )
    {
      p_interp[i+interp*m] = 
        ( ( t_data[right] - t                ) * p_data[i+left*m]   
        + (                 t - t_data[left] ) * p_data[i+right*m] ) 
        / ( t_data[right]     - t_data[left] );
    }
  }

  return p_interp;
}
//****************************************************************************80

double *interp_nearest ( int m, int data_num, double t_data[], double p_data[], 
  int interp_num, double t_interp[] )

//****************************************************************************80
//
//  Purpose:
//
//    INTERP_NEAREST: Nearest neighbor interpolation to a curve in M dimensions.
//
//  Discussion:
//
//    From a space of M dimensions, we are given a sequence of
//    DATA_NUM points, which are presumed to be successive samples
//    from a curve of points P.
//
//    We are also given a parameterization of this data, that is,
//    an associated sequence of DATA_NUM values of a variable T.
//
//    Thus, we have a sequence of values P(T), where T is a scalar,
//    and each value of P is of dimension M.
//
//    We are then given INTERP_NUM values of T, for which values P
//    are to be produced, by nearest neighbor interpolation.
//
//    The user may request extrapolation.  This occurs whenever
//    a T_INTERP value is less than the minimum T_DATA or greater than the
//    maximum T_DATA.  In that case, extrapolation is used.
//
//    The resulting interpolant is piecewise constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Input, double T_DATA[DATA_NUM], the value of the
//    independent variable at the sample points.
//
//    Input, double P_DATA[M*DATA_NUM], the value of the
//    dependent variables at the sample points.
//
//    Input, int INTERP_NUM, the number of points
//    at which interpolation is to be done.
//
//    Input, double T_INTERP[INTERP_NUM], the value of the
//    independent variable at the interpolation points.
//
//    Output, double INTERP_NEAREST[M*DATA_NUM], the interpolated
//    values of the dependent variables at the interpolation points.
//
{
  int i;
  int jd;
  int ji;
  double *p_interp;
//
//  For each interpolation point, find the index of the nearest data point.
//
  p_interp = new double[m*interp_num];

  for ( ji = 0; ji < interp_num; ji++ )
  {
    jd = r8vec_sorted_nearest0 ( data_num, t_data, t_interp[ji] );
    for ( i = 0; i < m; i++ )
    {
      p_interp[i+ji*m] = p_data[i+jd*m];
    }
  }

  return p_interp;
}
//****************************************************************************80

double *lagrange_value ( int data_num, double t_data[], int interp_num, 
  double t_interp[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAGRANGE_VALUE evaluates the Lagrange polynomials.
//
//  Discussion:
//
//    Given DATA_NUM distinct abscissas, T_DATA(1:DATA_NUM),
//    the I-th Lagrange polynomial L(I)(T) is defined as the polynomial of
//    degree DATA_NUM - 1 which is 1 at T_DATA(I) and 0 at the DATA_NUM - 1
//    other abscissas.
//
//    A formal representation is:
//
//      L(I)(T) = Product ( 1 <= J <= DATA_NUM, I /= J )
//       ( T - T(J) ) / ( T(I) - T(J) )
//
//    This routine accepts a set of INTERP_NUM values at which all the Lagrange
//    polynomials should be evaluated.
//
//    Given data values P_DATA at each of the abscissas, the value of the
//    Lagrange interpolating polynomial at each of the interpolation points
//    is then simple to compute by matrix multiplication:
//
//      P_INTERP(1:INTERP_NUM) =
//        P_DATA(1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
//
//    or, in the case where P is multidimensional:
//
//      P_INTERP(1:M,1:INTERP_NUM) =
//        P_DATA(1:M,1:DATA_NUM) * L_INTERP(1:DATA_NUM,1:INTERP_NUM)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the number of data points.
//    DATA_NUM must be at least 1.
//
//    Input, double T_DATA[DATA_NUM], the data points.
//
//    Input, int INTERP_NUM, the number of
//    interpolation points.
//
//    Input, double T_INTERP[INTERP_NUM], the
//    interpolation points.
//
//    Output, double LAGRANGE_VALUE[DATA_NUM*INTERP_NUM], the values
//    of the Lagrange polynomials at the interpolation points.
//
{
  int i;
  int i1;
  int i2;
  int j;
  double *l_interp;

  l_interp = new double[data_num*interp_num];
//
//  Evaluate the polynomial.
//
  for ( j = 0; j < interp_num; j++ )
  {
    for ( i = 0; i < data_num; i++ )
    {
      l_interp[i+j*data_num] = 1.0;
    }
  }

  for ( i1 = 0; i1 < data_num; i1++ )
  {
    for ( i2 = 0; i2 < data_num; i2++ )
    {
      if ( i1 != i2 )
      {
        for ( j = 0; j < interp_num; j++ )
        {
          l_interp[i1+j*data_num] = l_interp[i1+j*data_num] 
            * ( t_interp[j] - t_data[i2] ) / ( t_data[i1] - t_data[i2] );
        }
      }
    }
  }

  return l_interp;
}
//****************************************************************************80

double *ncc_abscissas ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_ABSCISSAS computes the Newton Cotes Closed abscissas.
//
//  Discussion:
//
//    The interval is [ -1, 1 ].
//
//    The abscissas are the equally spaced points between -1 and 1,
//    including the endpoints.
//
//    If N is 1, however, the single abscissas is X = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double NCC_ABSCISSAS[N], the abscissas.
//
{
  int i;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "NCC_ABSCISSAS - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.0;
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * ( -1.0 )   
           + ( double ) (     i     ) * ( +1.0 ) ) 
           / ( double ) ( n     - 1 ); 
  }

  return x;
}
//****************************************************************************80

double *ncc_abscissas_ab ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    NCC_ABSCISSAS_AB computes the Newton Cotes Closed abscissas for [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the order of the rule.
//
//    Output, double NCC_ABSCISSAS_AB[N], the abscissas.
//
{
  int i;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "NCC_ABSCISSAS_AB - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  if ( n == 1 )
  {
    x[0] = 0.5 * ( b + a );
    return x;
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i - 1 ) * ( a )   
           + ( double ) (     i     ) * ( b ) ) 
           / ( double ) ( n     - 1 ); 
  }

  return x;
}
//****************************************************************************80

double *nco_abscissas ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_ABSCISSAS computes the Newton Cotes Open abscissas.
//
//  Discussion:
//
//    The interval is [ -1, 1 ].
//
//    The abscissas are the equally spaced points between -1 and 1,
//    not including the endpoints.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the rule.
//
//    Output, double NCO_ABSCISSAS[N], the abscissas.
//
{
  int i;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "NCO_ABSCISSAS - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * ( -1.0 )   
           + ( double ) (     i + 1 ) * ( +1.0 ) ) 
           / ( double ) ( n     + 1 );
  }

  return x;
}
//****************************************************************************80

double *nco_abscissas_ab ( double a, double b, int n )

//****************************************************************************80
//
//  Purpose:
//
//    NCO_ABSCISSAS_AB computes the Newton Cotes Open abscissas for [A,B].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the endpoints of the interval.
//
//    Input, int N, the order of the rule.
//
//    Output, double NCO_ABSCISSAS_AB[N], the abscissas.
//
{
  int i;
  double *x;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "NCO_ABSCISSAS_AB - Fatal error\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( double ) ( n - i     ) * ( a )   
           + ( double ) (     i + 1 ) * ( b ) ) 
           / ( double ) ( n     + 1 );
  }

  return x;
}
//****************************************************************************80

double *parameterize_arc_length ( int m, int data_num, double p_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    PARAMETERIZE_ARC_LENGTH parameterizes data by pseudo-arclength.
//
//  Discussion:
//
//    A parameterization is required for the interpolation.
//
//    This routine provides a parameterization by computing the
//    pseudo-arclength of the data, that is, the Euclidean distance
//    between successive points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Input, double P_DATA[M*DATA_NUM], the data values.
//
//    Output, double PARAMETERIZE_ARC_LENGTH[DATA_NUM], parameter values
//    assigned to the data.
//
{
  int i;
  int j;
  double t;
  double *t_data;

  t_data = new double[data_num];

  t_data[0] = 0.0;
  for ( j = 1; j < data_num; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m; i++ )
    {
      t = t + pow ( p_data[i+j*m] - p_data[i+(j-1)*m], 2 );
    }
    t_data[j] = t_data[j-1] + sqrt ( t ); 
  }

  return t_data;
}
//****************************************************************************80

double *parameterize_index ( int m, int data_num, double p_data[] )

//****************************************************************************80
//
//  Purpose:
//
//    PARAMETERIZE_INDEX parameterizes data by its index.
//
//  Discussion:
//
//    A parameterization is required for the interpolation.
//
//    This routine provides a naive parameterization by vector index.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 March 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int DATA_NUM, the number of data points.
//
//    Input, double P_DATA[M*DATA_NUM], the data values.
//
//    Output, double PARAMETERIZE_INDEX[DATA_NUM], parameter values
//    assigned to the data.
//
{
  int j;
  double *t_data;

  t_data = new double[data_num];

  for ( j = 0; j < data_num; j++ )
  {
    t_data[j] = ( double ) j;
  }

  return t_data;
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
    value = + x;
  }
  else
  {
    value = - x;
  }
  return value;
}
//****************************************************************************80

double *r8mat_expand_linear2 ( int m, int n, double a[], int m2, int n2 )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_EXPAND_LINEAR2 expands an R8MAT by linear interpolation.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    In this version of the routine, the expansion is indicated
//    by specifying the dimensions of the expanded array.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns in A.
//
//    Input, double A(M,N), a "small" M by N array.
//
//    Input, int M2, N2, the number of rows and columns in A2.
//
//    Output, double R8MAT_EXPAND_LINEAR2[M2*N2], the expanded array,
//    which contains an interpolated version of the data in A.
//
{
  double *a2;
  int i;
  int i1;
  int i2;
  int j;
  int j1;
  int j2;
  double r;
  double r1;
  double r2;
  double s;
  double s1;
  double s2;

  a2 = new double[m2*n2];

  for ( i = 1; i <= m2; i++ )
  {
    if ( m2 == 1 )
    {
      r = 0.5;
    }
    else
    {
      r = ( double ) ( i - 1 ) / ( double ) ( m2 - 1 );
    }

    i1 = 1 + ( int ) ( r * ( double ) ( m - 1 ) );
    i2 = i1 + 1;

    if ( m < i2 )
    {
      i1 = m - 1;
      i2 = m;
    }

    r1 = ( double ) ( i1 - 1 ) / ( double ) ( m - 1 );
    r2 = ( double ) ( i2 - 1 ) / ( double ) ( m - 1 );

    for ( j = 1; j <= n2; j++ )
    {
      if ( n2 == 1 )
      {
        s = 0.5;
      }
      else
      {
        s = ( double ) ( j - 1 ) / ( double ) ( n2 - 1 );
      }

      j1 = 1 + ( int ) ( s * ( double ) ( n - 1 ) );
      j2 = j1 + 1;

      if ( n < j2 )
      {
        j1 = n - 1;
        j2 = n;
      }

      s1 = ( double ) ( j1 - 1 ) / ( double ) ( n - 1 );
      s2 = ( double ) ( j2 - 1 ) / ( double ) ( n - 1 );

      a2[i-1+(j-1)*m2] =
        ( ( r2 - r ) * ( s2 - s ) * a[i1-1+(j1-1)*m]
        + ( r - r1 ) * ( s2 - s ) * a[i2-1+(j1-1)*m]
        + ( r2 - r ) * ( s - s1 ) * a[i1-1+(j2-1)*m]
        + ( r - r1 ) * ( s - s1 ) * a[i2-1+(j2-1)*m] )
        / ( ( r2 - r1 ) * ( s2 - s1 ) );
    }
  }

  return a2;
}
//****************************************************************************80

double *r8mat_mm_new ( int n1, int n2, int n3, double a[], double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
//
{
  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for ( i = 0; i < n1; i++ )
  {
    for ( j = 0; j < n3; j++ )
    {
      c[i+j*n1] = 0.0;
      for ( k = 0; k < n2; k++ )
      {
        c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}
//****************************************************************************80

bool r8vec_ascends_strictly ( int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ASCENDS_STRICTLY determines if an R8VEC is strictly ascending.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    Notice the effect of entry number 6 in the following results:
//
//      X = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.4, 9.8 )
//      Y = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.5, 9.8 )
//      Z = ( -8.1, 1.3, 2.2, 3.4, 7.5, 7.6, 9.8 )
//
//      R8VEC_ASCENDS_STRICTLY ( X ) = FALSE
//      R8VEC_ASCENDS_STRICTLY ( Y ) = FALSE
//      R8VEC_ASCENDS_STRICTLY ( Z ) = TRUE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the size of the array.
//
//    Input, double X[N], the array to be examined.
//
//    Output, bool R8VEC_ASCENDS_STRICTLY, is TRUE if the
//    entries of X strictly ascend.
//
{
  int i;
  bool value;

  for ( i = 0; i < n - 1; i++ )
  {
    if ( x[i+1] <= x[i] )
    {
      value = false;
      return value;
    }
  }
  value = true;

  return value;
}
//****************************************************************************80

void r8vec_bracket0 ( int n, double x[], double xval, int &left,
  int &right )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET0 searches a sorted array for successive brackets of a value.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    If the values in the vector are thought of as defining intervals
//    on the real line, then this routine searches for the interval
//    nearest to or containing the given value.
//
//    It is always true that RIGHT = LEFT+1.
//
//    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
//      XVAL   < X[0] < X[1];
//    If X(1) <= XVAL < X[N-1], then
//      X[LEFT-1] <= XVAL < X[RIGHT-1];
//    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
//      X[LEFT-1] <= X[RIGHT-1] <= XVAL.
//
//    For consistency, this routine computes indices RIGHT and LEFT
//    that are 1-based, although it would be more natural in C and
//    C++ to use 0-based values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, length of input array.
//
//    Input, double X[N], an array that has been sorted into ascending order.
//
//    Input, double XVAL, a value to be bracketed.
//
//    Output, int &LEFT, &RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      left = i - 2;
      right = i - 1;
      return;
    }

   }

  left = n - 2;
  right = n - 1;

  return;
}
//****************************************************************************80

double *r8vec_expand_linear ( int n, double x[], int fat )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EXPAND_LINEAR linearly interpolates new data into an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of input data values.
//
//    Input, double X[N], the original data.
//
//    Input, int FAT, the number of data values to interpolate
//    between each pair of original data values.
//
//    Output, double R8VEC_EXPAND_LINEAR[(N-1)*(FAT+1)+1], the "fattened" data.
//
{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = new double[(n-1)*(fat+1)+1];

  k = 0;

  for ( i = 0; i < n-1; i++ )
  {
    xfat[k] = x[i];
    k = k + 1;

    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] )
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;

  return xfat;
}
//****************************************************************************80

double *r8vec_expand_linear2 ( int n, double x[], int before, int fat, 
  int after )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EXPAND_LINEAR2 linearly interpolates new data into an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    This routine starts with a vector of data.
//
//    The intent is to "fatten" the data, that is, to insert more points
//    between successive values of the original data.
//
//    There will also be extra points placed BEFORE the first original
//    value and AFTER that last original value.
//
//    The "fattened" data is equally spaced between the original points.
//
//    The BEFORE data uses the spacing of the first original interval,
//    and the AFTER data uses the spacing of the last original interval.
//
//  Example:
//
//    N = 3
//    BEFORE = 3
//    FAT = 2
//    AFTER = 1
//
//    X    = (/                   0.0,           6.0,             7.0       /)
//    XFAT = (/ -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 6.33, 6.66, 7.0, 7.66 /)
//            3 "BEFORE's"        Old  2 "FATS"  Old    2 "FATS"  Old  1 "AFTER"
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of input data values.
//    N must be at least 2.
//
//    Input, double X[N], the original data.
//
//    Input, int BEFORE, the number of "before" values.
//
//    Input, int FAT, the number of data values to interpolate
//    between each pair of original data values.
//
//    Input, int AFTER, the number of "after" values.
//
//    Output, double R8VEC_EXPAND_LINEAR2[BEFORE+(N-1)*(FAT+1)+1+AFTER], the
//    "fattened" data.
//
{
  int i;
  int j;
  int k;
  double *xfat;

  xfat = new double[before+(n-1)*(fat+1)+1+after];

  k = 0;
//
//  Points BEFORE.
//
  for ( j = 1 - before + fat; j <= fat; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * ( x[0] - ( x[1] - x[0] ) ) 
              + ( double ) (       j     ) *   x[0]          ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }
//
//  Original points and FAT points.
//
  for ( i = 0; i < n - 1; i++ )
  {
    xfat[k] = x[0];
    k = k + 1;
    for ( j = 1; j <= fat; j++ )
    {
      xfat[k] = ( ( double ) ( fat - j + 1 ) * x[i]
                + ( double ) (       j     ) * x[i+1] ) 
                / ( double ) ( fat     + 1 );
      k = k + 1;
    }
  }

  xfat[k] = x[n-1];
  k = k + 1;
//
//  Points AFTER.
//
  for ( j = 1; j <= after; j++ )
  {
    xfat[k] = ( ( double ) ( fat - j + 1 ) * x[n-1]
              + ( double ) (       j     ) * ( x[n-1] + ( x[n-1] - x[n-2] ) ) ) 
              / ( double ) ( fat     + 1 );
    k = k + 1;
  }

  return xfat;
}
//****************************************************************************80

int r8vec_sorted_nearest0 ( int n, double a[], double value )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORTED_NEAREST0 returns the nearest element in a sorted R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Input, double A[N], a sorted vector.
//
//    Input, double VALUE, the value whose nearest vector entry is sought.
//
//    Output, int R8VEC_SORTED_NEAREST0, the index of the nearest
//    entry in the vector.
//
{
  int hi;
  int lo;
  int mid;

  if ( n < 1 )
  {
    return (-1);
  }

  if ( n == 1 )
  {
    return 0;
  }

  if ( a[0] < a[n-1] )
  {
    if ( value < a[0] )
    {
      return 0;
    }
    else if ( a[n-1] < value )
    {
      return n - 1;
    }
//
//  Seek an interval containing the value.
//
    lo = 1;
    hi = n;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid - 1;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid;
      }
      else
      {
        lo = mid;
      }
    }
//
//  Take the nearest.
//
    if ( r8_abs ( value - a[lo-1] ) < r8_abs ( value - a[hi-1] ) )
    {
      return lo - 1;
    }
    else
    {
      return hi - 1;
    }
  }
//
//  A descending sorted vector A.
//
  else
  {
    if ( value < a[n-1] )
    {
      return n - 1;
    }
    else if ( a[0] < value )
    {
      return 0;
    }
//
//  Seek an interval containing the value.
//
    lo = n;
    hi = 1;

    while ( lo < hi - 1 )
    {
      mid = ( lo + hi ) / 2;

      if ( value == a[mid-1] )
      {
        return mid - 1;
      }
      else if ( value < a[mid-1] )
      {
        hi = mid - 1;
      }
      else
      {
        lo = mid - 1;
      }
    }
//
//  Take the nearest.
//
    if ( r8_abs ( value - a[lo-1] ) < r8_abs ( value - a[hi-1] ) )
    {
      return lo - 1;
    }
    else
    {
      return hi - 1;
    }
  }
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
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
