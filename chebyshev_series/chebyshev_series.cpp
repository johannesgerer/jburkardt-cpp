# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cassert>
# include <cmath>
# include <ctime>

using namespace std;

# include "chebyshev_series.hpp"

//****************************************************************************80

double echebser0 ( double x, double coef[], int nc )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEBSER0 evaluates a Chebyshev series.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ECHEBSER0, the value of the Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];

  x2 = 2.0 * x;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double echebser1 ( double x, double coef[], int nc, double &y1 )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEBSER1 evaluates a Chebyshev series and first derivative.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ECHEBSER1, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];
  c0 = coef[nc-1];

  x2 = 2.0 * x;

  for ( i = nc - 2; 0 <= i; i-- ) 
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;

    if ( 0 < i ) 
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
  }
  y1 = c0 - c2;

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double echebser2 ( double x, double coef[], int nc, double &y1, double &y2 )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEBSER2 evaluates a Chebyshev series and two derivatives.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ECHEBSER2, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
//    Output, double &Y2, the value of the 2nd derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double d0;
  double d1 = 0.0;
  double d2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];
  c0 = coef[nc-1];
  d0 = coef[nc-1];

  x2 = 2.0 * x;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;

    if ( 0 < i )
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
    if ( 1 < i ) 
    {
      d2 = d1;
      d1 = d0;
      d0 = c0 - d2 + x2 * d1;
    }
  }

  y2 = ( d0 - d2 ) * 4.0;
  y1 = c0 - c2;

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double echebser3 ( double x, double coef[], int nc, double &y1, double &y2,
   double &y3 )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEBSER3 evaluates a Chebyshev series and three derivatives.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ECHEBSER3, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
//    Output, double &Y2, the value of the 2nd derivative of the
//    Chebyshev series at X.
//
//    Output, double &Y3, the value of the 3rd derivative of the
//    Chebyshev series at X.
//
{
  assert 
  ( 
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double d0;
  double d1 = 0.0;
  double d2 = 0.0;
  double e0;
  double e1 = 0.0;
  double e2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];
  c0 = coef[nc-1];
  d0 = coef[nc-1];
  e0 = coef[nc-1];

  x2 = 2.0 * x;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;

    if ( 0 < i ) 
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
    if ( 1 < i )
    {
      d2 = d1;
      d1 = d0;
      d0 = c0 - d2 + x2 * d1;
    }
    if ( 2 < i )
    {
      e2 = e1;
      e1 = e0;
      e0 = d0 - e2 + x2 * e1;
    }
  }

  y3 = ( e0 - e2 ) * 24.0;
  y2 = ( d0 - d2 ) * 4.0;
  y1 = c0 - c2;

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double echebser4 ( double x, double coef[], int nc, double &y1, double &y2,
   double &y3, double &y4 )

//****************************************************************************80
//
//  Purpose:
//
//    ECHEBSER4 evaluates a Chebyshev series and four derivatives.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ECHEBSER3, the value of the Chebyshev series at X.
//
//    Output, double &Y1, &Y2, &Y3, &Y4, the value of the 1st derivative of the
//    Chebyshev series at X.
//
{
  assert 
  ( 
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double d0;
  double d1 = 0.0;
  double d2 = 0.0;
  double e0;
  double e1 = 0.0;
  double e2 = 0.0;
  double f0;
  double f1 = 0.0;
  double f2 = 0.0;
  int i;
  double value;
  double y0;

  b0 = coef[nc-1];
  c0 = coef[nc-1];
  d0 = coef[nc-1];
  e0 = coef[nc-1];
  f0 = coef[nc-1];

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + 2.0 * x * b1;

    if ( 0 < i ) 
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + 2.0 * x * c1;
    }
    if ( 1 < i )
    {
      d2 = d1;
      d1 = d0;
      d0 = c0 - d2 + 2.0 * x * d1;
    }
    if ( 2 < i )
    {
      e2 = e1;
      e1 = e0;
      e0 = d0 - e2 + 2.0 * x * e1;
    }
    if ( 3 < i )
    {
      f2 = f1;
      f1 = f0;
      f0 = e0 - f2 + 2.0 * x * f1;
    }
  }

  y0 = ( b0 - b2 )        / 2.0;
  y1 =   c0 - c2;
  y2 = ( d0 - d2 ) *  2.0 * 2.0;
  y3 = ( e0 - e2 ) *  6.0 * 4.0;
  y4 = ( f0 - f2 ) * 24.0 * 8.0;

  return y0;
}
//****************************************************************************80

double evenchebser0 ( double x, double coef[], int nc )

//****************************************************************************80
//
//  Purpose:
//
//    EVENCHEBSER0 evaluates an even Chebyshev series.
//
//  Discussion:
//
//    This function implements Clenshaw's modification of his
//    algorithm for even series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double EVENCHEBSER0, the value of the Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;
  }

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double evenchebser1 ( double x, double coef[], int nc, double &y1 )

//****************************************************************************80
//
//  Purpose:
//
//    EVENCHEBSER1 evaluates an even Chebyshev series and first derivative.
//
//  Discussion:
//
//    This function implements a modification and extension of 
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference 
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double EVENCHEBSER1, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];
  c0 = coef[nc-1];

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;
    if ( 0 < i )
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
  }

  y1 = ( c0 - c2 ) * 4.0 * x;

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double evenchebser2 ( double x, double coef[], int nc, double &y1, double &y2 )

//****************************************************************************80
//
//  Purpose:
//
//    EVENCHEBSER2 evaluates an even Chebyshev series and first two derivatives.
//
//  Discussion:
//
//    This function implements a modification and extension of
//    Maess's algorithm.  Table 6.5.1 on page 164 of the reference
//    gives an example for treating the first derivative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double EVENCHEBSER2, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
//    Output, double &Y2, the value of the 2nd derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc &&
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double d0;
  double d1 = 0.0;
  double d2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];
  c0 = coef[nc-1];
  d0 = coef[nc-1];

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;
    if ( 0 < i )
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
    if ( 1 < i )
    {
      d2 = d1;
      d1 = d0;
      d0 = c0 - d2 + x2 * d1;
    }
  }

  y2 = ( d0 - d2 ) * 64.0 * x * x + ( c0 - c2 ) * 4.0;
  y1 = ( c0 - c2 ) * 4.0 * x;

  value = 0.5 * ( b0 - b2 );

  return value;
}
//****************************************************************************80

double oddchebser0 ( double x, double coef[], int nc )

//****************************************************************************80
//
//  Purpose:
//
//    ODDCHEBSER0 evaluates an odd Chebyshev series.
//
//  Discussion:
//
//    This function implements Clenshaw's modification of  his algorithm
//    for odd series.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ODDCHEBSER0, the value of the Chebyshev series at X.
//
{
  assert
  (
    0 < nc && 
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  int i;
  double value;
  double x2;

  b0 = coef[nc-1];

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    b0 = coef[i] - b2 + x2 * b1;
  }

  value = ( b0 - b1 ) * x;

  return value;
}
//****************************************************************************80

double oddchebser1 ( double x, double coef[], int nc, double &y1 )

//****************************************************************************80
//
//  Purpose:
//
//    ODDCHEBSER1 evaluates an odd Chebyshev series and the first derivative.
//
//  Discussion:
//
//    This function implements a modification and extension of
//    Clenshaw's algorithm. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ODDCHEBSER1, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc &&
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double coefi;
  int i;
  double value;
  double x2;

  coefi = 2.0 * coef[nc-1];
  b0 = coefi;
  c0 = coefi;

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    coefi = 2.0 * coef[i] - coefi;
    b0 = coefi - b2 + x2 * b1;
    if ( 0 < i )
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
  }
  y1 = ( c0 - c2 ) * 4.0 * x * x + ( b0 - b2 ) * 0.5;
  value = ( b0 - b2 ) * 0.5 * x;

  return value;
}
//****************************************************************************80

double oddchebser2 ( double x, double coef[], int nc, double &y1, double &y2 )

//****************************************************************************80
//
//  Purpose:
//
//    ODDCHEBSER2 evaluates an odd Chebyshev series and first two derivatives.
//
//  Discussion:
//
//    This function implements a modification and extension of
//    Clenshaw's algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 January 2014
//
//  Author:
//
//    Manfred Zimmer
//
//  Reference:
//
//    Charles Clenshaw,
//    Mathematical Tables, Volume 5,
//    Chebyshev series for mathematical functions,
//    London, 1962.
//
//    Gerhard Maess,
//    Vorlesungen ueber Numerische Mathematik II, Analysis,
//    Berlin, Akademie_Verlag, 1984-1988,
//    ISBN: 978-3764318840,
//    LC: QA297.M325.  
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//    -1 <= X <= +1.
//
//    Input, double COEF[NC], the Chebyshev series.
//
//    Input, int NC, the number of terms in the Chebyshev series.
//    0 < NC.
//
//    Output, double ODDCHEBSER2, the value of the Chebyshev series at X.
//
//    Output, double &Y1, the value of the 1st derivative of the
//    Chebyshev series at X.
//
//    Output, double &Y2, the value of the 2nd derivative of the
//    Chebyshev series at X.
//
{
  assert
  (
    0 < nc &&
    fabs ( x ) <= 1.0
  );

  double b0;
  double b1 = 0.0;
  double b2 = 0.0;
  double c0;
  double c1 = 0.0;
  double c2 = 0.0;
  double d0;
  double d1 = 0.0;
  double d2 = 0.0;
  double coefi;
  int i;
  double value;
  double x2;

  coefi = 2.0 * coef[nc-1];
  b0 = coefi;
  c0 = coefi;
  d0 = coefi;

  x2 = 4.0 * x * x - 2.0;

  for ( i = nc - 2; 0 <= i; i-- )
  {
    b2 = b1;
    b1 = b0;
    coefi = 2.0 * coef[i] - coefi;
    b0 = coefi - b2 + x2 * b1;
    if ( 0 < i )
    {
      c2 = c1;
      c1 = c0;
      c0 = b0 - c2 + x2 * c1;
    }
    if ( 1 < i )
    {
      d2 = d1;
      d1 = d0;
      d0 = c0 - d2 + x2 * d1;
    }
  }
  value = ( b0 - b2 ) * 0.5 * x;
  
  x2 = x * x;
  y1 = ( c0 - c2 ) * 4.0 * x2  +  ( b0 - b2 ) * 0.5;
  y2 = (( d0 - d2 ) * 64.0 * x2 + ( c0 - c2 ) * 12.0) * x;

  return value;
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
  const struct tm *tm_ptr;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm_ptr = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}

