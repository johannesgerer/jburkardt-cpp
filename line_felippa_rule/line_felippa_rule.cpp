# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "line_felippa_rule.hpp"

//****************************************************************************80

double line_monomial ( double a, double b, int expon )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_MONOMIAL: monomial integral over a line segment in 1D.
//
//  Discussion:
//
//    This function returns the integral of X^EXPON.
//
//    The integration region is:
//    A <= X <= B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Input, double A, B, the lower and upper limits.
//
//    Input, int EXPON, the exponent of X.  The exponent must not be -1.
//
//    Output, double LINE_MONOMIAL, the integral of X^EXPON.
//
{
  double value;

  if ( expon == - 1 )
  {
    cerr << "\n";
    cerr << "LINE_MONOMIAL - Fatal error!\n";
    cerr << "  Exponent = -1 is not a legal input.\n";
    exit ( 1 );
  }
  
  value = ( pow ( b, expon + 1 ) - pow ( a, expon + 1 ) ) 
    / ( double ) ( expon + 1 );

  return value;
}
//****************************************************************************80

void line_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_MONOMIAL_TEST tests LINE_MONOMIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  double a = 0.0;
  double b = 1.0;
  int expon;
  double value;

  cout << "\n";
  cout << "LINE_MONOMIAL_TEST\n";
  cout << "  For a line segment in 1D,\n";
  cout << "  LINE_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^EXPON\n";
  cout << "\n";
  cout << "  Volume = " << line_volume ( a, b ) << "\n";
  cout << "\n";
  cout << "     EXPON      INTEGRAL\n";
  cout << "\n";

  for ( expon = 0; expon <= degree_max; expon++ )
  {
    value = line_monomial ( a, b, expon );

    cout << "  " << setw(8)  << expon
         << "  " << setw(14) << value << "\n";
  }

  return;
}
//****************************************************************************80

void line_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_QUAD_TEST tests the rules for a line segment in 1D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  double a = 0.0;
  double b = 1.0;
  int expon;
  int j;
  int order;
  double quad;
  double *v;
  double *w;
  double *x;

  cout << "\n";
  cout << "LINE_QUAD_TEST\n";
  cout << "  For a line segment in 1D,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  LINE_UNIT_O01, a 1 point rule.\n";
  cout << "  LINE_UNIT_O02, a 2 point rule.\n";
  cout << "  LINE_UNIT_O03, a 3 point rule.\n";
  cout << "  LINE_UNIT_O04, a 4 point rule.\n";
  cout << "  LINE_UNIT_O05, a 5 point rule.\n";

  for ( expon = 0; expon <= degree_max; expon++ )
  {
    cout << "\n";
    cout << "  Monomial exponent:   " << expon << "\n";
    cout << "\n";

    for ( order = 1; order <= 5; order++ )
    {
      v = new double[order];
      w = new double[order];
      x = new double[order];

      line_rule ( a, b, order, w, x );
      for ( j = 0; j < order; j++ )
      {
        v[j] = pow ( x[j], expon );
      }
      quad = r8vec_dot_product ( order, w, v );
      cout << setw(8) << order << "  "
           << setw(14) << quad << "\n";
      delete [] v;
      delete [] w;
      delete [] x;
    }
    cout << "\n";
    quad = line_monomial ( a, b, expon );
    cout << "   Exact  " << setw(14) << quad << "\n";
  }

  return;
}
//****************************************************************************80

void line_rule ( double a, double b, int order, double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_RULE returns a quadrature rule for a line segment in 1D.
//
//  Discussion:
//
//    The integration region is:
//      A <= X <= B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, double A, B, the lower and upper limits.
//
//    Input, int ORDER, the order of the rule.
//
//    Output, double W[ORDER], the weights.
//
//    Output, double X[ORDER], the abscissas.
//
{
  int j;

  if ( order == 1 )
  {
    line_unit_o01 ( w, x );
  }
  else if ( order == 2 )
  {
    line_unit_o02 ( w, x );
  }
  else if ( order == 3 )
  {
    line_unit_o03 ( w, x );
  }
  else if ( order == 4 )
  {
    line_unit_o04 ( w, x );
  }
  else if ( order == 5 )
  {
    line_unit_o05 ( w, x );
  }
  else
  {
    cerr << "\n";
    cerr << "LINE_RULE - Fatal error!\n";
    cerr << "  Illegal value of ORDER.\n";
    exit ( 1 );
  }
//
//  Transform from [-1,+1] to [A,B]
//
  for ( j = 0; j < order; j++ )
  {
    w[j] = w[j] * ( b - a ) / 2.0;
    x[j] = ( ( 1.0 - x[j] ) * a   
           + ( 1.0 + x[j] ) * b ) 
           /   2.0;
  }

  return;
}
//****************************************************************************80

void line_unit_o01 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_O01 returns a 1 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[1], the weights.
//
//    Output, double X[1], the abscissas.
//
{
  int i;
  int order = 1;
  double w_save[1] = { 2.0 };
  double x_save[1] = { 0.0 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_unit_o02 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_O02 returns a 2 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[2], the weights.
//
//    Output, double X[2], the abscissas.
//
{
  int i;
  int order = 2;
  double w_save[2] = { 
    1.0000000000000000000,
    1.0000000000000000000 };
  double x_save[2] = { 
    -0.57735026918962576451,
     0.57735026918962576451 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_unit_o03 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_O03 returns a 3 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[3], the weights.
//
//    Output, double X[3], the abscissas.
//
{
  int i;
  int order = 3;
  double w_save[3] = { 
    0.55555555555555555556,
    0.88888888888888888889,
    0.55555555555555555556 };
  double x_save[3] = { 
    -0.77459666924148337704,
     0.00000000000000000000,
     0.77459666924148337704 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_unit_o04 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_O04 returns a 4 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[4], the weights.
//
//    Output, double X[4], the abscissas.
//
{
  int i;
  int order = 4;
  double w_save[4] = { 
    0.34785484513745385737,
    0.65214515486254614263,
    0.65214515486254614263,
    0.34785484513745385737 };
  double x_save[4] = { 
    -0.86113631159405257522,
    -0.33998104358485626480,
     0.33998104358485626480,
     0.86113631159405257522 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

void line_unit_o05 ( double w[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_O05 returns a 5 point quadrature rule for the unit line.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Arthur Stroud,
//    Approximate Calculation of Multiple Integrals,
//    Prentice Hall, 1971,
//    ISBN: 0130438936,
//    LC: QA311.S85.
//
//  Parameters:
//
//    Output, double W[5], the weights.
//
//    Output, double X[5], the abscissas.
//
{
  int i;
  int order = 5;
  double w_save[5] = { 
    0.23692688505618908751,
    0.47862867049936646804,
    0.56888888888888888889,
    0.47862867049936646804,
    0.23692688505618908751 };
  double x_save[5] = { 
    -0.90617984593866399280,
    -0.53846931010568309104,
     0.00000000000000000000,
     0.53846931010568309104,
     0.90617984593866399280 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( order, x_save, x );

  return;
}
//****************************************************************************80

double line_volume ( double a, double b )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_VOLUME: volume of a line segment in 1D.
//
//  Discussion:
//
//    The integration region is:
//    A <= X <= B
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 September 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the lower and upper limits.
//
//    Output, double LINE_VOLUME, the volume of the line.
//
{
  double volume;

  volume = b - a;

  return volume;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Output, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double r8vec_dot_product ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
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

