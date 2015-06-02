# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "triangle_felippa_rule.hpp"

//****************************************************************************80

void comp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    COMP_NEXT computes the compositions of the integer N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to N.  The compositions (1,2,1)
//    and (1,1,2) are considered to be distinct.
//
//    The routine computes one composition on each call until there are no more.
//    For instance, one composition of 6 into 3 parts is
//    3+2+1, another would be 6+0+0.
//
//    On the first call to this routine, set MORE = FALSE.  The routine
//    will compute the first element in the sequence of compositions, and
//    return it, as well as setting MORE = TRUE.  If more compositions
//    are desired, call again, and again.  Each time, the routine will
//    return with a new composition.
//
//    However, when the LAST composition in the sequence is computed 
//    and returned, the routine will reset MORE to FALSE, signaling that
//    the end of the sequence has been reached.
//
//    This routine originally used a SAVE statement to maintain the
//    variables H and T.  I have decided that it is safer
//    to pass these variables as arguments, even though the user should
//    never alter them.  This allows this routine to safely shuffle
//    between several ongoing calculations.
//
//
//    There are 28 compositions of 6 into three parts.  This routine will
//    produce those compositions in the following order:
//
//     I         A
//     -     ---------
//     1     6   0   0
//     2     5   1   0
//     3     4   2   0
//     4     3   3   0
//     5     2   4   0
//     6     1   5   0
//     7     0   6   0
//     8     5   0   1
//     9     4   1   1
//    10     3   2   1
//    11     2   3   1
//    12     1   4   1
//    13     0   5   1
//    14     4   0   2
//    15     3   1   2
//    16     2   2   2
//    17     1   3   2
//    18     0   4   2
//    19     3   0   3
//    20     2   1   3
//    21     1   2   3
//    22     0   3   3
//    23     2   0   4
//    24     1   1   4
//    25     0   2   4
//    26     1   0   5
//    27     0   1   5
//    28     0   0   6
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the integer whose compositions are desired.
//
//    Input, int K, the number of parts in the composition.
//
//    Input/output, int A[K], the parts of the composition.
//
//    Input/output, bool *MORE.
//    Set MORE = FALSE on first call.  It will be reset to TRUE on return
//    with a new composition.  Each new call returns another composition until
//    MORE is set to FALSE when the last composition has been computed
//    and returned.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;

  if ( !( *more ) )
  {
    *t = n;
    *h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < *t )
    {
      *h = 0;
    }
    *h = *h + 1;
    *t = a[*h-1];
    a[*h-1] = 0;
    a[0] = *t - 1;
    a[*h] = a[*h] + 1;
  }

  *more = ( a[k-1] != n );

  return;
}
//****************************************************************************80

double *monomial_value ( int m, int n, int e[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    This routine evaluates a monomial of the form
//
//      product ( 1 <= i <= m ) x(i)^e(i)
//
//    The combination 0.0^0 is encountered is treated as 1.0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of evaluation points.
//
//    Input, int E[M], the exponents.
//
//    Input, double X[M*N], the point coordinates.
//
//    Output, double MONOMIAL_VALUE[N], the monomial values.
//
{
  int i;
  int j;
  double *v;

  v = new double[n];
  for ( j = 0; j < n; j++)
  {
    v[j] = 1.0;
  }
//v = r8vec_ones_new ( n );

  for ( i = 0; i < m; i++ )
  {
    if ( 0 != e[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        v[j] = v[j] * pow ( x[i+j*m], e[i] );
      }
    }
  }

  return v;
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

void subcomp_next ( int n, int k, int a[], bool *more, int *h, int *t )

//****************************************************************************80
//
//  Purpose:
//
//    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
//
//  Discussion:
//
//    A composition of the integer N into K parts is an ordered sequence
//    of K nonnegative integers which sum to a value of N.
//
//    A subcomposition of the integer N into K parts is a composition
//    of M into K parts, where 0 <= M <= N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 July 2008
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the integer whose subcompositions are desired.
//
//    Input, int K, the number of parts in the subcomposition.
//
//    Input/output, int A[K], the parts of the subcomposition.
//
//    Input/output, bool *MORE, set by the user to start the computation,
//    and by the routine to terminate it.
//
//    Input/output, int *H, *T, two internal parameters needed for the
//    computation.  The user should allocate space for these in the calling
//    program, include them in the calling sequence, but never alter them!
//
{
  int i;
  static bool more2 = false;
  static int n2 = 0;
//
//  The first computation.
//
  if ( !( *more ) )
  {
    n2 = 0;

    for ( i = 0; i < k; i++ )
    {
      a[i] = 0;
    }
    more2 = false;
    *h = 0;
    *t = 0;

    *more = true;
  }
//
//  Do the next element at the current value of N.
//
  else if ( more2 )
  {
    comp_next ( n2, k, a, &more2, h, t );
  }
  else
  {
    more2 = false;
    n2 = n2 + 1;

    comp_next ( n2, k, a, &more2, h, t );
  }
//
//  Termination occurs if MORE2 = FALSE and N2 = N.
//
  if ( !more2 && n2 == n )
  {
    *more = false;
  }

  return;
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
//****************************************************************************80

double triangle_unit_monomial ( int expon[2] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_MONOMIAL integrates a monomial over the unit triangle.
//
//  Discussion:
//
//    This routine integrates a monomial of the form
//
//      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//    Integral ( over unit triangle ) x^m y^n dx dy = m! * n! / ( m + n + 2 )!
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON[2], the exponents.
//
//    Output, double TRIANGLE_UNIT_MONOMIAL, the integral of the monomial.
//
{
  int i;
  int k;
  double value;
//
//  The first computation ends with VALUE = 1.0;
//
  value = 1.0;

// k = 0;
//
//  The first loop simply computes 1 so we short circuit it!
//
// for ( i = 1; i <= expon[0]; i++ )
// {
//   k = k + 1;
//   value = value * ( double ) ( i ) / ( double ) ( k );
// }

  k = expon[0];

  for ( i = 1; i <= expon[1]; i++ )
  {
    k = k + 1;
    value = value * ( double ) ( i ) / ( double ) ( k );
  }

  k = k + 1;
  value = value / ( double ) ( k );

  k = k + 1;
  value = value / ( double ) ( k );

  return value;
}
//****************************************************************************80

void triangle_unit_o01 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O01 returns a 1 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 1.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[1], the weights.
//
//    Output, double XY[2*1], the abscissas.
//
{
  int order = 1;

  double w_save[1] = { 
    1.0 };
  double xy_save[2*1] = { 
    0.33333333333333333333,  0.33333333333333333333 }; 

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o03 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O03 returns a 3 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 2.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[3], the weights.
//
//    Output, double XY[2*3], the abscissas.
//
{
  int order = 3;

  double w_save[3] = { 
    0.33333333333333333333, 
    0.33333333333333333333, 
    0.33333333333333333333 };
  double xy_save[2*3] = { 
    0.66666666666666666667,  0.16666666666666666667, 
    0.16666666666666666667,  0.66666666666666666667, 
    0.16666666666666666667,  0.16666666666666666667 }; 

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o03b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O03B returns a 3 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 2.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[3], the weights.
//
//    Output, double XY[2*3], the abscissas.
//
{
  int order = 3;

  double w_save[3] = { 
    0.33333333333333333333, 
    0.33333333333333333333, 
    0.33333333333333333333 };
  double xy_save[2*3] = { 
    0.0,  0.5, 
    0.5,  0.0, 
    0.5,  0.5 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o06 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O06 returns a 6 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 4.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[6], the weights.
//
//    Output, double XY[2*6], the abscissas.
//
{
  int order = 6;

  double w_save[6] = { 
    0.22338158967801146570, 
    0.22338158967801146570, 
    0.22338158967801146570, 
    0.10995174365532186764, 
    0.10995174365532186764, 
    0.10995174365532186764 };
  double xy_save[2*6] = { 
    0.10810301816807022736,  0.44594849091596488632, 
    0.44594849091596488632,  0.10810301816807022736, 
    0.44594849091596488632,  0.44594849091596488632, 
    0.81684757298045851308,  0.091576213509770743460, 
    0.091576213509770743460,  0.81684757298045851308, 
    0.091576213509770743460,  0.091576213509770743460 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o06b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O06B returns a 6 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 3.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[6], the weights.
//
//    Output, double XY[2*6], the abscissas.
//
{
  int order = 6;

  double w_save[6] = { 
    0.30000000000000000000, 
    0.30000000000000000000, 
    0.30000000000000000000, 
    0.033333333333333333333, 
    0.033333333333333333333, 
    0.033333333333333333333 };
  double xy_save[2*6] = { 
    0.66666666666666666667,  0.16666666666666666667, 
    0.16666666666666666667,  0.66666666666666666667, 
    0.16666666666666666667,  0.16666666666666666667, 
    0.0,  0.5, 
    0.5,  0.0, 
    0.5,  0.5 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o07 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O07 returns a 7 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 5.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 April 2009
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
//    Output, double W[7], the weights.
//
//    Output, double XY[2*7], the abscissas.
//
{
  int order = 7;

  double w_save[7] = { 
    0.12593918054482715260, 
    0.12593918054482715260, 
    0.12593918054482715260, 
    0.13239415278850618074, 
    0.13239415278850618074, 
    0.13239415278850618074, 
    0.22500000000000000000 };
  double xy_save[2*7] = { 
    0.79742698535308732240,  0.10128650732345633880, 
    0.10128650732345633880,  0.79742698535308732240, 
    0.10128650732345633880,  0.10128650732345633880, 
    0.059715871789769820459,  0.47014206410511508977, 
    0.47014206410511508977,  0.059715871789769820459, 
    0.47014206410511508977,  0.47014206410511508977, 
    0.33333333333333333333,  0.33333333333333333333 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

void triangle_unit_o12 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_O12 returns a 12 point quadrature rule for the unit triangle.
//
//  Discussion:
//
//    This rule is precise for monomials through degree 6.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2009
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
//    Output, double W[12], the weights.
//
//    Output, double XY[2*12], the abscissas.
//
{
  int order = 12;

  double w_save[12] = { 
     0.050844906370206816921,
     0.050844906370206816921,
     0.050844906370206816921,
     0.11678627572637936603,
     0.11678627572637936603,
     0.11678627572637936603,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194,
     0.082851075618373575194 };
  double xy_save[2*12] = { 
    0.87382197101699554332,  0.063089014491502228340, 
    0.063089014491502228340,  0.87382197101699554332, 
    0.063089014491502228340,  0.063089014491502228340, 
    0.50142650965817915742,  0.24928674517091042129, 
    0.24928674517091042129,  0.50142650965817915742, 
    0.24928674517091042129,  0.24928674517091042129, 
    0.053145049844816947353,  0.31035245103378440542, 
    0.31035245103378440542,  0.053145049844816947353, 
    0.053145049844816947353,  0.63650249912139864723, 
    0.31035245103378440542,  0.63650249912139864723, 
    0.63650249912139864723,  0.053145049844816947353, 
    0.63650249912139864723,  0.31035245103378440542 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 2*order, xy_save, xy );

  return;
}
//****************************************************************************80

double triangle_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
//
//  Discussion:
//
//    The "volume" of a triangle is usually called its area.
//
//    The integration region is:
//
//      0 <= X,
//      0 <= Y, 
//      X + Y <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TRIANGLE_UNIT_VOLUME, the volume.
//
{
  double volume;

  volume = 1.0 / 2.0;

  return volume;
}

