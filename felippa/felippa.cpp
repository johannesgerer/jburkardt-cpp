# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "felippa.hpp"

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

double hexa_unit_monomial ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    HEXA_UNIT_MONOMIAL integrates a monomial over the unit hexahedron.
//
//  Discussion:
//
//    This routine integrates a monomial of the form
//
//      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
//
//    The combination 0^0 should be treated as 1.
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
//    - 1.0 <= Z <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON[3], the exponents.
//
//    Output, double HEXA_UNIT_MONOMIAL, the integral of the monomial.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 0; i < 3; i++ )
  {
    if ( ( expon[i] % 2 ) == 1 )
    {
      value = 0.0;
    }
    else if ( expon[i] == -1 )
    {
      cout << "\n";
      cout << "HEXA_UNIT_MONOMIAL - Fatal error!\n";
      cout << "  Exponent of -1 encountered.\n";
      exit ( 1 );
    }
    else
    {
      value = value * 2.0 / ( double ) ( expon[i] + 1 );
    }
  }
  return value;
}
//****************************************************************************80

void hexa_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    HEXA_UNIT_MONOMIAL_TEST tests HEXA_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[3];
  int gamma;
  double value;

  cout << "\n";
  cout << "HEXA_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit hexahedron,\n";
  cout << "  HEXA_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << hexa_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = hexa_unit_monomial ( expon );

        cout << "  " << setw(8)  << expon[0]
             << "  " << setw(8)  << expon[1]
             << "  " << setw(8)  << expon[2]
             << "  " << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void hexa_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    HEXA_UNIT_QUAD_TEST tests the rules for the unit hexahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2008
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
# define DIM_NUM 3

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  int k;
  bool more;
  int order;
  int order_1d[DIM_NUM];
  double quad;
  int t;
  double *v;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "HEXA_UNIT_QUAD_TEST\n";
  cout << "  For the unit hexahedron,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  HEXA_UNIT_RULE, which returns N1 by N2 by N3 point rules..\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( i4vec_odd_any ( dim_num, expon ) )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    for ( k = 1; k <= 5; k++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        order_1d[dim] = k;
      }
      order = i4vec_product ( dim_num, order_1d );
      w = new double[order];
      xyz = new double[dim_num*order];
      hexa_unit_rule ( order_1d, w, xyz );
      v = monomial_value ( dim_num, order, expon, xyz );
      quad = hexa_unit_volume ( ) * r8vec_dot_product ( order, w, v );
      cout << "  " << setw(6) << order_1d[0]
           << "  " << setw(6) << order_1d[1]
           << "  " << setw(6) << order_1d[2]
           << "  " << setw(14) << quad << "\n";
      delete [] v;
      delete [] w;
      delete [] xyz;
    }
//
//  Try a rule of mixed orders.
//
    order_1d[0] = 3;
    order_1d[1] = 5;
    order_1d[2] = 2;
    order = i4vec_product ( dim_num, order_1d );
    w = new double[order];
    xyz = new double[dim_num*order];
    hexa_unit_rule ( order_1d, w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = hexa_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order_1d[0]
         << "  " << setw(6) << order_1d[1]
         << "  " << setw(6) << order_1d[2]
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    cout << "\n";
    quad = hexa_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << "      "
         << "  " << "      "
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void hexa_unit_rule ( int order_1d[], double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    HEXA_UNIT_RULE returns a quadrature rule for the unit hexahedron.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
//    - 1.0 <= Z <= 1.0
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
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, int ORDER_1D[3], the order of the rule in 
//    each dimension.  1 <= ORDER_1D(I) <= 5.
//
//    Output, double W[ORDER_1D[0]*ORDER_1D[1]*ORDER_1D[2]], the weights.
//
//    Output, double XYZ[3*ORDER_1D[0]*ORDER_1D[1]*ORDER_1D[2]], the abscissas.
//
{
  int dim;
  int dim_num = 3;
  int order;
  double *w_1d;
  double *x_1d;

  order = i4vec_product ( dim_num, order_1d );

  for ( dim = 0; dim < dim_num; dim++ )
  {
    w_1d = new double[order_1d[dim]];
    x_1d = new double[order_1d[dim]];

    if ( order_1d[dim] == 1 )
    {
      line_unit_o01 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 2 )
    {
      line_unit_o02 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 3 )
    {
      line_unit_o03 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 4 )
    {
      line_unit_o04 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 5 )
    {
      line_unit_o05 ( w_1d, x_1d );
    }
    else
    {
      cout << "\n";
      cout << "HEXA_UNIT_RULE - Fatal error!\n";
      cout << "  Illegal value of ORDER_1D[*].\n";
      exit ( 1 );
    }

    r8vec_direct_product ( dim, order_1d[dim], x_1d, 
      dim_num, order, xyz );

    r8vec_direct_product2 ( dim, order_1d[dim], w_1d, 
      dim_num, order, w );

    delete [] w_1d;
    delete [] x_1d;
  }

  return;
}
//****************************************************************************80

double hexa_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    HEXA_UNIT_VOLUME: volume of a unit hexahedron.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
//    - 1.0 <= Z <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double HEXA_UNIT_VOLUME, the volume.
//
{
  double value;

  value = 8.0;

  return value;
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

bool i4vec_odd_any ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_ODD_ANY is TRUE if any entry of an I4VEC is odd.
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
//    17 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector.
//
//    Output, bool I4VEC_ODD_ANY, TRUE if any entry is odd.
//
{
  int i;
  bool value;

  value = false;

  for ( i = 0; i < n; i++ )
  {
    if ( ( a[i] % 2 ) == 1 )
    {
      value = true;
      return value;
    }
  }
  return value;
}
//****************************************************************************80

int i4vec_product ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRODUCT multiplies the entries of an I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    Input:
//
//      A = ( 1, 2, 3, 4 )
//
//    Output:
//
//      I4VEC_PRODUCT = 24
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, int A[N], the vector
//
//    Output, int I4VEC_PRODUCT, the product of the entries of A.
//
{
  int i;
  int product;

  product = 1;
  for ( i = 0; i < n; i++ )
  {
    product = product * a[i];
  }

  return product;
}
//****************************************************************************80

double line_unit_monomial ( int expon[] )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_MONOMIAL: monomial integral in a unit line.
//
//  Discussion:
//
//    This function returns the integral of X^ALPHA over the unit line.
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
//    16 April 2009
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
//    Input, int EXPON[1], the exponent of X.  The exponent must not be -1.
//
//    Output, double LINE_UNIT_MONOMIAL, the integral of X^ALPHA.
//
{
  double value;

  if ( expon[0] == - 1 )
  {
    cerr << "\n";
    cerr << "LINE_UNIT_MONOMIAL\n";
    cerr << "  Exponent = -1 is not a legal input.\n";
    exit ( 1 );
  }
  else if ( ( expon[0] % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 / ( double ) ( expon[0] + 1 );
  }

  return value;
}
//****************************************************************************80

void line_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_MONOMIAL_TEST tests LINE_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int expon[1];
  double value;

  cout << "\n";
  cout << "LINE_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit line,\n";
  cout << "  LINE_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA\n";
  cout << "\n";
  cout << "  Volume = " << line_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    value = line_unit_monomial ( expon );

    cout << "  " << setw(8)  << expon[0]
         << "  " << setw(14) << value << "\n";
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

  for ( i = 0; i < order; i++ )
  {
    w[i] = w[i] / line_unit_volume ( );
  }

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

  for ( i = 0; i < order; i++ )
  {
    w[i] = w[i] / line_unit_volume ( );
  }

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

  for ( i = 0; i < order; i++ )
  {
    w[i] = w[i] / line_unit_volume ( );
  }

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

  for ( i = 0; i < order; i++ )
  {
    w[i] = w[i] / line_unit_volume ( );
  }

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

  for ( i = 0; i < order; i++ )
  {
    w[i] = w[i] / line_unit_volume ( );
  }

  return;
}
//****************************************************************************80

void line_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_QUAD_TEST tests the rules for the unit line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2008
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
# define DIM_NUM 1

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  bool more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *x;

  cout << "\n";
  cout << "LINE_UNIT_QUAD_TEST\n";
  cout << "  For the unit line,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  LINE_UNIT_O01, a 1 point rule.\n";
  cout << "  LINE_UNIT_O02, a 2 point rule.\n";
  cout << "  LINE_UNIT_O03, a 3 point rule.\n";
  cout << "  LINE_UNIT_O04, a 4 point rule.\n";
  cout << "  LINE_UNIT_O05, a 5 point rule.\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( i4vec_odd_any ( dim_num, expon ) )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    order = 1;
    w = new double[order];
    x = new double[dim_num*order];
    line_unit_o01 ( w, x );
    v = monomial_value ( dim_num, order, expon, x );
    quad = line_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    order = 2;
    w = new double[order];
    x = new double[dim_num*order];
    line_unit_o02 ( w, x );
    v = monomial_value ( dim_num, order, expon, x );
    quad = line_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    order = 3;
    w = new double[order];
    x = new double[dim_num*order];
    line_unit_o03 ( w, x );
    v = monomial_value ( dim_num, order, expon, x );
    quad = line_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    order = 4;
    w = new double[order];
    x = new double[dim_num*order];
    line_unit_o04 ( w, x );
    v = monomial_value ( dim_num, order, expon, x );
    quad = line_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    order = 5;
    w = new double[order];
    x = new double[dim_num*order];
    line_unit_o05 ( w, x );
    v = monomial_value ( dim_num, order, expon, x );
    quad = line_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] x;

    cout << "\n";
    quad = line_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double line_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    LINE_UNIT_VOLUME: volume of a unit line.
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
//    03 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double LINE_UNIT_VOLUME, the volume of the line.
//
{
  double volume;

  volume = 2.0;

  return volume;
}
//****************************************************************************80

double *monomial_value ( int dim_num, int point_num, int expon[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    MONOMIAL_VALUE evaluates a monomial.
//
//  Discussion:
//
//    F(X) = product ( 1 <= DIM <= DIM_NUM ) X(I)^EXPON(I)
//
//    with the convention that 0^0 = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    15 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int EXPON[DIM_NUM], the exponents.
//
//    Input, double X[DIM_NUM*POINT_NUM], the evaluation points.
//
//    Output, double MONOMIAL_VALUE[POINT_NUM], the monomial values.
//
{
  int dim;
  int point;
  double *v;

  v = new double[point_num];

  for ( point = 0; point < point_num; point++ )
  {
    v[point] = 1.0;
  }

  for ( dim = 0; dim < dim_num; dim++ )
  {
    if ( expon[dim] != 0.0 )
    {
      for ( point = 0; point < point_num; point++ )
      {
        v[point] = v[point] * pow ( x[dim+point*dim_num], expon[dim] );
      }
    }
  }

  return v;
}
//****************************************************************************80

double pyra_unit_monomial ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_MONOMIAL: monomial integral in a unit pyramid.
//
//  Discussion:
//
//    This function returns the value of the integral of X^ALPHA Y^BETA Z^GAMMA
//    over the unit pyramid.
//
//    The integration region is:
//
//    - ( 1 - Z ) <= X <= 1 - Z
//    - ( 1 - Z ) <= Y <= 1 - Z
//              0 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
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
//    Input, int EXPON[3], the exponents.
//
//    Output, double PYRA_UNIT_MONOMIAL, the volume of the pyramid.
//
{
  int i;
  int i_hi;
  double value;

  value = 0.0;

  if ( ( expon[0] % 2 ) == 0 && ( expon[1] % 2 ) == 0 )
  {
    i_hi = 2 + expon[0] + expon[1];

    for ( i = 0; i <= i_hi; i++ )
    {
      value = value + r8_mop ( i ) * r8_choose ( i_hi, i ) 
      / ( double ) ( i + expon[2] + 1 );
    }

    value = value 
          * 2.0 / ( double ) ( expon[0] + 1 )
          * 2.0 / ( double ) ( expon[1] + 1 );
  }

  return value;
}
//****************************************************************************80

void pyra_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_MONOMIAL_TEST tests PYRA_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[3];
  int gamma;
  double value;

  cout << "\n";
  cout << "PYRA_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit pyramid,\n";
  cout << "  PYRA_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << pyra_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = pyra_unit_monomial ( expon );

        cout << "  " << setw(8)  << expon[0]
             << "  " << setw(8)  << expon[1]
             << "  " << setw(8)  << expon[2]
             << "  " << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void pyra_unit_o01 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O01 returns a 1 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double XYZ[3*1], the abscissas.
//
{
  int order = 1;

  double w_save[1] = { 1.0 };

  double xyz_save[3*1] = { 0.00, 0.00, 0.25 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o05 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O05 returns a 5 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[5], the weights.
//
//    Output, double XYZ[3*5], the abscissas.
//
{
  int order = 5;

  double w_save[5] = {
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.21093750000000000000, 
   0.15625000000000000000 };

  double xyz_save[3*5] = {
  -0.48686449556014765641,   -0.48686449556014765641,   0.16666666666666666667, 
   0.48686449556014765641,   -0.48686449556014765641,   0.16666666666666666667,
   0.48686449556014765641,    0.48686449556014765641,   0.16666666666666666667, 
  -0.48686449556014765641,    0.48686449556014765641,   0.16666666666666666667,
   0.00000000000000000000,    0.00000000000000000000,   0.70000000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o06 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O06 returns a 6 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double XYZ[3*6], the abscissas.
//
{
  int order = 6;

  double w_save[6] = {
   0.21000000000000000000,
   0.21000000000000000000,
   0.21000000000000000000,
   0.21000000000000000000,
   0.06000000000000000000,
   0.10000000000000000000 };

  double xyz_save[3*6] = {
  -0.48795003647426658968,  -0.48795003647426658968,   0.16666666666666666667,
   0.48795003647426658968,  -0.48795003647426658968,   0.16666666666666666667,
   0.48795003647426658968,   0.48795003647426658968,   0.16666666666666666667,
  -0.48795003647426658968,   0.48795003647426658968,   0.16666666666666666667,
   0.00000000000000000000,   0.00000000000000000000,   0.58333333333333333333,
   0.00000000000000000000,   0.00000000000000000000,   0.75000000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o08 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O08 returns an 8 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[8], the weights.
//
//    Output, double XYZ[3*8], the abscissas.
//
{
  int order = 8;

 double w_save[8] = {
   0.075589411559869072938,
   0.075589411559869072938,
   0.075589411559869072938,
   0.075589411559869072938,
   0.17441058844013092706,
   0.17441058844013092706,
   0.17441058844013092706,
   0.17441058844013092706 };

  double xyz_save[3*8] = {
  -0.26318405556971359557,  -0.26318405556971359557,   0.54415184401122528880,
   0.26318405556971359557,  -0.26318405556971359557,   0.54415184401122528880,
   0.26318405556971359557,   0.26318405556971359557,   0.54415184401122528880,
  -0.26318405556971359557,   0.26318405556971359557,   0.54415184401122528880,
  -0.50661630334978742377,  -0.50661630334978742377,   0.12251482265544137787,
   0.50661630334978742377,  -0.50661630334978742377,   0.12251482265544137787,
   0.50661630334978742377,   0.50661630334978742377,   0.12251482265544137787,
  -0.50661630334978742377,   0.50661630334978742377,   0.12251482265544137787 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o08b ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O08B returns an 8 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[8], the weights.
//
//    Output, double XYZ[3*8], the abscissas.
//
{
  int order = 8;

  double w_save[8] = {
   0.16438287736328777572,
   0.16438287736328777572,
   0.16438287736328777572,
   0.16438287736328777572,
   0.085617122636712224276,
   0.085617122636712224276,
   0.085617122636712224276,
   0.085617122636712224276 };

  double xyz_save[3*8] = {
  -0.51197009372656270107,  -0.51197009372656270107,   0.11024490204163285720,
   0.51197009372656270107,  -0.51197009372656270107,   0.11024490204163285720,
   0.51197009372656270107,   0.51197009372656270107,   0.11024490204163285720,
  -0.51197009372656270107,   0.51197009372656270107,   0.11024490204163285720,
  -0.28415447557052037456,  -0.28415447557052037456,   0.518326526529795714229,
   0.28415447557052037456,  -0.28415447557052037456,   0.518326526529795714229,
   0.28415447557052037456,   0.28415447557052037456,   0.518326526529795714229,
  -0.28415447557052037456,   0.28415447557052037456,   0.518326526529795714229 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o09 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O09 returns a 9 point quadrature rule for the unit pyramid.
///
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[9], the weights.
//
//    Output, double XYZ[3*9], the abscissas.
//
{
  int order = 9;

  double w_save[9] = {
   0.13073389672275944791,
   0.13073389672275944791,
   0.13073389672275944791,
   0.13073389672275944791,
   0.10989110327724055209,
   0.10989110327724055209,
   0.10989110327724055209,
   0.10989110327724055209,
   0.03750000000000000000 };

  double xyz_save[3*9] = {
  -0.52966422253852215131,  -0.52966422253852215131,   0.08176876558246862335,
   0.52966422253852215131,  -0.52966422253852215131,   0.08176876558246862335,
   0.52966422253852215131,   0.52966422253852215131,   0.08176876558246862335,
  -0.52966422253852215131,   0.52966422253852215131,   0.08176876558246862335,
  -0.34819753825720418039,  -0.34819753825720418039,   0.400374091560388519511,
   0.34819753825720418039,  -0.34819753825720418039,   0.400374091560388519511,
   0.34819753825720418039,   0.34819753825720418039,   0.400374091560388519511,
  -0.34819753825720418039,   0.34819753825720418039,   0.400374091560388519511,
   0.00000000000000000000,   0.00000000000000000000,   0.83333333333333333333 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o13 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O13 returns a 13 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[13], the weights.
//
//    Output, double XYZ[3*13], the abscissas.
//
{
  int order = 13;

  double w_save[13] = {
   0.063061594202898550725,
   0.063061594202898550725,
   0.063061594202898550725,
   0.063061594202898550725,
   0.042101946815575556199,
   0.042101946815575556199,
   0.042101946815575556199,
   0.042101946815575556199,
   0.13172030707666776585,
   0.13172030707666776585,
   0.13172030707666776585,
   0.13172030707666776585,
   0.05246460761943250889 };

  double xyz_save[3*13] = {
  -0.38510399211870384331,  -0.38510399211870384331,  0.428571428571428571429,
   0.38510399211870384331,  -0.38510399211870384331,  0.428571428571428571429,
   0.38510399211870384331,   0.38510399211870384331,  0.428571428571428571429,
  -0.38510399211870384331,   0.38510399211870384331,  0.428571428571428571429,
  -0.40345831960728204766,   0.00000000000000000000,  0.33928571428571428571,
   0.40345831960728204766,   0.00000000000000000000,  0.33928571428571428571,
   0.00000000000000000000,  -0.40345831960728204766,  0.33928571428571428571,
   0.00000000000000000000,   0.40345831960728204766,  0.33928571428571428571,
  -0.53157877436961973359,  -0.53157877436961973359,  0.08496732026143790850,
   0.53157877436961973359,  -0.53157877436961973359,  0.08496732026143790850,
   0.53157877436961973359,   0.53157877436961973359,  0.08496732026143790850,
  -0.53157877436961973359,   0.53157877436961973359,  0.08496732026143790850,
   0.00000000000000000000,   0.00000000000000000000,  0.76219701803768503595 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o18 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O18 returns an 18 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[18], the weights.
//
//    Output, double XYZ[3*18], the abscissas.
//
{
  int order = 18;

  double w_save[18] = {
   0.023330065296255886709,
   0.037328104474009418735,
   0.023330065296255886709,
   0.037328104474009418735,
   0.059724967158415069975,
   0.037328104474009418735,
   0.023330065296255886709,
   0.037328104474009418735,
   0.023330065296255886709,
   0.05383042853090460712,
   0.08612868564944737139,
   0.05383042853090460712,
   0.08612868564944737139,
   0.13780589703911579422,
   0.08612868564944737139,
   0.05383042853090460712,
   0.08612868564944737139,
   0.05383042853090460712 };

  double xyz_save[3*18] = {
  -0.35309846330877704481,  -0.35309846330877704481,  0.544151844011225288800,
   0.00000000000000000000,  -0.35309846330877704481,  0.544151844011225288800,
   0.35309846330877704481,  -0.35309846330877704481,  0.544151844011225288800,
  -0.35309846330877704481,   0.00000000000000000000,  0.544151844011225288800,
   0.00000000000000000000,   0.00000000000000000000,  0.544151844011225288800,
   0.35309846330877704481,   0.00000000000000000000,  0.544151844011225288800,
  -0.35309846330877704481,   0.35309846330877704481,  0.544151844011225288800,
   0.00000000000000000000,   0.35309846330877704481,  0.544151844011225288800,
   0.35309846330877704481,   0.35309846330877704481,  0.544151844011225288800,
  -0.67969709567986745790,  -0.67969709567986745790,  0.12251482265544137787,
   0.00000000000000000000,  -0.67969709567986745790,  0.12251482265544137787,
   0.67969709567986745790,  -0.67969709567986745790,  0.12251482265544137787,
  -0.67969709567986745790,   0.00000000000000000000,  0.12251482265544137787,
   0.00000000000000000000,   0.00000000000000000000,  0.12251482265544137787,
   0.67969709567986745790,   0.00000000000000000000,  0.12251482265544137787,
  -0.67969709567986745790,   0.67969709567986745790,  0.12251482265544137787,
   0.00000000000000000000,   0.67969709567986745790,  0.12251482265544137787,
   0.67969709567986745790,   0.67969709567986745790,  0.12251482265544137787 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o27 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O27 returns a 27 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[27], the weights.
//
//    Output, double XYZ[3*27], the abscissas.
//
{
  int order = 27;

  double w_save[27] = {
   0.036374157653908938268,
   0.05819865224625430123,
   0.036374157653908938268,
   0.05819865224625430123,
   0.09311784359400688197,
   0.05819865224625430123,
   0.036374157653908938268,
   0.05819865224625430123,
   0.036374157653908938268,
   0.033853303069413431019,
   0.054165284911061489631,
   0.033853303069413431019,
   0.054165284911061489631,
   0.08666445585769838341,
   0.054165284911061489631,
   0.033853303069413431019,
   0.054165284911061489631,
   0.033853303069413431019,
   0.006933033103838124540,
   0.011092852966140999264,
   0.006933033103838124540,
   0.011092852966140999264,
   0.017748564745825598822,
   0.011092852966140999264,
   0.006933033103838124540,
   0.011092852966140999264,
   0.006933033103838124540 };

  double xyz_save[3*27] = {
  -0.7180557413198889387,   -0.7180557413198889387,   0.07299402407314973216,
   0.00000000000000000000,  -0.7180557413198889387,   0.07299402407314973216,
   0.7180557413198889387,   -0.7180557413198889387,   0.07299402407314973216,
  -0.7180557413198889387,    0.00000000000000000000,  0.07299402407314973216,
   0.00000000000000000000,   0.00000000000000000000,  0.07299402407314973216,
   0.7180557413198889387,    0.00000000000000000000,  0.07299402407314973216,
  -0.7180557413198889387,    0.7180557413198889387,   0.07299402407314973216,
   0.00000000000000000000,   0.7180557413198889387,   0.07299402407314973216,
   0.7180557413198889387,    0.7180557413198889387,   0.07299402407314973216,
  -0.50580870785392503961,  -0.50580870785392503961,  0.34700376603835188472,
   0.00000000000000000000,  -0.50580870785392503961,  0.34700376603835188472,
   0.50580870785392503961,  -0.50580870785392503961,  0.34700376603835188472,
  -0.50580870785392503961,   0.00000000000000000000,  0.34700376603835188472,
   0.00000000000000000000,   0.00000000000000000000,  0.34700376603835188472,
   0.50580870785392503961,   0.00000000000000000000,  0.34700376603835188472,
  -0.50580870785392503961,   0.50580870785392503961,  0.34700376603835188472,
   0.00000000000000000000,   0.50580870785392503961,  0.34700376603835188472,
   0.50580870785392503961,   0.50580870785392503961,  0.34700376603835188472,
  -0.22850430565396735360,  -0.22850430565396735360,  0.70500220988849838312,
   0.00000000000000000000,  -0.22850430565396735360,  0.70500220988849838312,
   0.22850430565396735360,  -0.22850430565396735360,  0.70500220988849838312,
  -0.22850430565396735360,   0.00000000000000000000,  0.70500220988849838312,
   0.00000000000000000000,   0.00000000000000000000,  0.70500220988849838312,
   0.22850430565396735360,   0.00000000000000000000,  0.70500220988849838312,
  -0.22850430565396735360,   0.22850430565396735360,  0.70500220988849838312,
   0.00000000000000000000,   0.22850430565396735360,  0.70500220988849838312,
   0.22850430565396735360,   0.22850430565396735360,  0.70500220988849838312 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_o48 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_O48 returns a 48 point quadrature rule for the unit pyramid.
//
//  Discussion:
//
//    The integration region is:
//
//     - ( 1 - Z ) <= X <= 1 - Z
//     - ( 1 - Z ) <= Y <= 1 - Z
//               0 <= Z <= 1.
//
//    When Z is zero, the integration region is a square lying in the (X,Y) 
//    plane, centered at (0,0,0) with "radius" 1.  As Z increases to 1, the 
//    radius of the square diminishes, and when Z reaches 1, the square has 
//    contracted to the single point (0,0,1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 April 2009
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
//    Output, double W[48], the weights.
//
//    Output, double XYZ[3*48], the abscissas.
//
{
  int order = 48;

  double w_save[48] = {
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.01241939442682455E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  2.60351137043010779E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.24557795239745531E-002, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  1.87873998794808156E-003, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  4.32957927807745280E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  1.97463249834127288E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  5.60127223523590526E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.55462562927473852E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  2.67977366291788643E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  1.22218992265373354E-002, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  4.04197740453215038E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003, 
  1.84346316995826843E-003 };

  double xyz_save[3*48] = {
  0.88091731624450909E+00,  0.00000000000000000E+00,  4.85005494469969989E-02,
 -0.88091731624450909E+00,  0.00000000000000000E+00,  4.85005494469969989E-02,
  0.00000000000000000E+00,  0.88091731624450909E+00,  4.85005494469969989E-02,
  0.00000000000000000E+00, -0.88091731624450909E+00,  4.85005494469969989E-02,
  0.70491874112648223E+00,  0.00000000000000000E+00,  0.23860073755186201E+00,
 -0.70491874112648223E+00,  0.00000000000000000E+00,  0.23860073755186201E+00,
  0.00000000000000000E+00,  0.70491874112648223E+00,  0.23860073755186201E+00,
  0.00000000000000000E+00, -0.70491874112648223E+00,  0.23860073755186201E+00,
  0.44712732143189760E+00,  0.00000000000000000E+00,  0.51704729510436798E+00,
 -0.44712732143189760E+00,  0.00000000000000000E+00,  0.51704729510436798E+00,
  0.00000000000000000E+00,  0.44712732143189760E+00,  0.51704729510436798E+00,
  0.00000000000000000E+00, -0.44712732143189760E+00,  0.51704729510436798E+00,
  0.18900486065123448E+00,  0.00000000000000000E+00,  0.79585141789677305E+00,
 -0.18900486065123448E+00,  0.00000000000000000E+00,  0.79585141789677305E+00,
  0.00000000000000000E+00,  0.18900486065123448E+00,  0.79585141789677305E+00,
  0.00000000000000000E+00, -0.18900486065123448E+00,  0.79585141789677305E+00,
  0.36209733410322176E+00,  0.36209733410322176E+00,  4.85005494469969989E-02,
 -0.36209733410322176E+00,  0.36209733410322176E+00,  4.85005494469969989E-02,
 -0.36209733410322176E+00, -0.36209733410322176E+00,  4.85005494469969989E-02,
  0.36209733410322176E+00, -0.36209733410322176E+00,  4.85005494469969989E-02,
  0.76688932060387538E+00,  0.76688932060387538E+00,  4.85005494469969989E-02,
 -0.76688932060387538E+00,  0.76688932060387538E+00,  4.85005494469969989E-02,
 -0.76688932060387538E+00, -0.76688932060387538E+00,  4.85005494469969989E-02,
  0.76688932060387538E+00, -0.76688932060387538E+00,  4.85005494469969989E-02,
  0.28975386476618070E+00,  0.28975386476618070E+00,  0.23860073755186201E+00,
 -0.28975386476618070E+00,  0.28975386476618070E+00,  0.23860073755186201E+00,
 -0.28975386476618070E+00, -0.28975386476618070E+00,  0.23860073755186201E+00,
  0.28975386476618070E+00, -0.28975386476618070E+00,  0.23860073755186201E+00,
  0.61367241226233160E+00,  0.61367241226233160E+00,  0.23860073755186201E+00,
 -0.61367241226233160E+00,  0.61367241226233160E+00,  0.23860073755186201E+00,
 -0.61367241226233160E+00, -0.61367241226233160E+00,  0.23860073755186201E+00,
  0.61367241226233160E+00, -0.61367241226233160E+00,  0.23860073755186201E+00,
  0.18378979287798017E+00,  0.18378979287798017E+00,  0.51704729510436798E+00,
 -0.18378979287798017E+00,  0.18378979287798017E+00,  0.51704729510436798E+00,
 -0.18378979287798017E+00, -0.18378979287798017E+00,  0.51704729510436798E+00,
  0.18378979287798017E+00, -0.18378979287798017E+00,  0.51704729510436798E+00,
  0.38925011625173161E+00,  0.38925011625173161E+00,  0.51704729510436798E+00,
 -0.38925011625173161E+00,  0.38925011625173161E+00,  0.51704729510436798E+00,
 -0.38925011625173161E+00, -0.38925011625173161E+00,  0.51704729510436798E+00,
  0.38925011625173161E+00, -0.38925011625173161E+00,  0.51704729510436798E+00,
  7.76896479525748113E-02,  7.76896479525748113E-02,  0.79585141789677305E+00,
 -7.76896479525748113E-02,  7.76896479525748113E-02,  0.79585141789677305E+00,
 -7.76896479525748113E-02, -7.76896479525748113E-02,  0.79585141789677305E+00,
  7.76896479525748113E-02, -7.76896479525748113E-02,  0.79585141789677305E+00,
  0.16453962988669860E+00,  0.16453962988669860E+00,  0.79585141789677305E+00,
 -0.16453962988669860E+00,  0.16453962988669860E+00,  0.79585141789677305E+00,
 -0.16453962988669860E+00, -0.16453962988669860E+00,  0.79585141789677305E+00,
  0.16453962988669860E+00, -0.16453962988669860E+00,  0.79585141789677305E+00 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void pyra_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_QUAD_TEST tests the rules for the unit pyramid.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2008
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
# define DIM_NUM 3

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  bool more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "PYRA_UNIT_QUAD_TEST\n";
  cout << "  For the unit pyramid,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  PYRA_UNIT_O01,\n";
  cout << "  PYRA_UNIT_O05,\n";
  cout << "  PYRA_UNIT_O06,\n";
  cout << "  PYRA_UNIT_O08,\n";
  cout << "  PYRA_UNIT_O08b,\n";
  cout << "  PYRA_UNIT_O09,\n";
  cout << "  PYRA_UNIT_O13,\n";
  cout << "  PYRA_UNIT_O18,\n";
  cout << "  PYRA_UNIT_O27,\n";
  cout << "  PYRA_UNIT_O48,\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[0] % 2 ) == 1 || ( expon[1] % 2 ) == 1 )
    {
      continue;
    }

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    order = 1;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o01 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 5;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o05 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 6;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o06 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o08 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o08b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 9;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o09 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 13;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o13 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 18;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o18 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 27;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o27 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 48;
    w = new double[order];
    xyz = new double[dim_num*order];
    pyra_unit_o48 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = pyra_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    cout << "\n";
    quad = pyra_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double pyra_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    PYRA_UNIT_VOLUME: volume of a unit pyramid with square base.
//
//  Discussion:
//
//    The volume of this unit pyramid is 4/3.
//
//    The integration region is:
//
//      - ( 1 - Z ) <= X <= 1 - Z
//      - ( 1 - Z ) <= Y <= 1 - Z
//                0 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double PYRA_UNIT_VOLUME, the volume of the pyramid.
//
{
  double volume;

  volume = 4.0 / 3.0;

  return volume;
}
//****************************************************************************80

double quad_unit_monomial ( int expon[2] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_UNIT_MONOMIAL integrates a monomial over the unit quadrilateral.
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
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON[2], the exponents.
//
//    Output, double QUAD_UNIT_MONOMIAL, the integral of the monomial.
//
{
  int i;
  double value;

  value = 1.0;

  for ( i = 0; i < 2; i++ )
  {
    if ( ( expon[i] % 2 ) == 1 )
    {
      value = 0.0;
    }
    else if ( expon[i] == -1 )
    {
      cout << "\n";
      cout << "QUAD_UNIT_MONOMIAL - Fatal error!\n";
      cout << "  Exponent of -1 encountered.\n";
      exit ( 1 );
    }
    else
    {
      value = value * 2.0 / ( double ) ( expon[i] + 1 );
    }
  }

  return value;
}
//****************************************************************************80

void quad_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_UNIT_MONOMIAL_TEST tests QUAD_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[2];
  double value;

  cout << "\n";
  cout << "QUAD_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit quadrilateral,\n";
  cout << "  QUAD_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA\n";
  cout << "\n";
  cout << "  Volume = " << quad_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;

      value = quad_unit_monomial ( expon );

      cout << "  " << setw(8)  << expon[0]
           << "  " << setw(8)  << expon[1]
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void quad_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_UNIT_QUAD_TEST tests the rules for the unit quadrilateral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 April 2008
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
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  int k;
  bool more;
  int order;
  int order_1d[DIM_NUM];
  double quad;
  int t;
  double *v;
  double *w;
  double *xy;

  cout << "\n";
  cout << "QUAD_UNIT_QUAD_TEST\n";
  cout << "  For the unit quadrilateral,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  QUAD_UNIT_RULE, which returns M by N point rules..\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( i4vec_odd_any ( dim_num, expon ) )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    for ( k = 1; k <= 5; k++ )
    {
      for ( dim = 0; dim < dim_num; dim++ )
      {
        order_1d[dim] = k;
      }
      order = i4vec_product ( dim_num, order_1d );
      w = new double[order];
      xy = new double[dim_num*order];
      quad_unit_rule ( order_1d, w, xy );
      v = monomial_value ( dim_num, order, expon, xy );
      quad = quad_unit_volume ( ) * r8vec_dot_product ( order, w, v );
      cout << "  " << setw(6) << order_1d[0]
           << "  " << setw(6) << order_1d[1]
           << "  " << setw(14) << quad << "\n";
      delete [] v;
      delete [] w;
      delete [] xy;
    }
//
//  Try a rule of mixed orders.
//
    order_1d[0] = 3;
    order_1d[1] = 5;
    order = i4vec_product ( dim_num, order_1d );
    w = new double[order];
    xy = new double[dim_num*order];
    quad_unit_rule ( order_1d, w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = quad_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order_1d[0]
         << "  " << setw(6) << order_1d[1]
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    cout << "\n";
    quad = quad_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << "      "
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void quad_unit_rule ( int order_1d[], double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_UNIT_RULE returns a quadrature rule for the unit quadrilateral.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
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
//  Reference:
//
//    Carlos Felippa,
//    A compendium of FEM integration formulas for symbolic work,
//    Engineering Computation,
//    Volume 21, Number 8, 2004, pages 867-890.
//
//  Parameters:
//
//    Input, int ORDER_1D[2], the order of the rule in 
//    each dimension.  1 <= ORDER_1D(I) <= 5.
//
//    Output, double W[ORDER_1D[0]*ORDER_1D[1]], the weights.
//
//    Output, double XY[2*ORDER_1D[0]*ORDER_1D[1]], the abscissas.
//
{
  int dim;
  int dim_num = 2;
  int order;
  double *w_1d;
  double *x_1d;

  order = i4vec_product ( dim_num, order_1d );

  for ( dim = 0; dim < dim_num; dim++ )
  {
    w_1d = new double[order_1d[dim]];
    x_1d = new double[order_1d[dim]];

    if ( order_1d[dim] == 1 )
    {
      line_unit_o01 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 2 )
    {
      line_unit_o02 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 3 )
    {
      line_unit_o03 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 4 )
    {
      line_unit_o04 ( w_1d, x_1d );
    }
    else if ( order_1d[dim] == 5 )
    {
      line_unit_o05 ( w_1d, x_1d );
    }
    else
    {
      cout << "\n";
      cout << "QUAD_UNIT_RULE - Fatal error!\n";
      cout << "  Illegal value of ORDER_1D[*].\n";
      exit ( 1 );
    }

    r8vec_direct_product ( dim, order_1d[dim], x_1d, 
      dim_num, order, xy );

    r8vec_direct_product2 ( dim, order_1d[dim], w_1d, 
      dim_num, order, w );

    delete [] w_1d;
    delete [] x_1d;
  }

  return;
}
//****************************************************************************80

double quad_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    QUAD_UNIT_VOLUME: volume of a unit quadrilateral.
//
//  Discussion:
//
//    The integration region is:
//
//    - 1.0 <= X <= 1.0
//    - 1.0 <= Y <= 1.0
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
//  Parameters:
//
//    Output, double QUAD_UNIT_VOLUME, the volume.
//
{
  double value;

  value = 4.0;

  return value;
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

double r8_choose ( int n, int k )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
{
  int i;
  int mn;
  int mx;
  int value;

  mn = i4_min ( k, n - k );

  if ( mn < 0 )
  {
    value = 0.0;
  }
  else if ( mn == 0 )
  {
    value = 1.0;
  }
  else
  {
    mx = i4_max ( k, n - k );
    value = ( double ) ( mx + 1 );

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( double ) ( mx + i ) ) / ( double ) i;
    }
  }

  return value;
}
//****************************************************************************80

double r8_mop ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
{
  double value;

  if ( ( i % 2 ) == 0 )
  {
    value = 1.0;
  }
  else
  {
    value = -1.0;
  }

  return value;
}
//****************************************************************************80

double r8mat_det_4d ( double a[4*4] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_DET_4D computes the determinant of a 4 by 4 R8MAT.
//
//  Discussion:
//
//    The two dimensional array is stored as a one dimensional vector,
//    by COLUMNS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    10 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[4*4], the matrix whose determinant is desired.
//
//    Output, double R8MAT_DET_4D, the determinant of the matrix.
//
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}
//****************************************************************************80

void r8mat_write ( string output_filename, int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_WRITE writes an R8MAT file with no header.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double TABLE[M*N], the table data.
//
{
  int i;
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8MAT_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }
//
//  Write the data.
//  For greater precision, try
//
//    output << "  " << setw(24) << setprecision(16) << table[i+j*m];
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << "  " << setw(24) << setprecision(16) << table[i+j*m];
    }
    output << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
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

void r8vec_direct_product ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2, 
//      ..., 
//      item JK of 1D rule K.
//
//    In particular, 
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1: 
//      Order = 4
//      X(1:4) = ( 1, 2, 3, 4 )
//
//    Rule 2:
//      Order = 3
//      X(1:3) = ( 10, 20, 30 )
//
//    Rule 3:
//      Order = 2
//      X(1:2) = ( 100, 200 )
//
//    Product Rule:
//      Order = 24
//      X(1:24) = 
//        ( 1, 10, 100 )
//        ( 2, 10, 100 )
//        ( 3, 10, 100 )
//        ( 4, 10, 100 )
//        ( 1, 20, 100 )
//        ( 2, 20, 100 )
//        ( 3, 20, 100 )
//        ( 4, 20, 100 )
//        ( 1, 30, 100 )
//        ( 2, 30, 100 )
//        ( 3, 30, 100 )
//        ( 4, 30, 100 )
//        ( 1, 10, 200 )
//        ( 2, 10, 200 )
//        ( 3, 10, 200 )
//        ( 4, 10, 200 )
//        ( 1, 20, 200 )
//        ( 2, 20, 200 )
//        ( 3, 20, 200 )
//        ( 4, 20, 200 )
//        ( 1, 30, 200 )
//        ( 2, 30, 200 )
//        ( 3, 30, 200 )
//        ( 4, 30, 200 )
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
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values
//    for factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double X[FACTOR_NUM*POINT_NUM], the elements of the
//    direct product, which are built up gradually.  
//
//  Local Parameters:
//
//    Local, int START, the first location of a block of values to set.
//
//    Local, int CONTIG, the number of consecutive values to set.
//
//    Local, int SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, int REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < factor_num; i++ )
    {
      for ( j = 0; j < point_num; j++ )
      {
        x[i+j*factor_num] = 0.0;
      }
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        x[factor_index+i*factor_num] = factor_value[j];
      }
      start = start + skip;
    }
  }
  contig = contig * factor_order;

  return;
}
//****************************************************************************80

void r8vec_direct_product2 ( int factor_index, int factor_order, 
  double factor_value[], int factor_num, int point_num, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    To explain what is going on here, suppose we had to construct
//    a multidimensional quadrature rule as the product of K rules
//    for 1D quadrature.
//
//    The product rule will be represented as a list of points and weights.
//
//    The J-th item in the product rule will be associated with
//      item J1 of 1D rule 1,
//      item J2 of 1D rule 2, 
//      ..., 
//      item JK of 1D rule K.
//
//    In particular, 
//      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
//    and
//      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
//
//    So we can construct the quadrature rule if we can properly
//    distribute the information in the 1D quadrature rules.
//
//    This routine carries out that task for the weights W.
//
//    Another way to do this would be to compute, one by one, the
//    set of all possible indices (J1,J2,...,JK), and then index
//    the appropriate information.  An advantage of the method shown
//    here is that you can process the K-th set of information and
//    then discard it.
//
//  Example:
//
//    Rule 1: 
//      Order = 4
//      W(1:4) = ( 2, 3, 5, 7 )
//
//    Rule 2:
//      Order = 3
//      W(1:3) = ( 11, 13, 17 )
//
//    Rule 3:
//      Order = 2
//      W(1:2) = ( 19, 23 )
//
//    Product Rule:
//      Order = 24
//      W(1:24) =
//        ( 2 * 11 * 19 )
//        ( 3 * 11 * 19 )
//        ( 4 * 11 * 19 )
//        ( 7 * 11 * 19 )
//        ( 2 * 13 * 19 )
//        ( 3 * 13 * 19 )
//        ( 5 * 13 * 19 )
//        ( 7 * 13 * 19 )
//        ( 2 * 17 * 19 )
//        ( 3 * 17 * 19 )
//        ( 5 * 17 * 19 )
//        ( 7 * 17 * 19 )
//        ( 2 * 11 * 23 )
//        ( 3 * 11 * 23 )
//        ( 5 * 11 * 23 )
//        ( 7 * 11 * 23 )
//        ( 2 * 13 * 23 )
//        ( 3 * 13 * 23 )
//        ( 5 * 13 * 23 )
//        ( 7 * 13 * 23 )
//        ( 2 * 17 * 23 )
//        ( 3 * 17 * 23 )
//        ( 5 * 17 * 23 )
//        ( 7 * 17 * 23 )
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
//  Parameters:
//
//    Input, int FACTOR_INDEX, the index of the factor being processed.
//    The first factor processed must be factor 0.
//
//    Input, int FACTOR_ORDER, the order of the factor.
//
//    Input, double FACTOR_VALUE[FACTOR_ORDER], the factor values for
//    factor FACTOR_INDEX.
//
//    Input, int FACTOR_NUM, the number of factors.
//
//    Input, int POINT_NUM, the number of elements in the direct product.
//
//    Input/output, double W[POINT_NUM], the elements of the
//    direct product, which are built up gradually.  
//
//  Local Parameters:
//
//    Local, integer START, the first location of a block of values to set.
//
//    Local, integer CONTIG, the number of consecutive values to set.
//
//    Local, integer SKIP, the distance from the current value of START
//    to the next location of a block of values to set.
//
//    Local, integer REP, the number of blocks of values to set.
//
{
  static int contig = 0;
  int i;
  int j;
  int k;
  static int rep = 0;
  static int skip = 0;
  int start;

  if ( factor_index == 0 )
  {
    contig = 1;
    skip = 1;
    rep = point_num;
    for ( i = 0; i < point_num; i++ )
    {
      w[i] = 1.0;
    }
  }

  rep = rep / factor_order;
  skip = skip * factor_order;

  for ( j = 0; j < factor_order; j++ )
  {
    start = 0 + j * contig;

    for ( k = 1; k <= rep; k++ )
    {
      for ( i = start; i < start + contig; i++ )
      {
        w[i] = w[i] * factor_value[j];
      }
      start = start + skip;
    }
  }

  contig = contig * factor_order;

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

double tetr_unit_monomial ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_MONOMIAL integrates a monomial over the unit tetrahedron.
//
//  Discussion:
//
//    This routine integrates a monomial of the form
//
//      product ( 1 <= dim <= 3 ) x(dim)^expon(dim)
//
//    where the exponents are nonnegative integers.  Note that
//    if the combination 0^0 is encountered, it should be treated
//    as 1.
//
//    Integral ( over unit tetrahedron ) x^l y^m z^n dx dy = 
//    l! * m! * n! / ( m + n + 3 )!
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      0 <= X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int EXPON[3], the exponents.
//
//    Output, double TETR_UNIT_MONOMIAL, the integral of the monomial.
//
{
  int i;
  int k;
  double value;
//
//  The first computation ends with VALUE = 1.0;
//
  value = 1.0;
//
//  The first loop simply calculates 1, so we short circuit it.
//
// k = 0;
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

  for ( i = 1; i <= expon[2]; i++ )
  {
    k = k + 1;
    value = value * ( double ) ( i ) / ( double ) ( k );
  }

  k = k + 1;
  value = value / ( double ) ( k );

  k = k + 1;
  value = value / ( double ) ( k );

  k = k + 1;
  value = value / ( double ) ( k ) ;

  return value;
}
//****************************************************************************80

void tetr_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_MONOMIAL_TEST tests TETR_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[3];
  int gamma;
  double value;

  cout << "\n";
  cout << "TETR_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit tetrahedron,\n";
  cout << "  TETR_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << tetr_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = tetr_unit_monomial ( expon );

        cout << "  " << setw(8)  << expon[0]
             << "  " << setw(8)  << expon[1]
             << "  " << setw(8)  << expon[2]
             << "  " << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void tetr_unit_o01 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O01 returns a 1 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[1], the weights.
//
//    Output, double XYZ[3*1], the abscissas.
//
{
  int order = 1;

  double w_save[1] = {
    1.0000000000000000000 };
  double xyz_save[3*1] = { 
    0.25000000000000000000,  0.25000000000000000000,  0.25000000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o04 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O04 returns a 4 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[4], the weights.
//
//    Output, double XYZ[3*4], the abscissas.
//
{
  int order = 4;

  double w_save[4] = { 
    0.25000000000000000000, 
    0.25000000000000000000, 
    0.25000000000000000000, 
    0.25000000000000000000 };
  double xyz_save[3*4] = { 
    0.58541019662496845446,  0.13819660112501051518,  0.13819660112501051518, 
    0.13819660112501051518,  0.58541019662496845446,  0.13819660112501051518, 
    0.13819660112501051518,  0.13819660112501051518,  0.58541019662496845446, 
    0.13819660112501051518,  0.13819660112501051518,  0.13819660112501051518 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o08 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O08 returns an 8 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[8], the weights.
//
//    Output, double XYZ[3*8], the abscissas.
//
{
  int order = 8;

  double w_save[8] = {
    0.13852796651186214232, 
    0.13852796651186214232, 
    0.13852796651186214232, 
    0.13852796651186214232, 
    0.11147203348813785768, 
    0.11147203348813785768, 
    0.11147203348813785768, 
    0.11147203348813785768 };
  double xyz_save[3*8] = { 
    0.015835909865720057993,  0.32805469671142664734,  0.32805469671142664734, 
    0.32805469671142664734,  0.015835909865720057993,  0.32805469671142664734, 
    0.32805469671142664734,  0.32805469671142664734,  0.015835909865720057993, 
    0.32805469671142664734,  0.32805469671142664734,  0.32805469671142664734, 
    0.67914317820120795168,  0.10695227393293068277,  0.10695227393293068277, 
    0.10695227393293068277,  0.67914317820120795168,  0.10695227393293068277, 
    0.10695227393293068277,  0.10695227393293068277,  0.67914317820120795168, 
    0.10695227393293068277,  0.10695227393293068277,  0.10695227393293068277 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o08b ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O08B returns an 8 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[8], the weights.
//
//    Output, double XYZ[3*8], the abscissas.
//
{
  int order = 8;

  double w_save[8] = {
    0.025000000000000000000, 
    0.025000000000000000000, 
    0.025000000000000000000, 
    0.025000000000000000000, 
    0.22500000000000000000, 
    0.22500000000000000000, 
    0.22500000000000000000, 
    0.22500000000000000000 };
  double xyz_save[3*8] = { 
    1.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 
    0.00000000000000000000,  1.00000000000000000000,  0.00000000000000000000, 
    0.00000000000000000000,  0.00000000000000000000,  1.00000000000000000000, 
    0.00000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 
    0.00000000000000000000,  0.33333333333333333333,  0.33333333333333333333, 
    0.33333333333333333333,  0.00000000000000000000,  0.33333333333333333333, 
    0.33333333333333333333,  0.33333333333333333333,  0.00000000000000000000, 
    0.33333333333333333333,  0.33333333333333333333,  0.33333333333333333333 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o14 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O14 returns a 14 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[ORDER], the weights.
//
//    Output, double XYZ[3*ORDER], the abscissas.
//
{
  int order = 14;

  double w_save[14] = {
    0.073493043116361949544, 
    0.073493043116361949544, 
    0.073493043116361949544, 
    0.073493043116361949544, 
    0.11268792571801585080, 
    0.11268792571801585080, 
    0.11268792571801585080, 
    0.11268792571801585080, 
    0.042546020777081466438, 
    0.042546020777081466438, 
    0.042546020777081466438, 
    0.042546020777081466438, 
    0.042546020777081466438, 
    0.042546020777081466438 };
  double xyz_save[3*14] = { 
    0.72179424906732632079,  0.092735250310891226402,  0.092735250310891226402, 
    0.092735250310891226402,  0.72179424906732632079,  0.092735250310891226402, 
    0.092735250310891226402,  0.092735250310891226402,  0.72179424906732632079, 
    0.092735250310891226402,  0.092735250310891226402,  0.092735250310891226402, 
    0.067342242210098170608,  0.31088591926330060980,  0.31088591926330060980, 
    0.31088591926330060980,  0.067342242210098170608,  0.31088591926330060980, 
    0.31088591926330060980,  0.31088591926330060980,  0.067342242210098170608, 
    0.31088591926330060980,  0.31088591926330060980,  0.31088591926330060980, 
    0.045503704125649649492,  0.045503704125649649492,  0.45449629587435035051, 
    0.045503704125649649492,  0.45449629587435035051,  0.045503704125649649492, 
    0.045503704125649649492,  0.45449629587435035051,  0.45449629587435035051, 
    0.45449629587435035051,  0.045503704125649649492,  0.045503704125649649492, 
    0.45449629587435035051,  0.045503704125649649492,  0.45449629587435035051, 
    0.45449629587435035051,  0.45449629587435035051,  0.045503704125649649492 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o14b ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O14B returns a 14 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[14], the weights.
//
//    Output, double XYZ[3*14], the abscissas.
//
{
  int order = 14;

  double w_save[14] = {
    0.13283874668559071814, 
    0.13283874668559071814, 
    0.13283874668559071814, 
    0.13283874668559071814, 
    0.088589824742980710434, 
    0.088589824742980710434, 
    0.088589824742980710434, 
    0.088589824742980710434, 
    0.019047619047619047619, 
    0.019047619047619047619, 
    0.019047619047619047619, 
    0.019047619047619047619, 
    0.019047619047619047619, 
    0.019047619047619047619  };
  double xyz_save[3*14] = { 
    0.056881379520423421748,  0.31437287349319219275,  0.31437287349319219275, 
    0.31437287349319219275,  0.056881379520423421748,  0.31437287349319219275, 
    0.31437287349319219275,  0.31437287349319219275,  0.056881379520423421748, 
    0.31437287349319219275,  0.31437287349319219275,  0.31437287349319219275, 
    0.69841970432438656092,  0.10052676522520447969,  0.10052676522520447969, 
    0.10052676522520447969,  0.69841970432438656092,  0.10052676522520447969, 
    0.10052676522520447969,  0.10052676522520447969,  0.69841970432438656092, 
    0.10052676522520447969,  0.10052676522520447969,  0.10052676522520447969, 
    0.50000000000000000000,  0.50000000000000000000,  0.00000000000000000000, 
    0.50000000000000000000,  0.00000000000000000000,  0.50000000000000000000, 
    0.50000000000000000000,  0.00000000000000000000,  0.00000000000000000000, 
    0.00000000000000000000,  0.50000000000000000000,  0.50000000000000000000, 
    0.00000000000000000000,  0.50000000000000000000,  0.00000000000000000000, 
    0.00000000000000000000,  0.00000000000000000000,  0.50000000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o15 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O15 returns a 15 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[15], the weights.
//
//    Output, double XYZ[3*15], the abscissas.
//
{
  int order = 15;

  double w_save[15] = {
    0.071937083779018620010, 
    0.071937083779018620010, 
    0.071937083779018620010, 
    0.071937083779018620010, 
    0.069068207226272385281, 
    0.069068207226272385281, 
    0.069068207226272385281, 
    0.069068207226272385281, 
    0.052910052910052910053, 
    0.052910052910052910053, 
    0.052910052910052910053, 
    0.052910052910052910053, 
    0.052910052910052910053, 
    0.052910052910052910053, 
    0.11851851851851851852 };
  double xyz_save[3*15] = { 
    0.72408676584183090163,  0.091971078052723032789,  0.091971078052723032789, 
    0.091971078052723032789,  0.72408676584183090163,  0.091971078052723032789, 
    0.091971078052723032789,  0.091971078052723032789,  0.72408676584183090163, 
    0.091971078052723032789,  0.091971078052723032789,  0.091971078052723032789, 
    0.040619116511110274837,  0.31979362782962990839,  0.31979362782962990839, 
    0.31979362782962990839,  0.040619116511110274837,  0.31979362782962990839, 
    0.31979362782962990839,  0.31979362782962990839,  0.040619116511110274837, 
    0.31979362782962990839,  0.31979362782962990839,  0.31979362782962990839, 
    0.44364916731037084426,  0.44364916731037084426,  0.056350832689629155741, 
    0.44364916731037084426,  0.056350832689629155741,  0.44364916731037084426, 
    0.44364916731037084426,  0.056350832689629155741,  0.056350832689629155741, 
    0.056350832689629155741,  0.44364916731037084426,  0.44364916731037084426, 
    0.056350832689629155741,  0.44364916731037084426,  0.056350832689629155741, 
    0.056350832689629155741,  0.056350832689629155741,  0.44364916731037084426, 
    0.25000000000000000000,  0.25000000000000000000,  0.25000000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o15b ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O15B returns a 15 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[15], the weights.
//
//    Output, double XYZ[3*15], the abscissas.
//
{
  int order = 15;

  double w_save[15] = {
    0.036160714285714285714, 
    0.036160714285714285714, 
    0.036160714285714285714, 
    0.036160714285714285714, 
    0.069871494516173816465, 
    0.069871494516173816465, 
    0.069871494516173816465, 
    0.069871494516173816465, 
    0.065694849368318756074, 
    0.065694849368318756074, 
    0.065694849368318756074, 
    0.065694849368318756074, 
    0.065694849368318756074, 
    0.065694849368318756074, 
    0.18170206858253505484 };
  double xyz_save[3*15] = { 
    0.00000000000000000000,  0.33333333333333333333,  0.33333333333333333333, 
    0.33333333333333333333,  0.00000000000000000000,  0.33333333333333333333, 
    0.33333333333333333333,  0.33333333333333333333,  0.00000000000000000000, 
    0.33333333333333333333,  0.33333333333333333333,  0.33333333333333333333, 
    0.72727272727272727273,  0.090909090909090909091,  0.090909090909090909091, 
    0.090909090909090909091,  0.72727272727272727273,  0.090909090909090909091, 
    0.090909090909090909091,  0.090909090909090909091,  0.72727272727272727273, 
    0.090909090909090909091,  0.090909090909090909091,  0.090909090909090909091, 
    0.43344984642633570176,  0.43344984642633570176,  0.066550153573664298240, 
    0.43344984642633570176,  0.066550153573664298240,  0.43344984642633570176, 
    0.43344984642633570176,  0.066550153573664298240,  0.066550153573664298240, 
    0.066550153573664298240,  0.43344984642633570176,  0.43344984642633570176, 
    0.066550153573664298240,  0.43344984642633570176,  0.066550153573664298240, 
    0.066550153573664298240,  0.066550153573664298240,  0.43344984642633570176, 
    0.25000000000000000000,  0.25000000000000000000,  0.250000000000000000 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_o24 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_O24 returns a 24 point quadrature rule for the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      0 <= Z
//      X + Y + Z <= 1.
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
//    Output, double W[24], the weights.
//
//    Output, double XYZ[3*24], the abscissas.
//
{
  int order = 24;

  double w_save[24] = {
    0.039922750257869636194, 
    0.039922750257869636194, 
    0.039922750257869636194, 
    0.039922750257869636194, 
    0.010077211055345822612, 
    0.010077211055345822612, 
    0.010077211055345822612, 
    0.010077211055345822612, 
    0.055357181543927398338, 
    0.055357181543927398338, 
    0.055357181543927398338, 
    0.055357181543927398338, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286, 
    0.048214285714285714286 };
  double xyz_save[3*24] = { 
    0.35619138622025439121,  0.21460287125991520293,  0.21460287125991520293, 
    0.21460287125991520293,  0.35619138622025439121,  0.21460287125991520293, 
    0.21460287125991520293,  0.21460287125991520293,  0.35619138622025439121, 
    0.21460287125991520293,  0.21460287125991520293,  0.21460287125991520293, 
    0.87797812439616594065,  0.040673958534611353116,  0.040673958534611353116, 
    0.040673958534611353116,  0.87797812439616594065,  0.040673958534611353116, 
    0.040673958534611353116,  0.040673958534611353116,  0.87797812439616594065, 
    0.040673958534611353116,  0.040673958534611353116,  0.040673958534611353116, 
    0.032986329573173468968,  0.32233789014227551034,  0.32233789014227551034, 
    0.32233789014227551034,  0.032986329573173468968,  0.32233789014227551034, 
    0.32233789014227551034,  0.32233789014227551034,  0.032986329573173468968, 
    0.32233789014227551034,  0.32233789014227551034,  0.32233789014227551034, 
    0.60300566479164914137,  0.26967233145831580803,  0.063661001875017525299, 
    0.60300566479164914137,  0.063661001875017525299,  0.26967233145831580803, 
    0.60300566479164914137,  0.063661001875017525299,  0.063661001875017525299, 
    0.063661001875017525299,  0.60300566479164914137,  0.26967233145831580803, 
    0.063661001875017525299,  0.60300566479164914137,  0.063661001875017525299, 
    0.063661001875017525299,  0.063661001875017525299,  0.60300566479164914137, 
    0.26967233145831580803,  0.60300566479164914137,  0.063661001875017525299, 
    0.26967233145831580803,  0.063661001875017525299,  0.60300566479164914137, 
    0.26967233145831580803,  0.063661001875017525299,  0.063661001875017525299, 
    0.063661001875017525299,  0.26967233145831580803,  0.60300566479164914137, 
    0.063661001875017525299,  0.26967233145831580803,  0.063661001875017525299, 
    0.063661001875017525299,  0.063661001875017525299,  0.26967233145831580803 };

  r8vec_copy ( order, w_save, w );
  r8vec_copy ( 3 * order, xyz_save, xyz );

  return;
}
//****************************************************************************80

void tetr_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_QUAD_TEST tests the rules for the unit tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 April 2008
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
# define DIM_NUM 3

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  bool more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "TETR_UNIT_QUAD_TEST\n";
  cout << "  For the unit tetrahedron,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  TETR_UNIT_O01,\n";
  cout << "  TETR_UNIT_O04,\n";
  cout << "  TETR_UNIT_O08,\n";
  cout << "  TETR_UNIT_O08b,\n";
  cout << "  TETR_UNIT_O14,\n";
  cout << "  TETR_UNIT_O14b,\n";
  cout << "  TETR_UNIT_O15,\n";
  cout << "  TETR_UNIT_O15b,\n";
  cout << "  TETR_UNIT_O24,\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    order = 1;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o01 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 4;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o04 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o08 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 8;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o08b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 14;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o14 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 14;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o14b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 15;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o15 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 15;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o15b ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    order = 24;
    w = new double[order];
    xyz = new double[dim_num*order];
    tetr_unit_o24 ( w, xyz );
    v = monomial_value ( dim_num, order, expon, xyz );
    quad = tetr_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xyz;

    cout << "\n";
    quad = tetr_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double tetr_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TETR_UNIT_VOLUME returns the volume of the unit tetrahedron.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X,
//      0 <= Y,
//      0 <= Z, 
//      X + Y + Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double TETR_UNIT_VOLUME, the volume.
//
{
  double volume;

  volume = 1.0 / 6.0;

  return volume;
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

double trig_unit_monomial ( int expon[2] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_MONOMIAL integrates a monomial over the unit triangle.
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
//    Output, double TRIG_UNIT_MONOMIAL, the integral of the monomial.
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

void trig_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_MONOMIAL_TEST tests TRIG_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[2];
  double value;

  cout << "\n";
  cout << "TRIG_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit triangle,\n";
  cout << "  TRIG_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA\n";
  cout << "\n";
  cout << "  Volume = " << trig_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;

      value = trig_unit_monomial ( expon );

      cout << "  " << setw(8)  << expon[0]
           << "  " << setw(8)  << expon[1]
           << "  " << setw(14) << value << "\n";
    }
  }

  return;
}
//****************************************************************************80

void trig_unit_o01 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O01 returns a 1 point quadrature rule for the unit triangle.
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

void trig_unit_o03 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O03 returns a 3 point quadrature rule for the unit triangle.
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

void trig_unit_o03b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O03B returns a 3 point quadrature rule for the unit triangle.
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

void trig_unit_o06 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O06 returns a 6 point quadrature rule for the unit triangle.
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

void trig_unit_o06b ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O06B returns a 6 point quadrature rule for the unit triangle.
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

void trig_unit_o07 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O07 returns a 7 point quadrature rule for the unit triangle.
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

void trig_unit_o12 ( double w[], double xy[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_O12 returns a 12 point quadrature rule for the unit triangle.
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

void trig_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_QUAD_TEST tests the rules for the unit triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 April 2008
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
# define DIM_NUM 2

  int dim;
  int dim_num = DIM_NUM;
  int expon[DIM_NUM];
  int h;
  bool more;
  int order;
  double quad;
  int t;
  double *v;
  double *w;
  double *xy;

  cout << "\n";
  cout << "TRIG_UNIT_QUAD_TEST\n";
  cout << "  For the unit triangle,\n";
  cout << "  we approximate monomial integrals with:\n";
  cout << "  TRIG_UNIT_O01,\n";
  cout << "  TRIG_UNIT_O03,\n";
  cout << "  TRIG_UNIT_O03b,\n";
  cout << "  TRIG_UNIT_O06,\n";
  cout << "  TRIG_UNIT_O06b,\n";
  cout << "  TRIG_UNIT_O07,\n";
  cout << "  TRIG_UNIT_O12,\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    order = 1;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o01 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 3;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o03 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 3;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o03b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 6;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o06 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 6;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o06b ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 7;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o07 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    order = 12;
    w = new double[order];
    xy = new double[dim_num*order];
    trig_unit_o12 ( w, xy );
    v = monomial_value ( dim_num, order, expon, xy );
    quad = trig_unit_volume ( ) * r8vec_dot_product ( order, w, v );
    cout << "  " << setw(6) << order
         << "  " << setw(14) << quad << "\n";
    delete [] v;
    delete [] w;
    delete [] xy;

    cout << "\n";
    quad = trig_unit_monomial ( expon );
    cout << "  " << " Exact"
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

double trig_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRIG_UNIT_VOLUME returns the "volume" of the unit triangle in 2D.
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
//    Output, double TRIG_UNIT_VOLUME, the volume.
//
{
  double volume;

  volume = 1.0 / 2.0;

  return volume;
}
//****************************************************************************80

double wedg_unit_monomial ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_MONOMIAL: monomial integral in a unit wedge.
//
//  Discussion:
//
//    This routine returns the integral of
//
//      product ( 1 <= I <= 3 ) X(I)^EXPON(I)
//
//    over the unit wedge.
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 March 2008
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
//    Input, int EXPON[3], the exponents.
//
//    Output, double WEDG_UNIT_MONOMIAL, the integral of the monomial.
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
//
//  Now account for integration in Z.
//
  if ( expon[2] == - 1 )
  {
    cout << "\n";
    cout << "WEDG_UNIT_MONOMIAL - Fatal error!\n";
    cout << "  EXPON[2] = -1 is not a legal input.\n";
    exit ( 1 );
  }
  else if ( ( expon[2] % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = value * 2.0 / ( double ) ( expon[2] + 1 );
  }

  return value;
}
//****************************************************************************80

void wedg_unit_monomial_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_MONOMIAL_TEST tests WEDG_UNIT_MONOMIAL.
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
//  Parameters:
//
//    Input, int DEGREE_MAX, the maximum total degree of the
//    monomials to check.
//
{
  int alpha;
  int beta;
  int expon[3];
  int gamma;
  double value;

  cout << "\n";
  cout << "WEDG_UNIT_MONOMIAL_TEST\n";
  cout << "  For the unit wedge,\n";
  cout << "  WEDG_UNIT_MONOMIAL returns the exact value of the\n";
  cout << "  integral of X^ALPHA Y^BETA Z^GAMMA\n";
  cout << "\n";
  cout << "  Volume = " << wedg_unit_volume ( ) << "\n";
  cout << "\n";
  cout << "     ALPHA      BETA     GAMMA      INTEGRAL\n";
  cout << "\n";

  for ( alpha = 0; alpha <= degree_max; alpha++ )
  {
    expon[0] = alpha;
    for ( beta = 0; beta <= degree_max - alpha; beta++ )
    {
      expon[1] = beta;
      for ( gamma = 0; gamma <= degree_max - alpha - beta; gamma++ )
      {
        expon[2] = gamma;

        value = wedg_unit_monomial ( expon );

        cout << "  " << setw(8)  << expon[0]
             << "  " << setw(8)  << expon[1]
             << "  " << setw(8)  << expon[2]
             << "  " << setw(14) << value << "\n";
      }
    }
  }

  return;
}
//****************************************************************************80

void wedg_unit_quad_test ( int degree_max )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_QUAD_TEST tests the rules for the unit wedge.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2009
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
# define DIM_NUM 3
# define TEST_NUM 7

  int dim;
  int dim_num = DIM_NUM;
  int expon[dim_num];
  int h;
  int line_order;
  int line_order_array[TEST_NUM] = { 1, 2, 2, 3, 2, 3, 4 };
  bool more;
  int order;
  double quad;
  int t;
  int test;
  int trig_order;
  int trig_order_array[TEST_NUM] = { 1, 3, -3, 6, -6, 7, 12 };
  double *v;
  double *w;
  double *xyz;

  cout << "\n";
  cout << "WEDG_UNIT_QUAD_TEST\n";
  cout << "  For the unit wedge,\n";
  cout << "  we approximate monomial integrals with WEDG_UNIT_RULE.\n";

  more = false;

  for ( ; ; )
  {
    subcomp_next ( degree_max, dim_num, expon, &more, &h, &t );

    if ( ( expon[2] % 2 ) == 1 )
    {
      if ( !more )
      {
        break;
      }
      else
      {
        continue;
      }
    }

    cout << "\n";
    cout << "  Monomial exponents: ";
    for ( dim = 0; dim < dim_num; dim++ )
    {
      cout << "  " << setw(2) << expon[dim];
    }
    cout << "\n";
    cout << "\n";

    for ( test = 0; test < TEST_NUM; test++ )
    {
      line_order = line_order_array[test];
      trig_order = trig_order_array[test];

      order = line_order * abs ( trig_order );

      w = new double[order];
      xyz = new double[dim_num*order];
      wedg_unit_rule ( line_order, trig_order, w, xyz );
      v = monomial_value ( dim_num, order, expon, xyz );
      quad = wedg_unit_volume ( ) * r8vec_dot_product ( order, w, v );
      cout << "  " << setw(6) << trig_order
           << "  " << setw(6) << line_order
           << "  " << setw(6) << order
           << "  " << setw(14) << quad << "\n";
      delete [] v;
      delete [] w;
      delete [] xyz;
    }
    cout << "\n";
    quad = wedg_unit_monomial ( expon );
    cout << "   Exact"
         << "                " 
         << "  " << setw(14) << quad << "\n";

    if ( !more )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
//****************************************************************************80

void wedg_unit_rule ( int line_order, int trig_order, double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_RULE returns a quadrature rule for the unit wedge.
//
//  Discussion:
//
//    It is usually sensible to take LINE_ORDER and TRIG_ORDER so that
//    the line and triangle rules are roughly the same precision.  For that
//    criterion, we recommend the following combinations:
//
//      TRIG_ORDER  LINE_ORDER  Precision
//      ----------  ----------  ---------
//          1           1       1 x 1 
//          3           2       2 x 3
//         -3           2       2 x 3
//          6           3       4 x 5
//         -6           2       3 x 3
//          7           3       5 x 5
//         12           4       6 x 7
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    20 April 2009
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
//    Input, int LINE_ORDER, the index of the line rule.
//    The index of the rule is equal to the order of the rule.
//    1 <= LINE_ORDER <= 5.
//
//    Input, int TRIG_ORDER, the indes of the triangle rule.
//    The index of the rule is 1, 3, -3, 6, -6, 7 or 12.
//
//    Output, double W[LINE_ORDER*abs(TRIG_ORDER)], the weights.
//
//    Output, double XYZ[3*LINE_ORDER*abs(TRIG_ORDER)], the abscissas.
//
{
  int i;
  int j;
  int k;
  double *line_w;
  double *line_x;
  double *trig_w;
  double *trig_xy;

  line_w = new double[line_order];
  line_x = new double[line_order];

  if ( line_order == 1 )
  {
    line_unit_o01 ( line_w, line_x );
  }
  else if ( line_order == 2 )
  {
    line_unit_o02 ( line_w, line_x );
  }
  else if ( line_order == 3 )
  {
    line_unit_o03 ( line_w, line_x );
  }
  else if ( line_order == 4 )
  {
    line_unit_o04 ( line_w, line_x );
  }
  else if ( line_order == 5 )
  {
    line_unit_o05 ( line_w, line_x );
  }
  else
  {
    cout << "\n";
    cout << "WEDG_UNIT_RULE - Fatal error!\n";
    cout << "  Illegal value of LINE_ORDER.\n";
    exit ( 1 );
  }

  trig_w = new double[abs(trig_order)];
  trig_xy = new double[2*abs(trig_order)];

  if ( trig_order == 1 )
  {
    trig_unit_o01 ( trig_w, trig_xy );
  }
  else if ( trig_order == 3 )
  {
    trig_unit_o03 ( trig_w, trig_xy );
  }
  else if ( trig_order == - 3 )
  {
    trig_unit_o03b ( trig_w, trig_xy );
  }
  else if ( trig_order == 6 )
  {
    trig_unit_o06 ( trig_w, trig_xy );
  }
  else if ( trig_order == - 6 )
  {
    trig_unit_o06b ( trig_w, trig_xy );
  }
  else if ( trig_order == 7 )
  {
    trig_unit_o07 ( trig_w, trig_xy );
  }
  else if ( trig_order == 12 )
  {
    trig_unit_o12 ( trig_w, trig_xy );
  }
  else
  {
    cout << "\n";
    cout << "WEDG_UNIT_RULE - Fatal error!\n";
    cout << "  Illegal value of TRIG_ORDER.\n";
    exit ( 1 );
  }

  k = 0;
  for ( i = 0; i < line_order; i++ )
  {
    for ( j = 0; j < abs ( trig_order ); j++ )
    {
      w[k] = line_w[i] * trig_w[j];
      xyz[0+k*3] = trig_xy[0+j*2];
      xyz[1+k*3] = trig_xy[1+j*2];
      xyz[2+k*3] = line_x[i];
      k = k + 1;
    }
  }

  delete [] line_w;
  delete [] line_x;
  delete [] trig_w;
  delete [] trig_xy;

  return;
}
//****************************************************************************80

double wedg_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_VOLUME: volume of a unit wedge.
//
//  Discussion:
//
//    The integration region is:
//
//      0 <= X
//      0 <= Y
//      X + Y <= 1
//      -1 <= Z <= 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 April 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double WEDG_UNIT_VOLUME, the volume.
//
{
  double value;

  value = 1.0;

  return value;
}
//****************************************************************************80

void wedg_unit_write_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    WEDG_UNIT_WRITE_TEST writes out some rules for the unit wedge.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 July 2009
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define RULE_NUM 7

  int dim_num = DIM_NUM;
  int line_order;
  int line_order_array[RULE_NUM] = { 1, 2, 2, 3, 2, 3, 4 };
  int order;
  int rule;
  int trig_order;
  int trig_order_array[RULE_NUM] = { 1, 3, -3, 6, -6, 7, 12 };
  double *w;
  string w_filename;
  double *x;
  string x_filename;

  cout << "\n";
  cout << "WEDG_UNIT_WRITE_TEST\n";
  cout << "  For the unit wedge,\n";
  cout << "  write some rules to a file.\n";
  cout << "\n";
  cout << "   Rule  Trig    Line   Total  W_File X_File\n";
  cout << "         Order   Order  Order\n";
  cout << "\n";

  for ( rule = 0; rule < RULE_NUM; rule++ )
  {
    if ( rule == 0 )
    {
      w_filename = "wedge_felippa_1x1_w.txt";
      x_filename = "wedge_felippa_1x1_x.txt";
    }
    else if ( rule == 1 )
    {
      w_filename = "wedge_felippa_3x2_w.txt";
      x_filename = "wedge_felippa_3x2_x.txt";
    }
    else if ( rule == 2 )
    {
      w_filename = "wedge_felippa_3bx2_w.txt";
      x_filename = "wedge_felippa_3bx2_x.txt";
    }
    else if ( rule == 3 )
    {
      w_filename = "wedge_felippa_6x3_w.txt";
      x_filename = "wedge_felippa_6x3_x.txt";
    }
    else if ( rule == 4 )
    {
      w_filename = "wedge_felippa_6bx2_w.txt";
      x_filename = "wedge_felippa_6bx2_x.txt";
    }
    else if ( rule == 5 )
    {
      w_filename = "wedge_felippa_7x3_w.txt";
      x_filename = "wedge_felippa_7x3_x.txt";
    }
    else if ( rule == 6 )
    {
      w_filename = "wedge_felippa_12x4_w.txt";
      x_filename = "wedge_felippa_12x4_x.txt";
    }

    line_order = line_order_array[rule];
    trig_order = trig_order_array[rule];

    order = line_order * abs ( trig_order );

    w = new double[order];
    x = new double[dim_num*order];
    wedg_unit_rule ( line_order, trig_order, w, x );
    r8mat_write ( w_filename, 1, order, w );
    r8mat_write ( x_filename, dim_num, order, x );
    cout << "  " << setw(6) << rule
         << "  " << setw(6) << trig_order
         << "  " << setw(6) << line_order
         << "  " << setw(6) << order
         << "  " << w_filename 
         << "  " << x_filename << "\n";
    delete [] w;
    delete [] x;
  }

  return;
# undef DIM_NUM
# undef TEST_NUM
}
