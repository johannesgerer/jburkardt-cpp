# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cstring>
# include <cmath>

using namespace std;

# include "pyramid_felippa_rule.hpp"

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

double pyramid_unit_monomial ( int expon[3] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_MONOMIAL: monomial integral in a unit pyramid.
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
//    Output, double PYRAMID__UNIT_MONOMIAL, the volume of the pyramid.
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

void pyramid_unit_o01 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O01 returns a 1 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o05 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O05 returns a 5 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o06 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O06 returns a 6 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o08 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O08 returns an 8 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o08b ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O08B returns an 8 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o09 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O09 returns a 9 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o13 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O13 returns a 13 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o18 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O18 returns an 18 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o27 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O27 returns a 27 point quadrature rule for the unit pyramid.
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

void pyramid_unit_o48 ( double w[], double xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_O48 returns a 48 point quadrature rule for the unit pyramid.
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

double pyramid_unit_volume ( )

//****************************************************************************80
//
//  Purpose:
//
//    PYRAMID__UNIT_VOLUME: volume of a unit pyramid with square base.
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
//    Output, double PYRAMID__UNIT_VOLUME, the volume of the pyramid.
//
{
  double volume;

  volume = 4.0 / 3.0;

  return volume;
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

