# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "bins.hpp"

//****************************************************************************80

void bin_search_one_2d ( int bin[2], int nset, double pset[], int nbin[2],
  int bin_start[], int bin_next[], double ptest[2], bool *found_a_neighbor,
  int *i_min, double *d_min_sq, int *compares )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_SEARCH_ONE_2D searches one cell in a 2D array of bins.
//
//  Discussion:
//
//    The bins are presumed to have been set up by successive calls to:
//
//      R82VEC_BIN_EVEN2,
//      R82VEC_BINNED_REORDER, and
//      R82VEC_BINNED_SORT_A.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int BIN[2], the indices of the cell to be examined.
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[2*NSET], the coordinates of the points in the set.
//
//    Input, int NBIN[2], the number of cells in the horizontal and
//    vertical directions.
//
//    Input, int BIN_START[NBIN[1]*NBIN[2]], BIN_LAST(NBIN(1),NBIN(2)),
//    indicates the index of the first and last element in the bin, or -1
//    if none.
//
//    Input, int BIN_NEXT[NSET], the index of the next element of the bin
//    containing this element.
//
//    Input, double PTEST[2], the coordinates of the test point.
//
//    Input/output, bool *FOUND_A_NEIGHBOR, is set to TRUE if at least
//    one point of PSET is found in the current bin.  Otherwise, it retains its
//    input value.
//
//    Input/output, int *I_MIN, the index of the nearest neighbor in
//    PSET to PTEST, if at least one neighbor has been found.
//
//    Input/output, double *D_MIN_SQ, the square of the distance from the nearest
//    neighbor in PSET to PTEST, if at least one neighbor has been found.
//
//    Input/output, int *COMPARES, the number of elements of PSET whose
//    distance to PTEST has been computed.
//
{
# define NDIM 2

  double d_sq;
  int i;
  int node;

  node = bin_start[bin[0]+bin[1]*nbin[0]];

  while ( 0 < node )
  {
    *found_a_neighbor = true;

    d_sq = 0.0;
    for ( i = 0; i < NDIM; i++ )
    {
      d_sq = d_sq + ( ptest[i] - pset[i+node*NDIM] ) *
                    ( ptest[i] - pset[i+node*NDIM] );
    }
    *compares = *compares + 1;

    if ( d_sq < *d_min_sq )
    {
      *d_min_sq = d_sq;
      *i_min = node;
    }

    node = bin_next[node];

  }

  return;
# undef NDIM
}
//****************************************************************************80

void bin_to_r8_even ( int nbin, int bin, double a, double b, double *cmin,
  double *cmax )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R8_EVEN returns the limits for a given "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
//    An initial bin takes everything less than A, and a final bin takes
//    everything greater than B.
//
//  Example:
//
//    NBIN = 7, A = 10, B = 20
//
//    BIN      CMIN  CMAX
//
//    1         -HUGE 10
//    2         10    12
//    3         12    14
//    4         14    16
//    5         16    18
//    6         18    20
//    7         20    HUGE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.  NBIN is normally at least 3.
//    If NBIN is 1 or 2, then everything is assigned to bin 1.
//
//    Input, int BIN, the index of the bin to be considered.
//    If BIN is less than 1, or greater than NBIN, the user will get what
//    the user deserves.
//
//    Input, double A, B, the lower and upper limits of the bin interval.
//    While A is expected to be less than B, the code should return useful
//    results if A is actually greater than B.
//
//    Output, double *CMIN, *CMAX, the minimum and maximum limits on the bin.
//
{
//
//  Take care of special cases.
//
  if ( nbin <= 2 )
  {
    *cmin = -r8_huge ( );
    *cmax =  r8_huge ( );
    return;
  }

  if ( b == a )
  {
    *cmin = -r8_huge ( );
    *cmax =  r8_huge ( );
    return;
  }
//
//  Compute the bin limits.
//
  if ( bin == 1 )
  {
    *cmin = -r8_huge ( );
    *cmax = a;
  }
  else if ( bin < nbin )
  {
    *cmin = ( ( double ) ( nbin - bin     ) * a
            + ( double ) (        bin - 2 ) * b )
            / ( double ) ( nbin       - 2 );

    *cmax = ( ( double ) ( nbin - bin - 1 ) * a
            + ( double ) (        bin - 1 ) * b )
            / ( double ) ( nbin       - 2 );
  }
  else if ( bin == nbin )
  {
    *cmin = b;
    *cmax = r8_huge ( );
  }
  else
  {
    *cmin = -r8_huge ( );
    *cmax =  r8_huge ( );
  }

  return;
}
//****************************************************************************80

void bin_to_r8_even2 ( int nbin, int bin, double a, double b,
  double *cmin, double *cmax )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R8_EVEN2 returns the limits for a given "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN equal subintervals or bins.
//
//  Example:
//
//    NBIN = 5, A = 10, B = 20
//
//    BIN      CMIN  CMAX
//
//    1         10    12
//    2         12    14
//    3         14    16
//    4         16    18
//    5         18    20
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.
//
//    Input, int BIN, the index of the bin to be considered.
//    If BIN is less than 1, or greater than NBIN, the user will get what
//    the user deserves.
//
//    Input, double A, B, the lower and upper limits of the bin interval.
//    While A is expected to be less than B, the code should return useful
//    results if A is actually greater than B.
//
//    Output, double *CMIN, *CMAX, the minimum and maximum limits on the bin.
//
{
//
//  Compute the bin limits.
//
  if ( bin < 1 )
  {
    *cmin = - r8_huge ( );
    *cmax = a;
  }
  else if ( bin <= nbin )
  {
    *cmin = ( ( double ) ( nbin - bin + 1 ) * a + ( double ) ( bin - 1 ) * b )
      / ( double ) ( nbin );
    *cmax = ( ( double ) ( nbin - bin ) * a + ( double ) ( bin ) * b )
      / ( double ) ( nbin );
  }
  else if ( nbin < bin )
  {
    *cmin = b;
    *cmax = r8_huge ( );
  }

  return;
}
//****************************************************************************80

void bin_to_r82_even ( int nbin, int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R82_EVEN returns the limits for a given R82 "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
//    An initial bin takes everything less than A, and a final bin takes
//    everything greater than B.
//
//  Example:
//
//    NBIN = 7, A(1) = 5, B(1) = 15
//              A(2) = 0, B(2) = 20
//
//     BIN         CMIN      CMAX
//    ------   -----------  --------
//    1, 1     -HUGE -HUGE   5     0
//    2, 2       5     0     7     4
//    3, 3       7     4     9     8
//    4, 4       9     8    11    12
//    5, 5      11    12    13    16
//    6, 6      13    16    15    20
//    7, 7      15    20    HUGE HUGE
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.  NBIN is normally
//    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
//
//    Input, int BIN[2], the index of the bin to be considered.
//    If BIN(I) is less than 1, or greater than NBIN, the user will get what
//    the user deserves.
//
//    Input, double A[2], B[2], the lower and upper limits of the
//    bin interval.  While A(I) is expected to be less than B(I), the code
//    should return useful results if A(I) is actually greater than B(I).
//
//    Output, double CMIN[2], CMAX[2], the minimum and maximum
//    limits on the bin.
//
{
  double cmax_1d;
  double cmin_1d;
  int i;

  for ( i = 0; i < 2; i++ )
  {
    bin_to_r8_even ( nbin, bin[i], a[i], b[i], &cmin_1d, &cmax_1d );
    cmin[i] = cmin_1d;
    cmax[i] = cmax_1d;
  }

  return;
}
//****************************************************************************80

void bin_to_r82_even2 ( int nbin, int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R82_EVEN2 returns the limits for a given R82 "bin" in [A,B].
//
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN equal subintervals or bins.
//
//  Example:
//
//    NBIN = 5, A(1) = 5, B(1) = 15
//              A[2] = 0, B[2] = 20
//
//     BIN         CMIN      CMAX
//    ------   -----------  --------
//    1, 1       5     0     7     4
//    2, 2       7     4     9     8
//    3, 3       9     8    11    12
//    4, 4      11    12    13    16
//    5, 5      13    16    15    20
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.
//
//    Input, int BIN[2], the index of the bin to be considered.
//
//    Input, double A[2], B[2], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Output, double CMIN[2], CMAX[2], the minimum and maximum limits on the bin.
//
{
# define NDIM 2

  double ai;
  double bi;
  double bini;
  double cmaxi;
  double cmini;
  int i;

  for ( i = 0; i < NDIM; i++ )
  {
    bin_to_r8_even2 ( nbin, bin[i], a[i], b[i], &cmin[i], &cmax[i] );
  }

  return;
# undef NDIM
}
//****************************************************************************80

void bin_to_r82_even3 ( int nbin[2], int bin[2], double a[2], double b[2],
  double cmin[2], double cmax[2] )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R82_EVEN3 returns the limits for a given R82 "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A(I) to B(I) is divided into NBIN(I) equal
//    subintervals or bins.
//
//  Example:
//
//    NBIN = (/ 4, 5, /)
//
//    A(1) = 5, B(1) = 15
//    A[2] = 0, B[2] = 20
//
//     BIN         CMIN      CMAX
//    ------   -----------  --------
//    1, 1       5     0     7     4
//    2, 2       7     4     9     8
//    3, 3       9     8    11    12
//    4, 4      11    12    13    16
//    5, 5      13    16    15    20
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN[2], the number of bins in each dimension.
//
//    Input, int BIN[2], the index of the bin to be considered.
//
//    Input, double A[2], B[2], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Output, double CMIN[2], CMAX[2], the minimum and maximum limits on the bin.
//
{
# define NDIM 2

  int i;

  for ( i = 0; i < NDIM; i++ )
  {
    bin_to_r8_even2 ( nbin[i], bin[i], a[i], b[i], &cmin[i], &cmax[i] );
  }

  return;
# undef NDIM
}
//****************************************************************************80

void bin_to_r83_even2 ( int nbin, int bin[3], double a[3], double b[3],
  double cmin[3], double cmax[3] )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R83_EVEN2 returns the limits for a given R83 "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN equal subintervals or bins.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.
//
//    Input, int BIN[3], the index of the bin to be considered.
//
//    Input, double A[3], B[3], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Output, double CMIN[3], CMAX[3], the minimum and maximum limits on the bin.
//
{
# define NDIM 3

  int i;

  for ( i = 0; i < NDIM; i++ )
  {
    bin_to_r8_even2 ( nbin, bin[i], a[i], b[i], &cmin[i], &cmax[i] );
  }

  return;
# undef NDIM
}
//****************************************************************************80

void bin_to_r83_even3 ( int nbin[3], int bin[3], double a[3], double b[3],
  double cmin[3], double cmax[3] )

//****************************************************************************80
//
//  Purpose:
//
//    BIN_TO_R83_EVEN3 returns the limits for a given R83 "bin" in [A,B].
//
//  Discussion:
//
//    The interval from A(I) to B(I) is divided into NBIN(I) equal
//    subintervals or bins.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN[3], the number of bins in each dimension.
//
//    Input, int BIN[3], the index of the bin to be considered.
//
//    Input, double A[3], B[3], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Output, double CMIN[3], CMAX[3], the minimum and maximum limits on the bin.
//
{
# define NDIM 3

  int i;

  for ( i = 0; i < NDIM; i++ )
  {
    bin_to_r8_even2 ( nbin[i], bin[i], a[i], b[i], &cmin[i], &cmax[i] );
  }

  return;
# undef NDIM
}
//****************************************************************************80

int diaedg ( double x0, double y0, double x1, double y1, double x2, double y2,
  double x3, double y3 )

//****************************************************************************80
//
//  Purpose:
//
//    DIAEDG chooses a diagonal edge.
//
//  Discussion:
//
//    The routine determines whether 0--2 or 1--3 is the diagonal edge
//    that should be chosen, based on the circumcircle criterion, where
//    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
//    quadrilateral in counterclockwise order.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double X0, Y0, X1, Y1, X2, Y2, X3, Y3, the coordinates of the
//    vertices of a quadrilateral, given in counter clockwise order.
//
//    Output, int DIAEDG, chooses a diagonal:
//    +1, if diagonal edge 02 is chosen;
//    -1, if diagonal edge 13 is chosen;
//     0, if the four vertices are cocircular.
//
{
  double ca;
  double cb;
  double dx10;
  double dx12;
  double dx30;
  double dx32;
  double dy10;
  double dy12;
  double dy30;
  double dy32;
  double s;
  double tol;
  double tola;
  double tolb;
  int value;

  tol = 100.0 * r8_epsilon ( );

  dx10 = x1 - x0;
  dy10 = y1 - y0;
  dx12 = x1 - x2;
  dy12 = y1 - y2;
  dx30 = x3 - x0;
  dy30 = y3 - y0;
  dx32 = x3 - x2;
  dy32 = y3 - y2;

  tola = tol * r8_max ( fabs ( dx10 ),
               r8_max ( fabs ( dy10 ),
               r8_max ( fabs ( dx30 ), fabs ( dy30 ) ) ) );

  tolb = tol * r8_max ( fabs ( dx12 ),
               r8_max ( fabs ( dy12 ),
               r8_max ( fabs ( dx32 ), fabs ( dy32 ) ) ) );

  ca = dx10 * dx30 + dy10 * dy30;
  cb = dx12 * dx32 + dy12 * dy32;

  if ( tola < ca && tolb < cb )
  {
    value = -1;
  }
  else if ( ca < -tola && cb < -tolb )
  {
    value = 1;
  }
  else
  {
    tola = r8_max ( tola, tolb );
    s = ( dx10 * dy30 - dx30 * dy10 ) * cb
      + ( dx32 * dy12 - dx12 * dy32 ) * ca;

    if ( tola < s )
    {
      value = -1;
    }
    else if ( s < -tola )
    {
      value = 1;
    }
    else
    {
      value = 0;
    }

  }

  return value;
}
//****************************************************************************80

int get_seed ( )

//****************************************************************************80
//
//  Purpose:
//
//    GET_SEED returns a random seed for the random number generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int GET_SEED, a random seed value.
//
{
# define I4_MAX 2147483647

  time_t clock;
  int i;
  int ihour;
  int imin;
  int isec;
  struct tm *lt;
  int seed;
  static int seed_internal = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed_internal == 0 )
  {

    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Hours is 1, 2, ..., 12.
//
    ihour = lt->tm_hour;

    if ( 12 < ihour )
    {
      ihour = ihour - 12;
    }
//
//  Move Hours to 0, 1, ..., 11
//
    ihour = ihour - 1;

    imin = lt->tm_min;

    isec = lt->tm_sec;

    seed_internal = isec + 60 * ( imin + 60 * ihour );
//
//  We want values in [1,43200], not [0,43199].
//
    seed_internal = seed_internal + 1;
//
//  Remap SEED from [1,43200] to [1,I4_MAX].
//
    seed_internal = ( int ) (
      ( double ) seed_internal * ( double ) I4_MAX
      / ( 60.0 * 60.0 * 12.0 )
    );

  }
//
//  Never use a seed of 0.
//
  if ( seed_internal == 0 )
  {
    seed_internal = 1;
  }

  if ( seed_internal == I4_MAX )
  {
    seed_internal = I4_MAX - 1;
  }

  seed = seed_internal;

  return seed;
# undef I4_MAX
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
//
{
  if ( i2 < i1 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the smaller of two I4's.
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
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

int i4_modp ( int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Discussion:
//
//    If
//      NREM = I4_MODP ( I, J )
//      NMULT = ( I - NREM ) / J
//    then
//      I = J * NMULT + NREM
//    where NREM is always nonnegative.
//
//    The MOD function computes a result with the same sign as the
//    quantity being divided.  Thus, suppose you had an angle A,
//    and you wanted to ensure that it was between 0 and 360.
//    Then mod(A,360) would do, if A was positive, but if A
//    was negative, your result would be between -360 and 0.
//
//    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
//
//  Example:
//
//        I         J     MOD  I4_MODP   I4_MODP Factorization
//
//      107        50       7       7    107 =  2 *  50 + 7
//      107       -50       7       7    107 = -2 * -50 + 7
//     -107        50      -7      43   -107 = -3 *  50 + 43
//     -107       -50      -7      43   -107 =  3 * -50 + 43
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
{
  int value;

  if ( j == 0 )
  {
    cout << "\n";
    cout << "I4_MODP - Fatal error!\n";
    cout << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}
//****************************************************************************80

int i4_sign ( int i )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SIGN returns the sign of an I4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the integer whose sign is desired.
//
//    Output, int I4_SIGN, the sign of I.
{
  if ( i < 0 )
  {
    return (-1);
  }
  else
  {
    return 1;
  }

}
//****************************************************************************80

void i4_swap ( int *i, int *j )

//****************************************************************************80
//
//  Purpose:
//
//    I4_SWAP switches two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 January 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int *I, *J.  On output, the values of I and
//    J have been interchanged.
//
{
  int k;

  k = *i;
  *i = *j;
  *j = k;

  return;
}
//****************************************************************************80

int i4_uniform ( int a, int b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM returns a scaled pseudorandom I4.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 November 2006
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
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Peter Lewis, Allen Goodman, James Miller
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int k;
  float r;
  int value;

  if ( *seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( float ) ( *seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) ( i4_min ( a, b ) ) - 0.5 )
    +         r   * ( ( float ) ( i4_max ( a, b ) ) + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = r4_nint ( r );

  value = i4_max ( value, i4_min ( a, b ) );
  value = i4_min ( value, i4_max ( a, b ) );

  return value;
}
//****************************************************************************80*

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80*
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Example:
//
//    ILO = 4, IHI = 8
//
//    I  I4_WRAP
//
//    -2     8
//    -1     4
//     0     5
//     1     6
//     2     7
//     3     8
//     4     4
//     5     5
//     6     6
//     7     7
//     8     8
//     9     4
//    10     5
//    11     6
//    12     7
//    13     8
//    14     4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
{
  int jhi;
  int jlo;
  int value;
  int wide;

  jlo = i4_min ( ilo, ihi );
  jhi = i4_max ( ilo, ihi );

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}
//****************************************************************************80

void i4mat_print ( int m, int n, int a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT prints an I4MAT, with an optional title.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_PRINT_SOME prints some of an I4MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title for the matrix.
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to INCX) entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void i4mat_transpose_print ( int m, int n, int a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT prints an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, int A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  int i;
  int j;
  int jhi;
  int jlo;

  i4mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void i4mat_transpose_print_some ( int m, int n, int a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4MAT_TRANSPOSE_PRINT_SOME prints some of an I4MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, int A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title for the matrix.
{
# define INCX 10

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of INCX.
//
  for ( i2lo = ilo; i2lo <= ihi; i2lo = i2lo + INCX )
  {
    i2hi = i2lo + INCX - 1;
    i2hi = i4_min ( i2hi, m );
    i2hi = i4_min ( i2hi, ihi );

    cout << "\n";
//
//  For each row I in the current range...
//
//  Write the header.
//
    cout << "  Row:    ";
    for ( i = i2lo; i <= i2hi; i++ )
    {
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
//
//  Print out (up to INCX) entries in column J, that lie in the current strip.
//
      cout << setw(5) << j << "  ";
      for ( i = i2lo; i <= i2hi; i++ )
      {
        cout << setw(6) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void i4vec_heap_d ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_HEAP_D reorders an I4VEC into a descending heap.
//
//  Discussion:
//
//    A heap is an array A with the property that, for every index J,
//    A[J] >= A[2*J+1] and A[J] >= A[2*J+2], (as long as the indices
//    2*J+1 and 2*J+2 are legal).
//
//  Diagram:
//
//                  A(0)
//                /      \
//            A(1)         A(2)
//          /     \        /  \
//      A(3)       A(4)  A(5) A(6)
//      /  \       /   \
//    A(7) A(8)  A(9) A(10)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    A Nijenhuis and H Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the size of the input array.
//
//    Input/output, int A[N].
//    On input, an unsorted array.
//    On output, the array has been reordered into a heap.
//
{
  int i;
  int ifree;
  int key;
  int m;
//
//  Only nodes (N/2)-1 down to 0 can be "parent" nodes.
//
  for ( i = (n/2)-1; 0 <= i; i-- )
  {
//
//  Copy the value out of the parent node.
//  Position IFREE is now "open".
//
    key = a[i];
    ifree = i;

    for ( ;; )
    {
//
//  Positions 2*IFREE + 1 and 2*IFREE + 2 are the descendants of position
//  IFREE.  (One or both may not exist because they equal or exceed N.)
//
      m = 2 * ifree + 1;
//
//  Does the first position exist?
//
      if ( n <= m )
      {
        break;
      }
      else
      {
//
//  Does the second position exist?
//
        if ( m + 1 < n )
        {
//
//  If both positions exist, take the larger of the two values,
//  and update M if necessary.
//
          if ( a[m] < a[m+1] )
          {
            m = m + 1;
          }
        }
//
//  If the large descendant is larger than KEY, move it up,
//  and update IFREE, the location of the free position, and
//  consider the descendants of THIS position.
//
        if ( key < a[m] )
        {
          a[ifree] = a[m];
          ifree = m;
        }
        else
        {
          break;
        }

      }

    }
//
//  When you have stopped shifting items up, return the item you
//  pulled out back to the heap.
//
    a[ifree] = key;

  }

  return;
}
//****************************************************************************80

int *i4vec_indicator ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_INDICATOR sets an I4VEC to the indicator vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of A.
//
//    Output, int I4VEC_INDICATOR(N), the initialized array.
//
{
  int *a;
  int i;

  a = new int[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = i + 1;
  }

  return a;
}
//****************************************************************************80

void i4vec_print ( int n, int a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_PRINT prints an I4VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 February 2003
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
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i << "  " << setw(8) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void i4vec_sort_heap_a ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORT_HEAP_A ascending sorts an I4VEC using heap sort.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    A Nijenhuis and H Wilf,
//    Combinatorial Algorithms,
//    Academic Press, 1978, second edition,
//    ISBN 0-12-519260-6.
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, int A[N].
//    On input, the array to be sorted;
//    On output, the array has been sorted.
//
{
  int n1;
  int temp;

  if ( n <= 1 )
  {
    return;
  }
//
//  1: Put A into descending heap form.
//
  i4vec_heap_d ( n, a );
//
//  2: Sort A.
//
//  The largest object in the heap is in A[0].
//  Move it to position A[N-1].
//
  temp = a[0];
  a[0] = a[n-1];
  a[n-1] = temp;
//
//  Consider the diminished heap of size N1.
//
  for ( n1 = n-1; 2 <= n1; n1-- )
  {
//
//  Restore the heap structure of the initial N1 entries of A.
//
    i4vec_heap_d ( n1, a );
//
//  Take the largest object from A[0] and move it to A[N1-1].
//
    temp = a[0];
    a[0] = a[n1-1];
    a[n1-1] = temp;

  }

  return;
}
//****************************************************************************80

int i4vec_sorted_unique ( int n, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_SORTED_UNIQUE finds the unique elements in a sorted I4VEC.
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
//    24 August 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements in A.
//
//    Input/output, int A[N].  On input, the sorted
//    integer array.  On output, the unique elements in A.
//
//    Output, int I4VEC_SORTED_UNIQUE, the number of unique elements in A.
//
{
  int i;
  int unique_num;

  unique_num = 0;

  if ( n <= 0 )
  {
    return unique_num;
  }

  unique_num = 1;

  for ( i = 1; i < n; i++ )
  {
    if ( a[i] != a[unique_num-1] )
    {
      unique_num = unique_num + 1;
      a[unique_num-1] = a[i];
    }
  }

  return unique_num;
}
//****************************************************************************80

int i4vec2_compare ( int n, int a1[], int a2[], int i, int j )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_COMPARE compares pairs of integers stored in two I4VECs.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of data items.
//
//    Input, int A1[N], A2[N], contain the two components of each item.
//
//    Input, int I, J, the items to be compared.  These values will be
//    1-based indices for the arrays A1 and A2.
//
//    Output, int I4VEC2_COMPARE, the results of the comparison:
//    -1, item I < item J,
//     0, item I = item J,
//    +1, item J < item I.
//
{
  int isgn;

  isgn = 0;

       if ( a1[i-1] < a1[j-1] )
  {
    isgn = -1;
  }
  else if ( a1[i-1] == a1[j-1] )
  {
         if ( a2[i-1] < a2[j-1] )
    {
      isgn = -1;
    }
    else if ( a2[i-1] < a2[j-1] )
    {
      isgn = 0;
    }
    else if ( a2[j-1] < a2[i-1] )
    {
      isgn = +1;
    }
  }
  else if ( a1[j-1] < a1[i-1] )
  {
    isgn = +1;
  }

  return isgn;
}
//****************************************************************************80

void i4vec2_sort_a ( int n, int a1[], int a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
//
//  Discussion:
//
//    Each item to be sorted is a pair of integers (I,J), with the I
//    and J values stored in separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items of data.
//
//    Input/output, int A1[N], A2[N], the data to be sorted..
//
{
  int i;
  int indx;
  int isgn;
  int j;
  int temp;
//
//  Initialize.
//
  i = 0;
  indx = 0;
  isgn = 0;
  j = 0;
//
//  Call the external heap sorter.
//
  for ( ; ; )
  {
    sort_heap_external ( n, &indx, &i, &j, isgn );
//
//  Interchange the I and J objects.
//
    if ( 0 < indx )
    {
      temp    = a1[i-1];
      a1[i-1] = a1[j-1];
      a1[j-1] = temp;

      temp    = a2[i-1];
      a2[i-1] = a2[j-1];
      a2[j-1] = temp;
    }
//
//  Compare the I and J objects.
//
    else if ( indx < 0 )
    {
      isgn = i4vec2_compare ( n, a1, a2, i, j );
    }
    else if ( indx == 0 )
    {
      break;
    }

  }

  return;
}
//****************************************************************************80

void i4vec2_sorted_unique ( int n, int a1[], int a2[], int *nuniq )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC2_SORTED_UNIQUE finds unique elements in a sorted I4VEC2.
//
//  Discussion:
//
//    Item I is stored as the pair A1(I), A2(I).
//
//    The items must have been sorted, or at least it must be the
//    case that equal items are stored in adjacent vector locations.
//
//    If the items were not sorted, then this routine will only
//    replace a string of equal values by a single representative.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2000
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of items.
//
//    Input/output, int A1[N], A2[N].
//    On input, the array of N items.
//    On output, an array of NUNIQ unique items.
//
//    Output, int *NUNIQ, the number of unique items.
//
{
  int itest;

  *nuniq = 0;

  if ( n <= 0 )
  {
    return;
  }

  *nuniq = 1;

  for ( itest = 1; itest < n; itest++ )
  {
    if ( a1[itest] != a1[*nuniq-1] ||
         a2[itest] != a2[*nuniq-1] )
    {
      a1[*nuniq] = a1[itest];
      a2[*nuniq] = a2[itest];
      *nuniq = *nuniq + 1;
    }
  }

  return;
}
//****************************************************************************80

void index_box2_next_2d ( int n1, int n2, int ic, int jc, int *i, int *j,
  int *more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_2D produces indices on the surface of a box in 2D.
//
//  Discussion:
//
//    The box has center at (IC,JC), and has half-widths N1 and N2.
//    The indices are exactly those which are between (IC-N1,JC-N2) and
//    (IC+N1,JC+N2) with the property that at least one of I and J
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, the half-widths of the box, that is, the
//    maximum distance allowed between (IC,JC) and (I,J).
//
//    Input, int IC, JC, the central cell of the box.
//
//    Input/output, int *I, *J.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I and J are ignored.
//
//    Input/output, bool *MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !(*more) )
  {
    *more = true;
    *i = ic - n1;
    *j = jc - n2;
    return;
  }

  if ( *i == ic + n1 &&
       *j == jc + n2 )
  {
    *more = false;
    return;
  }
//
//  Increment J.
//
  *j = *j + 1;
//
//  Check J.
//
  if ( jc + n2 < *j )
  {
    *j = jc - n2;
    *i = *i + 1;
  }
  else if ( *j < jc + n2 && ( *i == ic - n1 || *i == ic + n1 ) )
  {
    return;
  }
  else
  {
    *j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc,
  int *i, int *j, int *k, bool *more )

//****************************************************************************80
//
//  Purpose:
//
//    INDEX_BOX2_NEXT_3D produces indices on the surface of a box in 3D.
//
//  Discussion:
//
//    The box has a central cell of (IC,JC,KC), with a half widths of
//    (N1,N2,N3).  The indices are exactly those between (IC-N1,JC-N2,KC-N3) and
//    (IC+N1,JC+N2,KC+N3) with the property that at least one of I, J, and K
//    is an "extreme" value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the "half widths" of the box, that is, the
//    maximum distances from the central cell allowed for I, J and K.
//
//    Input, int IC, JC, KC, the central cell of the box.
//
//    Input/output, int *I, *J, *K.  On input, the previous index set.
//    On output, the next index set.  On the first call, MORE should
//    be set to FALSE, and the input values of I, J, and K are ignored.
//
//    Input/output, bool *MORE.
//    On the first call for a given box, the user should set MORE to FALSE.
//    On return, the routine sets MORE to TRUE.
//    When there are no more indices, the routine sets MORE to FALSE.
//
{
  if ( !(*more) )
  {
    *more = true;
    *i = ic - n1;
    *j = jc - n2;
    *k = kc - n3;
    return;
  }

  if ( *i == ic + n1 &&
       *j == jc + n2 &&
       *k == kc + n3 )
  {
    *more = false;
    return;
  }
//
//  Increment K.
//
  *k = *k + 1;
//
//  Check K.
//
  if ( kc + n3 < *k )
  {
    *k = kc - n3;
    *j = *j + 1;
  }
  else if ( *k < kc + n3 &&
    ( *i == ic - n1 ||
      *i == ic + n1 ||
      *j == jc - n2 ||
      *j == jc + n2 ) )
  {
    return;
  }
  else
  {
    *k = kc + n3;
    return;
  }
//
//  Check J.
//
  if ( jc + n2 < *j )
  {
    *j = jc - n2;
    *i = *i + 1;
  }
  else if ( *j < jc + n2 &&
    ( *i == ic - n1 ||
      *i == ic + n1 ||
      *k == kc - n3 ||
      *k == kc + n3 ) )
  {
    return;
  }
  else
  {
    *j = jc + n2;
    return;
  }

  return;
}
//****************************************************************************80

int lrline ( double xu, double yu, double xv1, double yv1, double xv2,
  double yv2, double dv )

//****************************************************************************80
//
//  Purpose:
//
//    LRLINE determines where a point lies in relation to a directed line.
//
//  Discussion:
//
//    LRLINE determines whether a point is to the left of, right of,
//    or on a directed line parallel to a line through given points.
//
//  Modified:
//
//    22 June 2009
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
//    directed line is parallel to and at signed distance DV to the left of
//    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
//    which the position relative to the directed line is to be determined.
//
//    Input, double DV, the signed distance, positive for left.
//
//    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
//    to the right of, on, or left of the directed line.  LRLINE is 0 if
//    the line degenerates to a point.
//
{
  double dx;
  double dxu;
  double dy;
  double dyu;
  double t;
  double tol;
  double tolabs;
  int value;

  tol = 100.0 * r8_epsilon ( );

  dx = xv2 - xv1;
  dy = yv2 - yv1;
  dxu = xu - xv1;
  dyu = yu - yv1;

  tolabs = tol * r8_max ( fabs ( dx ),
                 r8_max ( fabs ( dy ),
                 r8_max ( fabs ( dxu ),
                 r8_max ( fabs ( dyu ), fabs ( dv ) ) ) ) );

  t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy );

  if ( tolabs < t )
  {
    value = 1;
  }
  else if ( -tolabs <= t )
  {
    value = 0;
  }
  else if ( t < -tolabs )
  {
    value = -1;
  }

  return value;
}
//****************************************************************************80

bool perm_check ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_CHECK checks that a vector represents a permutation.
//
//  Discussion:
//
//    The routine verifies that each of the integers from 1
//    to N occurs among the N entries of the permutation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries.
//
//    Input, int P[N], the array to check.
//
//    Output, bool PERM_CHECK, is TRUE if the permutation is OK.
//
{
  bool found;
  int i;
  int seek;

  for ( seek = 1; seek <= n; seek++ )
  {
    found = false;

    for ( i = 0; i < n; i++ )
    {
      if ( p[i] == seek )
      {
        found = true;
        break;
      }
    }

    if ( !found )
    {
      return false;
    }

  }

  return true;
}
//****************************************************************************80

void perm_inv ( int n, int p[] )

//****************************************************************************80
//
//  Purpose:
//
//    PERM_INV inverts a permutation "in place".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects being permuted.
//
//    Input/output, int P[N], the permutation, in standard index form.
//    On output, P describes the inverse permutation
//
{
  int i;
  int i0;
  int i1;
  int i2;
  int is;

  if ( n <= 0 )
  {
    cout << "\n";
    cout << "PERM_INV - Fatal error!\n";
    cout << "  Input value of N = " << n << "\n";
    exit ( 1 );
  }

  if ( !perm_check ( n, p ) )
  {
    cout << "\n";
    cout << "PERM_INV - Fatal error!\n";
    cout << "  The input array does not represent\n";
    cout << "  a proper permutation.\n";
    exit ( 1 );
  }

  is = 1;

  for ( i = 1; i <= n; i++ )
  {
    i1 = p[i-1];

    while ( i < i1 )
    {
      i2 = p[i1-1];
      p[i1-1] = -i2;
      i1 = i2;
    }

    is = - i4_sign ( p[i-1] );
    p[i-1] = i4_sign ( is ) * abs ( p[i-1] );
  }

  for ( i = 1; i <= n; i++ )
  {
    i1 = -p[i-1];

    if ( 0 <= i1 )
    {
      i0 = i;

      for ( ; ; )
      {
        i2 = p[i1-1];
        p[i1-1] = i0;

        if ( i2 < 0 )
        {
          break;
        }
        i0 = i1;
        i1 = i2;
      }
    }
  }

  return;
}
//****************************************************************************80

int points_nearest_point_naive_2d ( int nset, double pset[], double ptest[],
  double *d_min )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_NEAREST_POINT_NAIVE_2D finds the nearest point to a given point in 2D.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[2*NSET], the coordinates of the points in the set.
//
//    Input, double PTEST[2], the point whose nearest neighbor is sought.
//
//    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
//
//    Output, int POINTS_NEAREST_POINT_NAIVE_2D, the index of the nearest
//    point in PSET to P.
//
{
# define NDIM 2

  double d;
  int i;
  int j;
  int p_min;

  *d_min = r8_huge ( );
  p_min = 0;

  for ( j = 0; j < nset; j++ )
  {
    d = 0.0;
    for ( i = 0; i < NDIM; i++ )
    {
      d = d + ( ptest[i] - pset[i+j*NDIM] ) * ( ptest[i] - pset[i+j*NDIM] );
    }
    if ( d < *d_min )
    {
      *d_min = d;
      p_min = j;
    }
  }

  *d_min = sqrt ( *d_min );

  return p_min;

# undef NDIM
}
//****************************************************************************80

int points_nearest_point_naive_3d ( int nset, double pset[], double ptest[],
  double *d_min )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_NEAREST_POINT_NAIVE_3D finds the nearest point to a given point in 3D.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[3*NSET], the coordinates of the points in the set.
//
//    Input, double PTEST[3], the point whose nearest neighbor is sought.
//
//    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
//
//    Output, int POINTS_NEAREST_POINT_NAIVE_3D, the index of the nearest
//    point in PSET to P.
//
{
# define NDIM 3

  double d;
  int i;
  int j;
  int p_min;

  *d_min = r8_huge ( );
  p_min = 0;

  for ( j = 0; j < nset; j++ )
  {
    d = 0.0;
    for ( i = 0; i < NDIM; i++ )
    {
      d = d + ( ptest[i] - pset[i+j*NDIM] ) * ( ptest[i] - pset[i+j*NDIM] );
    }
    if ( d < *d_min )
    {
      *d_min = d;
      p_min = j;
    }
  }

  *d_min = sqrt ( *d_min );

  return p_min;

# undef NDIM
}
//****************************************************************************80

int points_nearest_point_naive_nd ( int ndim, int nset, double pset[],
  double ptest[], double *d_min )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_NEAREST_POINT_NAIVE_ND finds the nearest point to a given point in ND.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[NDIM*NSET], the coordinates of the points in the set.
//
//    Input, double PTEST[NDIM], the point whose nearest neighbor is sought.
//
//    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
//
//    Output, int POINTS_NEAREST_POINT_NAIVE_ND, the index of the nearest
//    point in PSET to P.
//
{
  double d;
  int i;
  int j;
  int p_min;

  *d_min = r8_huge ( );
  p_min = 0;

  for ( j = 0; j < nset; j++ )
  {
    d = 0.0;
    for ( i = 0; i < ndim; i++ )
    {
      d = d + ( ptest[i] - pset[i+j*ndim] ) * ( ptest[i] - pset[i+j*ndim] );
    }
    if ( d < *d_min )
    {
      *d_min = d;
      p_min = j;
    }
  }

  *d_min = sqrt ( *d_min );

  return p_min;
}
//****************************************************************************80

int *points_nearest_points_naive_2d ( int nset, double pset[], int ntest,
  double ptest[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_NEAREST_POINTS_NAIVE_2D finds the nearest point to given points in 2D.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[2*NSET], the coordinates of the points in the set.
//
//    Input, int NTEST, the number of test points.
//
//    Input, double PTEST[2*NTEST], the coordinates of the test points.
//
//    Output, int POINTS_NEAREST_POINTS_NAIVE_2D[NTEST], the index of the
//    nearest point in PSET to each point in PTEST.
//
{
# define NDIM 2

  double d;
  double d_min;
  int i;
  int *nearest;
  int set;
  int test;

  nearest = new int[ntest];

  for ( test = 0; test < ntest; test++ )
  {
    d_min = r8_huge ( );
    nearest[test] = -1;

    for ( set = 0; set < nset; set++ )
    {
      d = 0.0;
      for ( i = 0; i < NDIM; i++ )
      {
        d = d + ( ptest[i,test] - pset[i,set] )
              * ( ptest[i,test] - pset[i,set] );
      }

      if ( d < d_min )
      {
        d_min = d;
        nearest[test] = set;
      }
    }
  }

  return nearest;
# undef NDIM
}
//****************************************************************************80

int *points_nearest_points_naive_3d ( int nset, double pset[], int ntest,
  double ptest[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINTS_NEAREST_POINTS_NAIVE_3D finds the nearest point to given points in 3D.
//
//  Discussion:
//
//    A naive algorithm is used.  The distance to every point is calculated,
//    in order to determine the smallest.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NSET, the number of points in the set.
//
//    Input, double PSET[3*NSET], the coordinates of the points in the set.
//
//    Input, int NTEST, the number of test points.
//
//    Input, double PTEST[3*NTEST], the coordinates of the test points.
//
//    Output, int POINTS_NEAREST_POINTS_NAIVE_3D[NTEST], the index of the
//    nearest point in PSET to each point in PTEST.
//
{
# define NDIM 3

  double d;
  double d_min;
  int i;
  int *nearest;
  int set;
  int test;

  nearest = new int[ntest];

  for ( test = 0; test < ntest; test++ )
  {
    d_min = r8_huge ( );
    nearest[test] = -1;

    for ( set = 0; set < nset; set++ )
    {
      d = 0.0;
      for ( i = 0; i < NDIM; i++ )
      {
        d = d + ( ptest[i,test] - pset[i,set] )
              * ( ptest[i,test] - pset[i,set] );
      }

      if ( d < d_min )
      {
        d_min = d;
        nearest[test] = set;
      }
    }
  }

  return nearest;
# undef NDIM
}
//****************************************************************************80

int r4_nint ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NINT returns the nearest integer to an R4.
//
//  Example:
//
//        X         R4_NINT
//
//      1.3         1
//      1.4         1
//      1.5         1 or 2
//      1.6         2
//      0.0         0
//     -0.7        -1
//     -1.1        -1
//     -1.6        -2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the value.
//
//    Output, int R4_NINT, the nearest integer to X.
//
{
  int s;

  if ( x < 0.0 )
  {
    s = -1;
  }
  else
  {
    s = 1;
  }

  return ( s * ( int ) ( fabs ( x ) + 0.5 ) );
}
//****************************************************************************80

double r8_add ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_ADD adds two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the numbers to be added.
//
//    Output, double R8_ADD, the sum of X and Y.
//
{
  double value;

  value = x + y;

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
//    11 August 2010
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
  double one;
  double temp;
  double test;
  double value;

  one = ( double ) ( 1 );

  value = one;
  temp = value / 2.0;
  test = r8_add ( one, temp );

  while ( one < test )
  {
    value = temp;
    temp = value / 2.0;
    test = r8_add ( one, temp );
  }


  return value;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8 value.
//
//  Discussion:
//
//    HUGE_VAL is the largest representable legal double precision number, and is usually
//    defined in math.h, or sometimes in stdlib.h.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" double precision value.
//
{
  return HUGE_VAL;
}
//****************************************************************************80

double r8_max ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
//****************************************************************************80

int r8_to_bin_even ( int nbin, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_BIN_EVEN determines the appropriate "bin" for C in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN-2 equal subintervals or bins.
//    An initial bin takes everything less than A, and a final bin takes
//    everything greater than B.
//
//  Example:
//
//    NBIN = 7, A = 5, B = 15
//
//    C   BIN
//
//    1    1
//    3    1
//    4.9  1
//    5    2
//    6    2
//    7    3
//    8    3
//    9.5  4
//   13    6
//   14    6
//   15    6
//   15.1  7
//   99    7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 January 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.  NBIN is normally
//    at least 3.  If NBIN is 1 or 2, then everything is assigned to bin 1.
//
//    Input, double A, B, the lower and upper limits of the bin
//    interval.  While A is expected to be less than B, the code should
//    return useful results if A is actually greater than B.
//
//    Input, double C, a value to be placed in a bin.
//
//    Output, inte R8_TO_BIN_EVEN, the index of the bin to which C is
//    assigned.
//
{
  double a2;
  double b2;
  int bin;
  bool swap;
//
//  Take care of special cases.
//
  if ( nbin < 1 )
  {
    bin = 0;
    return bin;
  }
  else if ( nbin == 1 || nbin == 2 )
  {
    bin = 1;
    return bin;
  }

  if ( b == a )
  {
    bin = 0;
    return bin;
  }
//
//  If the limits are descending, then we switch them now, and
//  unswitch the results at the end.
//
  if ( a < b )
  {
    swap = false;
    a2 = a;
    b2 = b;
  }
  else
  {
    swap = true;
    a2 = b;
    b2 = a;
  }
//
//  Compute the bin.
//
  if ( c < a2 )
  {
    bin = 1;
  }
  else if ( c == a2 )
  {
    bin = 2;
  }
  else if ( c == b2 )
  {
    bin = nbin - 1;
  }
  else if ( b2 < c )
  {
    bin = nbin;
  }
  else
  {
    bin = 2 + ( int ) ( ( double ) ( nbin - 2 ) * ( c - a2 ) / ( b2 - a2 ) );
    bin = i4_max ( bin, 2 );
    bin = i4_min ( bin, nbin - 1 );
  }
//
//  Reverse the switching.
//
  if ( swap )
  {
    bin = nbin + 1 - bin;
  }
  return bin;
}
//****************************************************************************80

int r8_to_bin_even2 ( int nbin, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    R8_TO_BIN_EVEN2 determines the appropriate "bin" for C in [A,B].
//
//  Discussion:
//
//    The interval from A to B is divided into NBIN equal subintervals or bins.
//
//  Example:
//
//    NBIN = 5, A = 5, B = 15
//
//    <-1-+-2-+-3-+-4-+-5->
//    5   7   9  11  13  15
//
//    C   BIN
//
//    1    1
//    3    1
//    4.9  1
//    5    1
//    6    1
//    7.1  2
//    8    2
//    9.5  3
//   12    4
//   14    5
//   15    5
//   15.1  5
//   99    5
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins.
//
//    Input, double A, B, the lower and upper limits of the bin interval.
//    While A is expected to be less than B, the code should return useful
//    results if A is actually greater than B.
//
//    Input, double C, a value to be placed in a bin.
//
//    Output, int R_TO_BIN_EVEN_2, the index of the bin to which C is assigned.
//
{
  double a2;
  double b2;
  int bin;
  bool reorder;
//
//  Take care of special cases.
//
  if ( nbin < 1 )
  {
    bin = 0;
    return bin;
  }

  if ( nbin == 1 )
  {
    bin = 1;
    return bin;
  }

  if ( b == a )
  {
    bin = 1;
    return bin;
  }
//
//  If the limits are descending, then we switch them now, and
//  unswitch the results at the end.
//
  if ( a < b )
  {
    reorder = false;
    a2 = a;
    b2 = b;
  }
  else
  {
    reorder = true;
    a2 = b;
    b2 = a;
  }
//
//  Compute the bin.
//
  if ( c <= a2 )
  {
    bin = 1;
  }
  else if ( b2 <= c )
  {
    bin = nbin;
  }
  else
  {
    bin = 1 + ( int ) ( ( ( double ) nbin ) * ( c - a2 ) / ( b2 - a2 ) );
    bin = i4_max ( bin, 1 );
    bin = i4_min ( bin, nbin );
  }
//
//  Reverse the switching.
//
  if ( reorder )
  {
    bin = nbin + 1 - bin;
  }

  return bin;
}
//****************************************************************************80

double r8_uniform ( double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM returns a scaled pseudorandom R8.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A, B, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8_UNIFORM, a number strictly between A and B.
//
{
  double value;

  value = a + ( b - a ) * r8_uniform_01 ( seed );

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      r8_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

int *r82_to_bin_even2 ( int nbin, double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82_TO_BIN_EVEN2 determines the appropriate "bin" for an R82 value.
//
//  Discussion:
//
//    The intervals [A(1),B(1)] and [A(2),B(2)] are each divided into NBIN
//    equal subintervals or bins.  Boundary bins take care of extreme values.
//
//  Example:
//
//    NBIN = 5, A(1) = 5,  A(2) = 0,
//              B(1) = 15, B(2) = 20.
//
//   20 +    +    +    +    +    +
//        15 | 25 | 35 | 45 | 55
//   16 +----+----+----+----+----+
//        14 | 24 | 34 | 44 | 54
//   12 +----+----+----+----+----+
//        13 | 23 | 33 | 43 | 53
//    8 +----+----+----+----+----+
//        12 | 22 | 32 | 42 | 52
//    4 +----+----+----+----+----+
//        11 | 21 | 31 | 41 | 51
//    0 +    +    +    +    +    +
//      5    7    9   11   13   15
//
//      C      BIN
//   ------  ------
//    8 -2    2  1
//    0  1    1  1
//    6  9    1  3
//   10 11    3  3
//   14 23    5  5
//   25 13    5  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins in each dimension.
//    NBIN must be at least 1.
//
//    Input, double A[2], B[2], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Input, double C[2], a value to be placed in a bin.
//
//    Output, int R82_TO_BIN_EVEN2[2], the index of the bin to which C is assigned.
//
{
  int *bin;
  int i;

  bin = new int[2];

  for ( i = 0; i < 2; i++ )
  {
    bin[i] = r8_to_bin_even2 ( nbin, a[i], b[i], c[i] );
  }

  return bin;
}
//****************************************************************************80

int *r82_to_bin_even3 ( int nbin[], double a[], double b[], double c[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82_TO_BIN_EVEN3 determines the appropriate "bin" for an R82 value.
//
//  Discussion:
//
//    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
//    or bins.
//
//  Example:
//
//    NBIN = (/ 4, 5 /),
//
//      A(1) = 1,  A(2) = 0,
//      B(1) = 17, B(2) = 20.
//
//   20 +    +    +    +    +
//        15 | 25 | 35 | 45
//   16 +----+----+----+----+
//        14 | 24 | 34 | 44
//   12 +----+----+----+----+
//        13 | 23 | 33 | 43
//    8 +----+----+----+----+
//        12 | 22 | 32 | 42
//    4 +----+----+----+----+
//        11 | 21 | 31 | 41
//    0 +    +    +    +    +
//      1    5    9   13   17
//
//      C      BIN
//   ------  ------
//    8 -2    2  1
//    0  1    1  1
//    6  9    2  3
//   10 11    3  3
//   14 23    4  5
//   25 13    4  4
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN[2], the number of bins in each dimension.
//
//    Input, double A[2], B[2], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Input, double C[2], a value to be placed in a bin.
//
//    Output, int R82_TO_BIN_EVEN3[2], the bin assignment for the value.
//
{
  int *bin;
  int i;

  bin = new int[2];

  for ( i = 0; i < 2; i++ )
  {
    bin[i] = r8_to_bin_even2 ( nbin[i], a[i], b[i], c[i] );
  }

  return bin;
}
//****************************************************************************80

void r82_uniform ( double rlo[], double rhi[], int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82_UNIFORM returns a random R82 value in a given range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double RLO[2], RHI[2], the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R[2], the randomly chosen value.
//
{
  int i;

  for ( i = 0; i < 2; i++ )
  {
    r[i] = rlo[i] + ( rhi[i] - rlo[i] ) * r8_uniform_01 ( seed );
  }

  return;
}
//****************************************************************************80

void r82vec_part_quick_a ( int n, double a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PART_QUICK_A reorders an R82 vector as part of a quick sort.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 2, stored as a vector
//    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
//
//    The routine reorders the entries of A.  Using A(1:2,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Example:
//
//    Input:
//
//      N = 8
//
//      A = ( (2,4), (8,8), (6,2), (0,2), (10,6), (10,0), (0,6), (4,8) )
//
//    Output:
//
//      L = 2, R = 4
//
//      A = ( (0,2), (0,6), (2,4), (8,8), (6,2), (10,6), (10,0), (4,8) )
//             -----------          ----------------------------------
//             LEFT          KEY    RIGHT
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, double A[N*2].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:2,1).  Then
//    I <= L                 A(1:2,I) < KEY;
//         L < I < R         A(1:2,I) = KEY;
//                 R <= I    A(1:2,I) > KEY.
//
{
  int i;
  int j;
  double key[2];
  int ll;
  int m;
  int rr;
//
  if ( n < 1 )
  {
    cout << "\n";
    cout << "R82VEC_PART_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[2*0+0];
  key[1] = a[2*0+1];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 2, a+2*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 2, a+2*(rr-1), a+2*ll );
    }
    else if ( r8vec_eq ( 2, a+2*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 2, a+2*(m-1), a+2*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 2, a+2*ll, key ) )
    {
      ll = ll + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = a[2*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
//****************************************************************************80*

void r82vec_permute ( int n, double a[], int p[] )

//****************************************************************************80*
//
//  Purpose:
//
//    R82VEC_PERMUTE permutes an R82VEC in place.
//
//  Discussion:
//
//    This routine permutes an array of real "objects", but the same
//    logic can be used to permute an array of objects of any arithmetic
//    type, or an array of objects of any complexity.  The only temporary
//    storage required is enough to store a single object.  The number
//    of data movements made is N + the number of cycles of order 2 or more,
//    which is never more than N + N/2.
//
//  Example:
//
//    Input:
//
//      N = 5
//      P = (   2,    4,    5,    1,    3 )
//      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
//          (11.0, 22.0, 33.0, 44.0, 55.0 )
//
//    Output:
//
//      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
//             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of objects.
//
//    Input/output, double A[2*N], the array to be permuted.
//
//    Input, int P[N], the permutation.  P(I) = J means
//    that the I-th element of the output array should be the J-th
//    element of the input array.  P must be a legal permutation
//    of the integers from 1 to N, otherwise the algorithm will
//    fail catastrophically.
//
{
  double a_temp[2];
  int i;
  int iget;
  int iput;
  int istart;

  if ( !perm_check ( n, p ) )
  {
    cout << "\n";
    cout << "R82VEC_PERMUTE - Fatal error!\n";
    cout << "  The input array does not represent\n";
    cout << "  a proper permutation.\n";
    i4vec_print ( n, p, "  The faulty permutation:" );
    exit ( 1 );
  }
//
//  Search for the next element of the permutation that has not been used.
//
  for ( istart = 1; istart <= n; istart++ )
  {
    if ( p[istart-1] < 0 )
    {
      continue;
    }
    else if ( p[istart-1] == istart )
    {
      p[istart-1] = -p[istart-1];
      continue;
    }
    else
    {
      a_temp[0] = a[0+(istart-1)*2];
      a_temp[1] = a[1+(istart-1)*2];
      iget = istart;
//
//  Copy the new value into the vacated entry.
//
      for ( ; ; )
      {
        iput = iget;
        iget = p[iget-1];

        p[iput-1] = -p[iput-1];

        if ( iget < 1 || n < iget )
        {
          cout << "\n";
          cout << "R82VEC_PERMUTE - Fatal error!\n";
          cout << "  Entry IPUT = " << iput << " of the permutation has\n";
          cout << "  an illegal value IGET = " << iget << ".\n";
          exit ( 1 );
        }

        if ( iget == istart )
        {
          a[0+(iput-1)*2] = a_temp[0];
          a[1+(iput-1)*2] = a_temp[1];
          break;
        }
        a[0+(iput-1)*2] = a[0+(iget-1)*2];
        a[1+(iput-1)*2] = a[1+(iget-1)*2];
      }
    }
  }
//
//  Restore the signs of the entries.
//
  for ( i = 0; i < n; i++ )
  {
    p[i] = -p[i];
  }

  return;
}
//****************************************************************************80

void r82vec_print ( int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_PRINT prints an R82VEC.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 2, stored as a vector
//    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N*2], the vector to be printed.
//
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6)  << i        << "  "
         << setw(14) << a[2*i+0] << "  "
         << setw(14) << a[2*i+1] << "\n";
  }

  return;
}
//****************************************************************************80

int *r82vec_sort_heap_index_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R82VEC.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(1:2,INDX(I)), I = 1 to N is sorted,
//
//    or explicitly, by the call
//
//      call R82VEC_PERMUTE ( N, A, INDX )
//
//    after which A(1:2,I), I = 1 to N is sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[2*N], an array to be index-sorted.
//
//    Output, int R82VEC_SORT_HEAP_INDEX_A[N], the sort index.  The
//    I-th element of the sorted array is A(0:1,R82VEC_SORT_HEAP_INDEX_A(I-1)).
//
{
  double aval[2];
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;
//
  if ( n < 1 )
  {
    return NULL;
  }

  if ( n == 1 )
  {
    indx = new int[1];
    indx[0] = 1;
    return indx;
  }

  indx = i4vec_indicator ( n );

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval[0] = a[0+(indxt-1)*2];
      aval[1] = a[1+(indxt-1)*2];
    }
    else
    {
      indxt = indx[ir-1];
      aval[0] = a[0+(indxt-1)*2];
      aval[1] = a[1+(indxt-1)*2];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        break;
      }

    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if (   a[0+(indx[j-1]-1)*2] <  a[0+(indx[j]-1)*2] ||
             ( a[0+(indx[j-1]-1)*2] == a[0+(indx[j]-1)*2] &&
               a[1+(indx[j-1]-1)*2] <  a[1+(indx[j]-1)*2] ) )
        {
          j = j + 1;
        }
      }

      if (   aval[0] <  a[0+(indx[j-1]-1)*2] ||
           ( aval[0] == a[0+(indx[j-1]-1)*2] &&
             aval[1] <  a[1+(indx[j-1]-1)*2] ) )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }

  return indx;
}
//****************************************************************************80

void r82vec_sort_quick_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_QUICK_A ascending sorts an R82VEC using quick sort.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 2, stored as a vector
//    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N*2].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 25

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "R82VEC_SORT_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r82vec_part_quick_a ( n_segment, a+2*(base-1)+0, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cout << "\n";
        cout << "R82VEC_SORT_QUICK_A - Fatal error!\n";
        cout << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }

      }

    }

  }
  return;
# undef LEVEL_MAX
}
//****************************************************************************80

void r82vec_uniform ( int n, double alo[], double ahi[], int *seed,
  double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_UNIFORM returns a random R82VEC in a given range.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 2, stored as a vector
//    of rows: A(0,0), A(0,1), // A(1,0), A(1,1) // ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double ALO[2], AHI[2], the range allowed for the entries.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double A[N*2], the vector of randomly chosen integers.
//
{
  int i;
  int j;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < 2; j++ )
    {
      a[2*i+j] = r8_uniform ( alo[j], ahi[j], seed );
    }
  }

  return;
}
//****************************************************************************80

int *r83_to_bin_even2 ( int nbin, double a[3], double b[3], double c[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_TO_BIN_EVEN2 determines the appropriate "bin" for an R83.
//
//  Discussion:
//
//    The intervals [A(I),B(I)] are each divided into NBIN
//    equal subintervals or bins.  Boundary bins take care of extreme values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN, the number of bins in each dimension.
//    NBIN must be at least 1.
//
//    Input, double A[3], B[3], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Input, double C[3], a value to be placed in a bin.
//
//    Output, int R83_TO_BIN_EVEN2[3], the index of the bin to which C is assigned.
//
{
  int *bin;
  int i;

  bin = new int[3];

  for ( i = 0; i < 3; i++ )
  {
    bin[i] = r8_to_bin_even2 ( nbin, a[i], b[i], c[i] );
  }

  return bin;
}
//****************************************************************************80

int *r83_to_bin_even3 ( int nbin[3], double a[3], double b[3], double c[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R83_TO_BIN_EVEN3 determines the appropriate "bin" for an R83.
//
//  Discussion:
//
//    The interval [A(I),B(I)] is divided into NBIN(I) equal subintervals
//    or bins.
//
//  Example:
//
//    NBIN = (/ 4, 5, 2 /),
//
//      A(1) = 1,  A(2) = 0,  A(3) = 8
//      B(1) = 17, B(2) = 20, B(3) = 10
//
//
//            8 < Z < 9                    9 < Z < 10
//
//   20 +     +     +     +     +     20 +     +     +     +     +
//        151 | 251 | 351 | 451            152 | 252 | 352 | 452
//   16 +-----+-----+-----+-----+     16 +-----+-----+-----+-----+
//        141 | 241 | 341 | 441            142 | 242 | 342 | 442
//   12 +-----+-----+-----+-----+     12 +-----+-----+-----+-----+
//        131 | 231 | 331 | 431            132 | 232 | 332 | 432
//    8 +-----+-----+-----+-----+      8 +-----+-----+-----+-----+
//        121 | 221 | 321 | 421            122 | 222 | 322 | 422
//    4 +-----+-----+-----+-----+      4 +-----+-----+-----+-----+
//        111 | 211 | 311 | 411            112 | 212 | 312 | 412
//    0 +     +     +     +     +      0 +     +     +     +     +
//      1     5     9    13    17        1     5     9    13    17
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NBIN[3], the number of bins in each dimension.
//
//    Input, double A[3], B[3], the lower and upper limits of the bin interval.
//    While A(I) is expected to be less than B(I), the code should return useful
//    results if A(I) is actually greater than B(I).
//
//    Input, double C[3], a value to be placed in a bin.
//
//    Output, int R83_TO_BIN_EVEN3[3], the index of the bin to which C is assigned.
//
{
  int *bin;
  int i;

  bin = new int[3];

  for ( i = 0; i < 3; i++ )
  {
    bin[i] = r8_to_bin_even2 ( nbin[i], a[i], b[i], c[i] );
  }

  return bin;
}
//****************************************************************************80

void r83vec_part_quick_a ( int n, double a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R83VEC_PART_QUICK_A reorders an R83VEC as part of a quick sort.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 3, stored as a vector
//    of rows: A(0,0), A(0,1), A(0,2) // A(1,0), A(1,1), A(1,2) // ...
//
//    The routine reorders the entries of A.  Using A(1:3,1) as a
//    key, all entries of A that are less than or equal to the key will
//    precede the key, which precedes all entries that are greater than the key.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, double A[N*3].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int *L, *R, the indices of A that define the three segments.
//    Let KEY = the input value of A(1:3,1).  Then
//    I <= L                 A(1:3,I) < KEY;
//         L < I < R         A(1:3,I) = KEY;
//                 R <= I    A(1:3,I) > KEY.
//
{
  int i;
  int j;
  double key[3];
  int ll;
  int m;
  int rr;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "R83VEC_PART_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key[0] = a[3*0+0];
  key[1] = a[3*0+1];
  key[2] = a[3*0+2];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  ll = 1;
  rr = n + 1;

  for ( i = 2; i <= n; i++ )
  {
    if ( r8vec_gt ( 3, a+3*ll, key ) )
    {
      rr = rr - 1;
      r8vec_swap ( 3, a+3*(rr-1), a+3*ll );
    }
    else if ( r8vec_eq ( 3, a+3*ll, key ) )
    {
      m = m + 1;
      r8vec_swap ( 3, a+3*(m-1), a+3*ll );
      ll = ll + 1;
    }
    else if ( r8vec_lt ( 3, a+3*ll, key ) )
    {
      ll = ll + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 0; i < ll - m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = a[3*(i+m)+j];
    }
  }

  ll = ll - m;

  for ( i = ll; i < ll+m; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = key[j];
    }
  }

  *l = ll;
  *r = rr;

  return;
}
//****************************************************************************80

void r83vec_sort_quick_a ( int n, double a[] )

//****************************************************************************80*
//
//  Purpose:
//
//    R83VEC_SORT_QUICK_A ascending sorts an R83VEC using quick sort.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 3, stored as a vector
//    of rows: A(0,0), A(0,1), A(0,2) // A(1,0), A(1,1), A(1,2) // ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input/output, double A[N*3].
//    On input, the array to be sorted.
//    On output, the array has been sorted.
//
{
# define LEVEL_MAX 25

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cout << "\n";
    cout << "R83VEC_SORT_QUICK_A - Fatal error!\n";
    cout << "  N < 1.\n";
    exit ( 1 );
  }

  if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[level-1] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r83vec_part_quick_a ( n_segment, a+3*(base-1)+0, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {
      if ( LEVEL_MAX < level )
      {
        cout << "\n";
        cout << "R83VEC_SORT_QUICK_A - Fatal error!\n";
        cout << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( level <= 1 )
        {
          n_segment = 0;
          break;
        }

        base = rsave[level-1];
        n_segment = rsave[level-2] - rsave[level-1];
        level = level - 1;

        if ( 0 < n_segment )
        {
          break;
        }

      }

    }

  }
  return;
# undef LEVEL_MAX
}
//****************************************************************************80

void r83vec_uniform ( int n, double alo[], double ahi[], int *seed, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R83VEC_UNIFORM returns a random R83VEC in a given range.
//
//  Discussion:
//
//    A is a two dimensional array of order N by 3, stored as a vector
//    of rows: A(0,0), A(0,1), A(0,2), // A(1,0), A(1,1), A(1,2) // ...
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double ALO[3], AHI[3], the minimum and maximum values.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double A[N*3], the vector of randomly chosen values.
//
{
  int i;
  int j;

  for ( i = 0; i < n; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      a[3*i+j] = r8_uniform ( alo[j], ahi[j], seed );
    }
  }

  return;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    The doubly dimensioned array A is treated as a one dimensional vector,
//    stored by COLUMNS.  Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, char *TITLE, a title to be printed.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, char *TITLE, a title for the matrix.
{
# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }
//
//  Print the columns of the matrix, in strips of 5.
//
  for ( j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX )
  {
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min ( j2hi, n );
    j2hi = i4_min ( j2hi, jhi );

    cout << "\n";
//
//  For each column J in the current range...
//
//  Write the header.
//
    cout << "  Col:    ";
    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(7) << j << "       ";
    }
    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
//
//  Determine the range of the rows in this strip.
//
    i2lo = i4_max ( ilo, 1 );
    i2hi = i4_min ( ihi, m );

    for ( i = i2lo; i <= i2hi; i++ )
    {
//
//  Print out (up to) 5 entries in row I, that lie in the current strip.
//
      cout << setw(5) << i << "  ";
      for ( j = j2lo; j <= j2hi; j++ )
      {
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }

  }

  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8mat_transpose_print ( int m, int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT prints an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
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
//    Input, char *TITLE, an optional title.
//
{
  r8mat_transpose_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_transpose_print_some ( int m, int n, double a[], int ilo, int jlo,
  int ihi, int jhi, char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
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
//    Input, char *TITLE, an optional title.
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

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

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
      cout << setw(7) << i << "       ";
    }
    cout << "\n";
    cout << "  Col\n";
    cout << "\n";

    j2lo = i4_max ( jlo, 1 );
    j2hi = i4_min ( jhi, n );

    for ( j = j2lo; j <= j2hi; j++ )
    {
      cout << setw(5) << j << " ";
      for ( i2 = 1; i2 <= inc; i2++ )
      {
        i = i2lo - 1 + i2;
        cout << setw(14) << a[(i-1)+(j-1)*m];
      }
      cout << "\n";
    }
  }
  cout << "\n";

  return;
# undef INCX
}
//****************************************************************************80

void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_BRACKET searches a sorted array for successive brackets of a value.
//
//  Discussion:
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
//    Output, int *LEFT, *RIGHT, the results of the search.
//
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}
//****************************************************************************80

bool r8vec_eq ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_EQ is true if every pair of entries in two vectors is equal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], two vectors to compare.
//
//    Output, bool R8VEC_EQ.
//    R8VEC_EQ is TRUE if every pair of elements A1(I) and A2(I) are equal,
//    and FALSE otherwise.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] != a2[i] )
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

bool r8vec_gt ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_GT == ( A1 > A2 ) for real vectors.
//
//  Discussion:
//
//    The comparison is lexicographic.
//
//    A1 > A2  <=>                              A1(1) > A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) > A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) > A2(N)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double A1[N], A2[N], the vectors to be compared.
//
//    Output, bool R8VEC_GT, is TRUE if and only if A1 > A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a2[i] < a1[i] )
    {
       return true;
    }
    else if ( a1[i] < a2[i] )
    {
      return false;
    }
  }

  return false;
}
//****************************************************************************80

bool r8vec_lt ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LT == ( A1 < A2 ) for real vectors.
//
//  Discussion:
//
//    The comparison is lexicographic.
//
//    A1 < A2  <=>                              A1(1) < A2(1) or
//                 ( A1(1)     == A2(1)     and A1(2) < A2(2) ) or
//                 ...
//                 ( A1(1:N-1) == A2(1:N-1) and A1(N) < A2(N)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double A1[N], A2[N], the vectors to be compared.
//
//    Output, bool R8VEC_LT, is TRUE if and only if A1 < A2.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    if ( a1[i] < a2[i] )
    {
      return true;
    }
    else if ( a2[i] < a1[i] )
    {
      return false;
    }

  }

  return false;
}
//****************************************************************************80

void r8vec_part_quick_a ( int n, double a[], int *l, int *r )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PART_QUICK_A reorders an R8VEC as part of a quick sort.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The routine reorders the entries of A.  Using A[0] as a
//    key, all entries of A that are less than or equal to A[0] will
//    precede A[0] which precedes all entries that are greater than A[0].
//
//  Example:
//
//    Input:
//
//  N = 8
//
//  A = ( 6, 7, 3, 1, 6, 8, 2, 9 )
//
//    Output:
//
//  L = 3, R = 6
//
//  A = ( 3, 1, 2, 6, 6, 8, 9, 7 )
//        -------        -------
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, double A[N].  On input, the array to be checked.
//    On output, A has been reordered as described above.
//
//    Output, int L, R, the indices of A that define the three segments.
//    Let KEY = the input value of A[0].  Then
//    I <= L             A(I) < KEY;
//     L < I < R         A(I) = KEY;
//             R <= I    A(I) > KEY.
//
{
  int i;
  double key;
  int m;
  double temp;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_PART_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    *l = 0;
    *r = 2;
    return;
  }

  key = a[0];
  m = 1;
//
//  The elements of unknown size have indices between L+1 and R-1.
//
  *l = 1;
  *r = n + 1;

  for ( i = 2; i <= n; i++ )
  {

    if ( key < a[*l] )
    {
      *r = *r - 1;
      temp = a[*r-1];
      a[*r-1] = a[*l];
      a[*l] = temp;
    }
    else if ( a[*l] == key )
    {
      m = m + 1;
      temp = a[m-1];
      a[m-1] = a[*l];
      a[*l] = temp;
      *l = *l + 1;
    }
    else if ( a[*l] < key )
    {
      *l = *l + 1;
    }

  }
//
//  Now shift small elements to the left, and KEY elements to center.
//
  for ( i = 1; i <= *l -m; i++ )
  {
    a[i-1] = a[i+m-1];
  }

  *l = *l - m;

  for ( i = *l+1; i <= *l+m; i++ )
  {
    a[i-1] = key;
  }

  return;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], char *title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints a R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2003
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
//    Input, char *TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( s_len_trim ( title ) != 0 )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i << "  " << setw(10) << a[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec_sort_quick_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_QUICK_A ascending sorts an R8VEC using quick sort.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    Input:
//
//      N = 7
//
//      A = ( 6, 7, 3, 2, 9, 1, 8 )
//
//    Output:
//
//      A = ( 1, 2, 3, 6, 7, 8, 9 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 April 1999
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries of A.
//
//    Input/output, double A[N].  On input, the array to be sorted.
//    On output, A has been reordered into ascending order.
//
{
# define LEVEL_MAX 30

  int base;
  int l_segment;
  int level;
  int n_segment;
  int rsave[LEVEL_MAX];
  int r_segment;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "R8VEC_SORT_QUICK_A - Fatal error!\n";
    cerr << "  N < 1.\n";
    exit ( 1 );
  }
  else if ( n == 1 )
  {
    return;
  }

  level = 1;
  rsave[0] = n + 1;
  base = 1;
  n_segment = n;

  while ( 0 < n_segment )
  {
//
//  Partition the segment.
//
    r8vec_part_quick_a ( n_segment, a+base-1, &l_segment, &r_segment );
//
//  If the left segment has more than one element, we need to partition it.
//
    if ( 1 < l_segment )
    {

      if ( LEVEL_MAX < level )
      {
        cerr << "\n";
        cerr << "R8VEC_SORT_QUICK_A - Fatal error!\n";
        cerr << "  Exceeding recursion maximum of " << LEVEL_MAX << "\n";
        exit ( 1 );
      }

      level = level + 1;
      n_segment = l_segment;
      rsave[level-1] = r_segment + base - 1;
    }
//
//  The left segment and the middle segment are sorted.
//  Must the right segment be partitioned?
//
    else if ( r_segment < n_segment )
    {
      n_segment = n_segment + 1 - r_segment;
      base = base + r_segment - 1;
    }
//
//  Otherwise, we back up a level if there is an earlier one.
//
    else
    {
      for ( ; ; )
      {
        if ( 1 < level )
        {
          base = rsave[level-1];
          n_segment = rsave[level-2] - rsave[level-1];
          level = level - 1;
          if ( 0 < n_segment )
          {
            break;
          }
        }
        else
        {
          n_segment = 0;
          break;
        }
      }
    }
  }

  return;
# undef LEVEL_MAX
}
//****************************************************************************80

void r8vec_swap ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SWAP swaps the entries of two R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the arrays.
//
//    Input/output, double A1[N], A2[N], the vectors to swap.
//
{
  int i;
  double temp;

  for ( i = 0; i < n; i++ )
  {
    temp  = a1[i];
    a1[i] = a2[i];
    a2[i] = temp;
  }

  return;
}
//****************************************************************************80

double *r8vec_uniform ( int n, double a, double b, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM fills an R8VEC with scaled pseudorandom values.
//
//  Discussion:
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
//    30 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, L E Schrage,
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
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A, B, the lower and upper limits of the pseudorandom values.
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

    r[i] = a + ( b - a ) * ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

int s_len_trim ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn )

//****************************************************************************80
//
//  Purpose:
//
//    SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
//
//  Discussion:
//
//    The actual list is not passed to the routine.  Hence it may
//    consist of integers, reals, numbers, names, etc.  The user,
//    after each return from the routine, will be asked to compare or
//    interchange two items.
//
//    The current version of this code mimics the FORTRAN version,
//    so the values of I and J, in particular, are FORTRAN indices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2004
//
//  Author:
//
//    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the length of the input list.
//
//    Input/output, int *INDX.
//    The user must set INDX to 0 before the first call.
//    On return,
//      if INDX is greater than 0, the user must interchange
//      items I and J and recall the routine.
//      If INDX is less than 0, the user is to compare items I
//      and J and return in ISGN a negative value if I is to
//      precede J, and a positive value otherwise.
//      If INDX is 0, the sorting is done.
//
//    Output, int *I, *J.  On return with INDX positive,
//    elements I and J of the user's list should be
//    interchanged.  On return with INDX negative, elements I
//    and J are to be compared by the user.
//
//    Input, int ISGN. On return with INDX negative, the
//    user should compare elements I and J of the list.  If
//    item I is to precede item J, set ISGN negative,
//    otherwise set ISGN positive.
//
{
  static int i_save = 0;
  static int j_save = 0;
  static int k = 0;
  static int k1 = 0;
  static int n1 = 0;
//
//  INDX = 0: This is the first call.
//
  if ( *indx == 0 )
  {

    i_save = 0;
    j_save = 0;
    k = n / 2;
    k1 = k;
    n1 = n;
  }
//
//  INDX < 0: The user is returning the results of a comparison.
//
  else if ( *indx < 0 )
  {
    if ( *indx == -2 )
    {
      if ( isgn < 0 )
      {
        i_save = i_save + 1;
      }
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( 0 < isgn )
    {
      *indx = 2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      if ( n1 == 1 )
      {
        i_save = 0;
        j_save = 0;
        *indx = 0;
      }
      else
      {
        i_save = n1;
        j_save = 1;
        n1 = n1 - 1;
        *indx = 1;
      }
      *i = i_save;
      *j = j_save;
      return;
    }

    k = k - 1;
    k1 = k;

  }
//
//  0 < INDX: the user was asked to make an interchange.
//
  else if ( *indx == 1 )
  {
    k1 = k;
  }

  for ( ;; )
  {

    i_save = 2 * k1;

    if ( i_save == n1 )
    {
      j_save = k1;
      k1 = i_save;
      *indx = -1;
      *i = i_save;
      *j = j_save;
      return;
    }
    else if ( i_save <= n1 )
    {
      j_save = i_save + 1;
      *indx = -2;
      *i = i_save;
      *j = j_save;
      return;
    }

    if ( k <= 1 )
    {
      break;
    }

    k = k - 1;
    k1 = k;
  }

  if ( n1 == 1 )
  {
    i_save = 0;
    j_save = 0;
    *indx = 0;
    *i = i_save;
    *j = j_save;
  }
  else
  {
    i_save = n1;
    j_save = 1;
    n1 = n1 - 1;
    *indx = 1;
    *i = i_save;
    *j = j_save;
  }

  return;
}
//****************************************************************************80

int swapec ( int i, int *top, int *btri, int *bedg, int point_num,
  double point_xy[], int tri_num, int tri_vert[], int tri_nabe[],
  int stack[] )

//****************************************************************************80
//
//  Purpose:
//
//    SWAPEC swaps diagonal edges until all triangles are Delaunay.
//
//  Discussion:
//
//    The routine swaps diagonal edges in a 2D triangulation, based on
//    the empty circumcircle criterion, until all triangles are Delaunay,
//    given that I is the index of the new vertex added to the triangulation.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Parameters:
//
//    Input, int I, the index of the new vertex.
//
//    Input/output, int *TOP, the index of the top of the stack.
//    On output, TOP is zero.
//
//    Input/output, int *BTRI, *BEDG; on input, if positive, are the
//    triangle and edge indices of a boundary edge whose updated indices
//    must be recorded.  On output, these may be updated because of swaps.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the points.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input/output, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//    May be updated on output because of swaps.
//
//    Input/output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list;
//    negative values are used for links of the counter-clockwise linked
//    list of boundary edges;  May be updated on output because of swaps.
//
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Workspace, int STACK[MAXST]; on input, entries 1 through TOP
//    contain the indices of initial triangles (involving vertex I)
//    put in stack; the edges opposite I should be in interior;  entries
//    TOP+1 through MAXST are used as a stack.
//
//    Output, int SWAPEC, is set to 8 for abnormal return.
//
{
  int a;
  int b;
  int c;
  int e;
  int ee;
  int em1;
  int ep1;
  int f;
  int fm1;
  int fp1;
  int l;
  int r;
  int s;
  int swap;
  int t;
  int tt;
  int u;
  double x;
  double y;
//
//  Determine whether triangles in stack are Delaunay, and swap
//  diagonal edge of convex quadrilateral if not.
//
  x = point_xy[2*(i-1)+0];
  y = point_xy[2*(i-1)+1];

  for ( ; ; )
  {
    if ( *top <= 0 )
    {
      break;
    }

    t = stack[(*top)-1];
    *top = *top - 1;

    if ( tri_vert[3*(t-1)+0] == i )
    {
      e = 2;
      b = tri_vert[3*(t-1)+2];
    }
    else if ( tri_vert[3*(t-1)+1] == i )
    {
      e = 3;
      b = tri_vert[3*(t-1)+0];
    }
    else
    {
      e = 1;
      b = tri_vert[3*(t-1)+1];
    }

    a = tri_vert[3*(t-1)+e-1];
    u = tri_nabe[3*(t-1)+e-1];

    if ( tri_nabe[3*(u-1)+0] == t )
    {
      f = 1;
      c = tri_vert[3*(u-1)+2];
    }
    else if ( tri_nabe[3*(u-1)+1] == t )
    {
      f = 2;
      c = tri_vert[3*(u-1)+0];
    }
    else
    {
      f = 3;
      c = tri_vert[3*(u-1)+1];
    }

    swap = diaedg ( x, y,
      point_xy[2*(a-1)+0], point_xy[2*(a-1)+1],
      point_xy[2*(c-1)+0], point_xy[2*(c-1)+1],
      point_xy[2*(b-1)+0], point_xy[2*(b-1)+1] );

    if ( swap == 1 )
    {
      em1 = i4_wrap ( e - 1, 1, 3 );
      ep1 = i4_wrap ( e + 1, 1, 3 );
      fm1 = i4_wrap ( f - 1, 1, 3 );
      fp1 = i4_wrap ( f + 1, 1, 3 );

      tri_vert[3*(t-1)+ep1-1] = c;
      tri_vert[3*(u-1)+fp1-1] = i;
      r = tri_nabe[3*(t-1)+ep1-1];
      s = tri_nabe[3*(u-1)+fp1-1];
      tri_nabe[3*(t-1)+ep1-1] = u;
      tri_nabe[3*(u-1)+fp1-1] = t;
      tri_nabe[3*(t-1)+e-1] = s;
      tri_nabe[3*(u-1)+f-1] = r;

      if ( 0 < tri_nabe[3*(u-1)+fm1-1] )
      {
        *top = *top + 1;
        stack[(*top)-1] = u;
      }

      if ( 0 < s )
      {
        if ( tri_nabe[3*(s-1)+0] == u )
        {
          tri_nabe[3*(s-1)+0] = t;
        }
        else if ( tri_nabe[3*(s-1)+1] == u )
        {
          tri_nabe[3*(s-1)+1] = t;
        }
        else
        {
          tri_nabe[3*(s-1)+2] = t;
        }

        *top = *top + 1;

        if ( point_num < *top )
        {
          return 8;
        }

        stack[(*top)-1] = t;
      }
      else
      {
        if ( u == *btri && fp1 == *bedg )
        {
          *btri = t;
          *bedg = e;
        }

        l = - ( 3 * t + e - 1 );
        tt = t;
        ee = em1;

        while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
        {
          tt = tri_nabe[3*(tt-1)+ee-1];

          if ( tri_vert[3*(tt-1)+0] == a )
          {
            ee = 3;
          }
          else if ( tri_vert[3*(tt-1)+1] == a )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }

        }

        tri_nabe[3*(tt-1)+ee-1] = l;

      }

      if ( 0 < r )
      {
        if ( tri_nabe[3*(r-1)+0] == t )
        {
          tri_nabe[3*(r-1)+0] = u;
        }
        else if ( tri_nabe[3*(r-1)+1] == t )
        {
          tri_nabe[3*(r-1)+1] = u;
        }
        else
        {
          tri_nabe[3*(r-1)+2] = u;
        }
      }
      else
      {
        if ( t == *btri && ep1 == *bedg )
        {
          *btri = u;
          *bedg = f;
        }

        l = - ( 3 * u + f - 1 );
        tt = u;
        ee = fm1;

        while ( 0 < tri_nabe[3*(tt-1)+ee-1] )
        {
          tt = tri_nabe[3*(tt-1)+ee-1];

          if ( tri_vert[3*(tt-1)+0] == b )
          {
            ee = 3;
          }
          else if ( tri_vert[3*(tt-1)+1] == b )
          {
            ee = 1;
          }
          else
          {
            ee = 2;
          }
        }
        tri_nabe[3*(tt-1)+ee-1] = l;
      }
    }
  }
  return 0;
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

double triangle_area_2d ( double t[2*3] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Discussion:
//
//    If the triangle's vertices are given in counter clockwise order,
//    the area will be positive.  If the triangle's vertices are given
//    in clockwise order, the area will be negative!
//
//    An earlier version of this routine always returned the absolute
//    value of the computed area.  I am convinced now that that is
//    a less useful result!  For instance, by returning the signed
//    area of a triangle, it is possible to easily compute the area
//    of a nonconvex polygon as the sum of the (possibly negative)
//    areas of triangles formed by node 1 and successive pairs of vertices.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.
//
{
  double area;

  area = 0.5 * (
    t[0+0*2] * ( t[1+1*2] - t[1+2*2] ) +
    t[0+1*2] * ( t[1+2*2] - t[1+0*2] ) +
    t[0+2*2] * ( t[1+0*2] - t[1+1*2] ) );

  return area;
}
//****************************************************************************80

void triangle_sample ( double t[2*3], int n, int *seed, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE returns random points in a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
# define DIM_NUM 2

  double alpha;
  double beta;
  int j;
  double r;
  double p12[DIM_NUM];
  double p13[DIM_NUM];

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
    alpha = sqrt ( r );
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
    p12[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+1*2];
    p12[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+1*2];

    p13[0] = ( 1.0 - alpha ) * t[0+0*2] + alpha * t[0+2*2];;
    p13[1] = ( 1.0 - alpha ) * t[1+0*2] + alpha * t[1+2*2];;
//
//  Now choose, uniformly at random, a point on the line L.
//
    beta = r8_uniform_01 ( seed );

    p[0+j*2] = ( 1.0 - beta ) * p12[0] + beta * p13[0];
    p[1+j*2] = ( 1.0 - beta ) * p12[1] + beta * p13[1];
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangulation_print ( int node_num, int triangle_num,
  double node_xy[], int triangle_node[], int triangle_neighbor[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_PRINT prints information defining a triangulation.
//
//  Discussion:
//
//    Triangulations created by R8TRIS2 include extra information encoded
//    in the negative values of TRIANGLE_NEIGHBOR.
//
//    Because some of the nodes counted in NODE_NUM may not actually be
//    used in the triangulation, I needed to compute the true number
//    of vertices.  I added this calculation on 13 October 2001.
//
//    Ernest Fasse pointed out an error in the indexing of VERTEX_LIST,
//    which was corrected on 19 February 2004.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up
//    the triangles.
//
//    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbors
//    on each side.  If there is no triangle neighbor on a particular side,
//    the value of TRIANGLE_NEIGHBOR should be negative.  If the
//    triangulation data was created by R8TRIS2, then there is more
//    information encoded in the negative values.
//
{
# define DIM_NUM 2

  int boundary_num;
  int i;
  int j;
  int k;
  int n1;
  int n2;
  int s;
  int s1;
  int s2;
  bool skip;
  int t;
  int *vertex_list;
  int vertex_num;

  cout << "\n";
  cout << "TRIANGULATION_PRINT\n";
  cout << "  Information defining a triangulation.\n";
  cout << "\n";
  cout << "  The number of nodes is " << node_num << "\n";

  r8mat_transpose_print ( DIM_NUM, node_num, node_xy, "  Node coordinates" );

  cout << "\n";
  cout << "  The number of triangles is " << triangle_num << "\n";
  cout << "\n";
  cout << "  Sets of three nodes are used as vertices of\n";
  cout << "  the triangles.  For each triangle, the nodes\n";
  cout << "  are listed in counterclockwise order.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_node, "  Triangle nodes" );

  cout << "\n";
  cout << "  On each side of a given triangle, there is either\n";
  cout << "  another triangle, or a piece of the convex hull.\n";
  cout << "  For each triangle, we list the indices of the three\n";
  cout << "  neighbors, or (if negative) the codes of the\n";
  cout << "  segments of the convex hull.\n";

  i4mat_transpose_print ( 3, triangle_num, triangle_neighbor,
    "  Triangle neighbors" );
//
//  Determine VERTEX_NUM, the number of vertices.
//
  vertex_list = new int[3*triangle_num];

  k = 0;
  for ( t = 0; t < triangle_num; t++ )
  {
    for ( s = 0; s < 3; s++ )
    {
      vertex_list[k] = triangle_node[s+t*3];
      k = k + 1;
    }
  }

  i4vec_sort_heap_a ( 3*triangle_num, vertex_list );

  vertex_num = i4vec_sorted_unique ( 3*triangle_num, vertex_list );

  delete [] vertex_list;
//
//  Determine the number of boundary points.
//
  boundary_num = 2 * vertex_num - triangle_num - 2;

  cout << "\n";
  cout << "  The number of boundary points is " << boundary_num << "\n";
  cout << "\n";
  cout << "  The segments that make up the convex hull can be\n";
  cout << "  determined from the negative entries of the triangle\n";
  cout << "  neighbor list.\n";
  cout << "\n";
  cout << "     #   Tri  Side    N1    N2\n";
  cout << "\n";

  skip = false;

  k = 0;

  for ( i = 0; i < triangle_num; i++ )
  {
    for ( j = 0; j < 3; j++ )
    {
      if ( triangle_neighbor[j+i*3] < 0 )
      {
        s = -triangle_neighbor[j+i*3];
        t = s / 3;

        if ( t < 1 || triangle_num < t )
        {
          cout << "\n";
          cout << "  Sorry, this data does not use the R8TRIS2\n";
          cout << "  convention for convex hull segments.\n";
          skip = true;
          break;
        }

        s1 = ( s % 3 ) + 1;
        s2 = i4_wrap ( s1+1, 1, 3 );
        k = k + 1;
        n1 = triangle_node[s1-1+(t-1)*3];
        n2 = triangle_node[s2-1+(t-1)*3];
        cout                  << "  "
             << setw(4) << k  << "  "
             << setw(4) << t  << "  "
             << setw(4) << s1 << "  "
             << setw(4) << n1 << "  "
             << setw(4) << n2 << "\n";
      }
    }

    if ( skip )
    {
      break;
    }
  }

  return;
# undef DIM_NUM
}
//****************************************************************************80

void triangulation_sample ( int node_num, double node_xy[],
  int triangle_num, int triangle_node[], int num_ran, int *seed,
  double xd[], int td[] )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGULATION_SAMPLE returns random points in a triangulation.
//
//  Discussion:
//
//    It is assumed that the triangulation consists of a set of non-overlapping
//    triangles.
//
//    The point is chosen uniformly in the area covered by the triangulation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
//    triangles.
//
//    Input, int NUM_RAN, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double XD[2*NUM_RAN], the sample points.
//
//    Output, int TD[NUM_RAN], the triangle to which each sample point
//    belongs.
//
{
  double area;
  double *area_cum;
  double area_total;
  int i;
  int i1;
  int i2;
  int i3;
  int left;
  double r;
  int right;
  double t[2*3];
//
//  Compute the areas of the triangles.
//  Build a cumulative area vector.
//  Convert it to a relative cumulative area vector.
//
  area_cum = new double[triangle_num+1];
  area_cum[0] = 0.0;

  for ( i = 0; i < triangle_num; i++ )
  {
    i1 = triangle_node[0+i*3];
    t[0+0*2] = node_xy[0+i1*2];
    t[1+0*2] = node_xy[1+i1*2];

    i2 = triangle_node[1+i*3];
    t[0+1*2] = node_xy[0+i2*2];
    t[1+1*2] = node_xy[1+i2*2];

    i3 = triangle_node[2+i*3];
    t[0+2*2] = node_xy[0+i3*2];
    t[1+2*2] = node_xy[1+i3*2];

    area_cum[i+1] = area_cum[i] + triangle_area_2d ( t );
  }

  area_total = area_cum[triangle_num];

  for ( i = 0; i <= triangle_num; i++ )
  {
    area_cum[i] = area_cum[i] / area_total;
  }
//
//  Pick random values.  A random value R indicates the corresponding triangle
//  whose cumulative relative area contains R.
//
//  Bracket the random value in the cumulative relative areas,
//  indicating a triangle.
//
//  Pick a random point in the triangle.
//
  for ( i = 0; i < num_ran; i++ )
  {
    r = r8_uniform_01 ( seed );

    r8vec_bracket ( triangle_num+1, area_cum, r, &left, &right );

    td[i] = right - 1;

    i1 = triangle_node[0+(td[i]-1)*3];
    t[0+0*2] = node_xy[0+i1*2];
    t[1+0*2] = node_xy[1+i1*2];

    i2 = triangle_node[1+(td[i]-1)*3];
    t[0+1*2] = node_xy[0+i2*2];
    t[1+1*2] = node_xy[1+i2*2];

    i3 = triangle_node[2+(td[i]-1)*3];
    t[0+2*2] = node_xy[0+i3*2];
    t[1+2*2] = node_xy[1+i3*2];

    triangle_sample ( t, 1, seed, xd+i*2 );
  }

  delete [] area_cum;

  return;
}
//****************************************************************************80

void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT2 computes the next element of an integer tuple space.
//
//  Discussion:
//
//    The elements X are N vectors.
//
//    Each entry X(I) is constrained to lie between XMIN(I) and XMAX(I).
//
//    The elements are produced one at a time.
//
//    The first element is
//      (XMIN(1), XMIN(2), ..., XMIN(N)),
//    the second is (probably)
//      (XMIN(1), XMIN(2), ..., XMIN(N)+1),
//    and the last element is
//      (XMAX(1), XMAX(2), ..., XMAX(N))
//
//    Intermediate elements are produced in a lexicographic order, with
//    the first index more important than the last, and the ordering of
//    values at a fixed index implicitly defined by the sign of
//    XMAX(I) - XMIN(I).
//
//  Example:
//
//    N = 2,
//    XMIN = (/ 1, 10 /)
//    XMAX = (/ 3,  8 /)
//
//    RANK    X
//    ----  -----
//      1   1 10
//      2   1  9
//      3   1  8
//      4   2 10
//      5   2  9
//      6   2  8
//      7   3 10
//      8   3  9
//      9   3  8
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
//    Input, int N, the number of components.
//
//    Input, int XMIN[N], XMAX[N], the "minimum" and "maximum" entry values.
//    These values are minimum and maximum only in the sense of the lexicographic
//    ordering.  In fact, XMIN(I) may be less than, equal to, or greater
//    than XMAX(I).
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
//    Input/output, int *RANK, the rank of the item.  On first call,
//    set RANK to 0 to start up the sequence.  On return, if RANK is zero,
//    there are no more items in the sequence.
//
{
  int i;
  int test;

  if ( *rank < 0 )
  {
    cout << "\n";
    cout << "TUPLE_NEXT2 - Fatal error!\n";
    cout << "  Illegal value of RANK = " << *rank << "\n";
    exit ( 1 );
  }

  test = 1;
  for ( i = 0; i < n; i++ )
  {
    test = test * ( 1 + abs ( xmax[i] - xmin[i] ) );
  }

  if ( test < *rank )
  {
    cout << "\n";
    cout << "TUPLE_NEXT2 - Fatal error!\n";
    cout << "  Illegal value of RANK = " << *rank << "\n";
    exit ( 1 );
  }

  if ( *rank == 0 )
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = xmin[i];
    }
    *rank = 1;
    return;
  }

  *rank = *rank + 1;
  i = n - 1;

  for ( ; ; )
  {
    if ( x[i] != xmax[i] )
    {
      if ( xmin[i] < xmax[i] )
      {
        x[i] = x[i] + 1;
      }
      else
      {
        x[i] = x[i] - 1;
      }
      break;
    }

    x[i] = xmin[i];

    if ( i == 0 )
    {
      *rank = 0;
      break;
    }

    i = i - 1;

  }

  return;
}
//****************************************************************************80

void vbedg ( double x, double y, int point_num, double point_xy[], int tri_num,
  int tri_vert[], int tri_nabe[], int *ltri, int *ledg, int *rtri, int *redg )

//****************************************************************************80
//
//  Purpose:
//
//    VBEDG determines which boundary edges are visible to a point.
//
//  Discussion:
//
//    The point (X,Y) is assumed to be outside the convex hull of the
//    region covered by the 2D triangulation.
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Barry Joe,
//    GEOMPACK - a software package for the generation of meshes
//    using geometric algorithms,
//    Advances in Engineering Software,
//    Volume 13, pages 325-331, 1991.
//
//  Modified:
//
//    02 September 2003
//
//  Parameters:
//
//    Input, double X, Y, the coordinates of a point outside the convex hull
//    of the current triangulation.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//
//    Input, int TRI_NUM, the number of triangles.
//
//    Input, int TRI_VERT[TRI_NUM*3], the triangle incidence list.
//
//    Input, int TRI_NABE[TRI_NUM*3], the triangle neighbor list; negative
//    values are used for links of a counter clockwise linked list of boundary
//    edges;
//      LINK = -(3*I + J-1) where I, J = triangle, edge index.
//
//    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
//    assumed to be already computed and are not changed, else they are updated.
//    On output, LTRI is the index of boundary triangle to the left of the
//    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
//    edge of triangle LTRI to the left of the leftmost boundary edge visible
//    from (X,Y).  1 <= LEDG <= 3.
//
//    Input/output, int *RTRI.  On input, the index of the boundary triangle
//    to begin the search at.  On output, the index of the rightmost boundary
//    triangle visible from (X,Y).
//
//    Input/output, int *REDG, the edge of triangle RTRI that is visible
//    from (X,Y).  1 <= REDG <= 3.
//
{
  int a;
  double ax;
  double ay;
  int b;
  double bx;
  double by;
  bool done;
  int e;
  int l;
  int lr;
  int t;
//
//  Find the rightmost visible boundary edge using links, then possibly
//  leftmost visible boundary edge using triangle neighbor information.
//
  if ( *ltri == 0 )
  {
    done = false;
    *ltri = *rtri;
    *ledg = *redg;
  }
  else
  {
    done = true;
  }

  for ( ; ; )
  {
    l = -tri_nabe[3*((*rtri)-1)+(*redg)-1];
    t = l / 3;
    e = 1 + l % 3;
    a = tri_vert[3*(t-1)+e-1];

    if ( e <= 2 )
    {
      b = tri_vert[3*(t-1)+e];
    }
    else
    {
      b = tri_vert[3*(t-1)+0];
    }

    ax = point_xy[2*(a-1)+0];
    ay = point_xy[2*(a-1)+1];

    bx = point_xy[2*(b-1)+0];
    by = point_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

    *rtri = t;
    *redg = e;

  }

  if ( done )
  {
    return;
  }

  t = *ltri;
  e = *ledg;

  for ( ; ; )
  {
    b = tri_vert[3*(t-1)+e-1];
    e = i4_wrap ( e-1, 1, 3 );

    while ( 0 < tri_nabe[3*(t-1)+e-1] )
    {
      t = tri_nabe[3*(t-1)+e-1];

      if ( tri_vert[3*(t-1)+0] == b )
      {
        e = 3;
      }
      else if ( tri_vert[3*(t-1)+1] == b )
      {
        e = 1;
      }
      else
      {
        e = 2;
      }

    }

    a = tri_vert[3*(t-1)+e-1];
    ax = point_xy[2*(a-1)+0];
    ay = point_xy[2*(a-1)+1];

    bx = point_xy[2*(b-1)+0];
    by = point_xy[2*(b-1)+1];

    lr = lrline ( x, y, ax, ay, bx, by, 0.0 );

    if ( lr <= 0 )
    {
      break;
    }

  }

  *ltri = t;
  *ledg = e;

  return;
}
