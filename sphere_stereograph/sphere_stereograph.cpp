# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "sphere_stereograph.hpp"

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

void plane_normal_basis_3d ( double pp[3], double pn[3], double pq[3],
  double pr[3] )

//****************************************************************************80
//
//  Purpose:
//
//    PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
//
//  Discussion:
//
//    The normal form of a plane in 3D is:
//
//      PP is a point on the plane,
//      N is a normal vector to the plane.
//
//    The two vectors to be computed, PQ and PR, can be regarded as
//    the basis of a Cartesian coordinate system for points in the plane.
//    Any point in the plane can be described in terms of the "origin"
//    point PP plus a weighted sum of the two vectors PQ and PR:
//
//      P = PP + a * PQ + b * PR.
//
//    The vectors PQ and PR have unit length, and are perpendicular to N
//    and to each other.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double PP[3], a point on the plane.
//
//    Input, double PN[3], a normal vector to the plane.  The
//    vector must not have zero length, but it is not necessary for PN
//    to have unit length.
//
//    Output, double PQ[3], a vector of unit length, perpendicular
//    to the vector PN and the vector PR.
//
//    Output, double PR[3], a vector of unit length, perpendicular
//    to the vector PN and the vector PQ.
//
{
# define DIM_NUM 3

  int i;
  double normal_norm;
  double pr_norm;
  double *temp;
//
//  Compute the length of NORMAL.
//
  normal_norm = r8vec_norm ( DIM_NUM, pn );

  if ( normal_norm == 0.0 )
  {
    cerr << "\n";
    cerr << "PLANE_NORMAL_BASIS_3D - Fatal error!\n";
    cerr << "  The normal vector is 0.\n";
    exit ( 1 );
  }
//
//  Find a vector PQ that is normal to PN and has unit length.
//
  temp = r8vec_any_normal ( DIM_NUM, pn );
  r8vec_copy ( DIM_NUM, temp, pq );
  delete [] temp;
//
//  Now just take the cross product PN x PQ to get the PR vector.
//
  temp = r8vec_cross_product_3d ( pn, pq );

  pr_norm = r8vec_norm ( DIM_NUM, temp );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    pr[i] = temp[i] / pr_norm;
  }
  delete [] temp;

  return;
# undef DIM_NUM
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

  if ( m <= 0 || n <= 0 )
  {
    cout << "\n";
    cout << "  (None)\n";
    return;
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

double *r8mat_uniform_01_new ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01_NEW returns a unit pseudorandom R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8's,  stored as a vector
//    in column-major order.
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
//    03 October 2005
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
//    Philip Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0, otherwise the output value of SEED
//    will still be 0, and R8_UNIFORM will be 0.  On output, SEED has
//    been updated.
//
//    Output, double R8MAT_UNIFORM_01_NEW[M*N], a matrix of pseudorandom values.
//
{
  int i;
  int j;
  int k;
  double *r;

  r = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      k = *seed / 127773;

      *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

      if ( *seed < 0 )
      {
        *seed = *seed + 2147483647;
      }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number//
//
      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double *r8vec_any_normal ( int dim_num, double v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ANY_NORMAL returns some normal vector to V1.
//
//  Discussion:
//
//    If DIM_NUM < 2, then no normal vector can be returned.
//
//    If V1 is the zero vector, then any unit vector will do.
//
//    No doubt, there are better, more robust algorithms.  But I will take
//    just about ANY reasonable unit vector that is normal to V1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double V1[DIM_NUM], the vector.
//
//    Output, double R8VEC_ANY_NORMAL[DIM_NUM], a vector that is
//    normal to V2, and has unit Euclidean length.
//
{
  int i;
  int j;
  int k;
  double *v2;
  double vj;
  double vk;

  if ( dim_num < 2 )
  {
    cerr << "\n";
    cerr << "R8VEC_ANY_NORMAL - Fatal error!\n";
    cerr << "  Called with DIM_NUM < 2.\n";
    exit ( 1 );
  }

  v2 = new double[dim_num];

  if ( r8vec_norm ( dim_num, v1 ) == 0.0 )
  {
    r8vec_zero ( dim_num, v2 );
    v2[0] = 1.0;
    return v2;
  }
//
//  Seek the largest entry in V1, VJ = V1(J), and the
//  second largest, VK = V1(K).
//
//  Since V1 does not have zero norm, we are guaranteed that
//  VJ, at least, is not zero.
//
  j = -1;
  vj = 0.0;

  k = -1;
  vk = 0.0;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( r8_abs ( vk ) < r8_abs ( v1[i] ) || k == -1 )
    {
      if ( r8_abs ( vj ) < r8_abs ( v1[i] ) || j == -1 )
      {
        k = j;
        vk = vj;
        j = i;
        vj = v1[i];
      }
      else
      {
        k = i;
        vk = v1[i];
      }
    }
  }
//
//  Setting V2 to zero, except that V2(J) = -VK, and V2(K) = VJ,
//  will just about do the trick.
//
  r8vec_zero ( dim_num, v2 );

  v2[j] = -vk / sqrt ( vk * vk + vj * vj );
  v2[k] =  vj / sqrt ( vk * vk + vj * vj );

  return v2;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
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
//    Input, double A2[N], the copy of A1.
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

double *r8vec_cross_product_3d ( double v1[3], double v2[3] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_CROSS_PRODUCT_3D computes the cross product of two R8VEC's in 3D.
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
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_PRODUCT_3D[3], the cross product vector.
//
{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//****************************************************************************80

double r8vec_norm ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM returns the L2 norm of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The vector L2 norm is defined as:
//
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)^2 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in A.
//
//    Input, double A[N], the vector whose L2 norm is desired.
//
//    Output, double R8VEC_NORM, the L2 norm of A.
//
{
  int i;
  double v;

  v = 0.0;

  for ( i = 0; i < n; i++ )
  {
    v = v + a[i] * a[i];
  }
  v = sqrt ( v );

  return v;
}
//****************************************************************************80

double r8vec_norm_affine ( int n, double v0[], double v1[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORM_AFFINE returns the affine L2 norm of an R8VEC.
//
//  Discussion:
//
//    The affine vector L2 norm is defined as:
//
//      R8VEC_NORM_AFFINE(V0,V1)
//        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vectors.
//
//    Input, double V0[N], the base vector.
//
//    Input, double V1[N], the vector whose affine L2 norm is desired.
//
//    Output, double R8VEC_NORM_AFFINE, the affine L2 norm of V1.
//
{
  int i;
  double value;

  value = 0.0;

  for ( i = 0; i < n; i++ )
  {
    value = value + ( v1[i] - v0[i] ) * ( v1[i] - v0[i] );
  }
  value = sqrt ( value );

  return value;
}
//*****************************************************************************

void r8vec_normal_01 ( int n, int *seed, double x[] )

//*****************************************************************************
//
//  Purpose:
//
//    R8VEC_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
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
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double X[N], a sample of the standard normal PDF.
//
{
  int i;
  int m;
  double pi = 3.141592653589793;
  double *r;
  static int made = 0;
  static int saved = 0;
  int xhi;
  int xlo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the random number seed.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return;
  }
  else if ( n == 0 )
  {
    return;
  }
//
//  Record the range of X we need to fill in.
//
  xlo = 1;
  xhi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    xlo = 2;
  }
//
//  If we don't need any more values, return.
//
  if ( xhi - xlo + 1 == 0 )
  {
    return;
  }

  r = new double[n+1];
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( xhi - xlo + 1 == 1 )
  {
    r8vec_uniform_01 ( 2, seed, r );

    x[xhi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =        sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( ( xhi-xlo+1) % 2 ) == 0 )
  {
    m = ( xhi-xlo+1 ) / 2;

    r8vec_uniform_01 ( 2*m, seed, r );

    for ( i = 0; i < 2*m; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i]   = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    made = made + xhi - xlo + 1;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    xhi = xhi - 1;

    m = ( xhi-xlo+1 ) / 2 + 1;

    r8vec_uniform_01 ( 2*m, seed, r );

    for ( i = 0; i < 2*m-2; i = i + 2 )
    {
      x[xlo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[xlo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    x[n-1] = sqrt ( -2.0 * log ( r[2*m-2] ) ) * cos ( 2.0 * pi * r[2*m-1] );
    y =      sqrt ( -2.0 * log ( r[2*m-2] ) ) * sin ( 2.0 * pi * r[2*m-1] );

    saved = 1;

    made = made + xhi - xlo + 2;

  }

  delete [] r;

  return;
}
//****************************************************************************80

void r8vec_transpose_print ( int n, double x[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT prints a vector on one line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, double X[N], the vector.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << title;
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << x[i];
  }
  cout << "\n";

  return;
}
//****************************************************************************80

void r8vec_uniform_01 ( int n, int *seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 fills a double precision vector with pseudorandom values.
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
//    19 August 2004
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
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int k;

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

  return;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
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
  int i4_huge = 2147483647;
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
      *seed = *seed + i4_huge;
    }

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void r8vec_zero ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO zeroes an R8VEC.
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
//    Input, int N, the number of entries in the vector.
//
//    Output, double A[N], a vector of zeroes.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return;
}
//****************************************************************************80

double *sphere_stereograph ( int m, int n, double p[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_STEREOGRAPH computes the stereographic image of points on a sphere.
//
//  Discussion:
//
//    We start with a sphere of radius 1 and center (0,0,0).
//
//    The north pole N = (0,0,1) is the point of tangency to the sphere
//    of a plane, and the south pole S = (0,0,-1) is the focus for the
//    stereographic projection.
//
//    For any point P on the sphere, the stereographic projection Q of the
//    point is defined by drawing the line from S through P, and computing
//    Q as the intersection of this line with the plane.
//
//    Actually, we allow the spatial dimension M to be arbitrary.  Values
//    of M make sense starting with 2.  The north and south poles are
//    selected as the points (0,0,...,+1) and (0,0,...,-1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    C F Marcus,
//    The stereographic projection in vector notation,
//    Mathematics Magazine,
//    Volume 39, Number 2, March 1966, pages 100-102.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double P[M*N], a set of points on the unit sphere.
//
//    Output, double SPHERE_STEREOGRAPH[M*N], the coordinates of the
//    image points.
//
{
  int i;
  int j;
  double *q;

  q = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      q[i+j*m] = 2.0 * p[i+j*m] / ( 1.0 + p[m-1+j*m] );
    }
    q[m-1+j*m] = 1.0;
  }

  return q;
}
//****************************************************************************80

double *sphere_stereograph_inverse ( int m, int n, double q[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_STEREOGRAPH_INVERSE computes stereographic preimages of points.
//
//  Discussion:
//
//    We start with a sphere of radius 1 and center (0,0,0).
//
//    The north pole N = (0,0,1) is the point of tangency to the sphere
//    of a plane, and the south pole S = (0,0,-1) is the focus for the
//    stereographic projection.
//
//    For any point Q on the plane, the stereographic inverse projection
//    P of the point is defined by drawing the line from S through Q, and
//    computing P as the intersection of this line with the sphere.
//
//    Actually, we allow the spatial dimension M to be arbitrary.  Values
//    of M make sense starting with 2.  The north and south poles are
//    selected as the points (0,0,...,+1) and (0,0,...,-1).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    C F Marcus,
//    The stereographic projection in vector notation,
//    Mathematics Magazine,
//    Volume 39, Number 2, March 1966, pages 100-102.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Q[M*N], the points, which are presumed to lie
//    on the plane Z = 1.
//
//    Output, double SPHERE_STEREOGRAPH_INVERSE[M*N], the stereographic
//    inverse projections of the points.
//
{
  int i;
  int j;
  double *p;
  double qn;

  p = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    qn = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      qn = qn + pow ( q[i+j*m], 2 );
    }
    for ( i = 0; i < m - 1; i++ )
    {
      p[i+j*m] = 4.0 * q[i+j*m] / ( 4.0 + qn );
    }
    p[m-1+j*m] = ( 4.0 - qn ) / ( 4.0 + qn );
  }

  return p;
}
//****************************************************************************80

double *sphere_stereograph2 ( int m, int n, double p[], double focus[],
  double center[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_STEREOGRAPH2 computes the stereographic image of points on a sphere.
//
//  Discussion:
//
//    We start with a sphere of center C.
//
//    F is a point on the sphere which is the focus of the mapping,
//    and the antipodal point 2*C-F is the point of tangency
//    to the sphere of a plane.
//
//    For any point P on the sphere, the stereographic projection Q of the
//    point is defined by drawing the line from F through P, and computing
//    Q as the intersection of this line with the plane.
//
//    The spatial dimension M is arbitrary, but should be at least 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    C F Marcus,
//    The stereographic projection in vector notation,
//    Mathematics Magazine,
//    Volume 39, Number 2, March 1966, pages 100-102.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double P[M*N], a set of points on the unit sphere.
//
//    Input, double FOCUS[M], the coordinates of the focus point.
//
//    Input, double CENTER[M], the coordinates of the center of the sphere.
//
//    Output, double SPHERE_STEREOGRAPH2[M*N], the coordinates of the
//    image points,
//
{
  double cf_dot_pf;
  double cf_normsq;
  int i;
  int j;
  double *q;
  double s;

  q = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    cf_normsq = 0.0;
    cf_dot_pf = 0.0;
    for ( i = 0; i < m; i++ )
    {
      cf_normsq = cf_normsq + pow ( center[i] - focus[i], 2 );
      cf_dot_pf = cf_dot_pf + ( center[i] - focus[i] ) * ( p[i+j*m] - focus[i] );
    }
    s = 2.0 * cf_normsq / cf_dot_pf;
    for ( i = 0; i < m; i++ )
    {
      q[i+j*m] = s * p[i+j*m] + ( 1.0 - s ) * focus[i];
    }
  }
  return q;
}
//****************************************************************************80

double *sphere_stereograph2_inverse ( int m, int n, double q[], double focus[],
  double center[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_STEREOGRAPH2_INVERSE computes stereographic preimages of points.
//
//  Discussion:
//
//    We start with a sphere of center C.
//
//    F is a point on the sphere which is the focus of the mapping,
//    and the antipodal point 2*C-F is the point of tangency
//    to the sphere of a plane.
//
//    For any point Q on the plane, the stereographic inverse projection
//    P of the point is defined by drawing the line from F through Q, and
//    computing P as the intersection of this line with the sphere.
//
//    The spatial dimension M is arbitrary, but should be at least 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    C F Marcus,
//    The stereographic projection in vector notation,
//    Mathematics Magazine,
//    Volume 39, Number 2, March 1966, pages 100-102.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Q[M*N], the points, which are presumed to lie
//    on the plane.
//
//    Input, double FOCUS[M], the coordinates of the focus point.
//
//    Input, double CENTER[M], the coordinates of the center of the sphere.
//
//    Output, double SPHERE_STEREOGRAPH2_INVERSE[M*N], the stereographic
//    inverse projections of the points.
//
{
  double cf_dot_qf;
  int i;
  int j;
  double *p;
  double qf_normsq;
  double s;

  p = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    cf_dot_qf = 0.0;
    qf_normsq = 0.0;
    for ( i = 0; i < m; i++ )
    {
      cf_dot_qf = cf_dot_qf + ( center[i] - focus[i] ) * ( q[i+j*m] - focus[i] );
      qf_normsq = qf_normsq + pow ( q[i+j*m] - focus[i], 2 );
    }
    s = 2.0 * cf_dot_qf / qf_normsq;
    for ( i = 0; i < m; i++ )
    {
      p[i+j*m] = s * q[i+j*m] + ( 1.0 - s ) * focus[i];
    }
  }
  return p;
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

double *uniform_on_sphere01_map ( int dim_num, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_MAP maps uniform points onto the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
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
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 168.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Krieger, 1992,
//    ISBN: 0894647644,
//    LC: QA298.R79.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_MAP[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  double norm;
  double *u;
  double *x;

  u = new double[dim_num];
  x = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    r8vec_normal_01 ( dim_num, seed, u );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      norm = norm + u[i] * u[i];
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = u[i] / norm;
    }

  }

  delete [] u;
  return x;
}
