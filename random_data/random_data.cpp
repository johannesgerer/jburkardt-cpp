# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "random_data.hpp"

//****************************************************************************80

double arc_cosine ( double c )

//****************************************************************************80
//
//  Purpose:
//
//    ARC_COSINE computes the arc cosine function, with argument truncation.
//
//  Discussion:
//
//    If you call your system ACOS routine with an input argument that is
//    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
//    This routine truncates arguments outside the range.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double C, the argument, the cosine of an angle.
//
//    Output, double ARC_COSINE, an angle whose cosine is C.
//
{
  double angle;
  double pi = 3.141592653589793;

  if ( c <= -1.0 )
  {
    angle = pi;
  }
  else if ( 1.0 <= c )
  {
    angle = 0.0;
  }
  else
  {
    angle = acos ( c );
  }
  return angle;
}
//****************************************************************************80

double *bad_in_simplex01 ( int dim_num, int point_num, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    BAD_IN_SIMPLEX01 is a "bad" (nonuniform) sampling of the unit simplex.
//
//  Discussion:
//
//    The interior of the unit DIM_NUM-dimensional simplex is the set of
//    points X(1:DIM_NUM) such that each X(I) is nonnegative, and
//    sum(X(1:DIM_NUM)) <= 1.
//
//    Any point in the unit simplex CAN be chosen by this algorithm.
//
//    However, the points that are chosen tend to be clustered near
//    the centroid.
//
//    This routine is supplied as an example of "bad" sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double BAD_IN_SIMPLEX01[DIM_NUM*POINT_NUM], the points.
//
{
  double *e;
  double e_sum;
  int i;
  int j;
  double *x;

  x = new double[dim_num*point_num];

  for ( j = 0; j < point_num; j++ )
  {
    e = r8vec_uniform_01_new ( dim_num + 1, seed );

    e_sum = r8vec_sum ( dim_num + 1, e );

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = e[i] / e_sum;
    }

    delete [] e;
  }

  return x;
}
//****************************************************************************80

double *brownian ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    BROWNIAN creates Brownian motion points.
//
//  Discussion:
//
//    A starting point is generated at the origin.  The next point
//    is generated at a uniformly random angle and a (0,1) normally
//    distributed distance from the previous point.
//
//    It is up to the user to rescale the data, if desired.
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
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, int &SEED, a seed for the random number generator.
//
//    Output, double BROWNIAN[DIM_NUM*N], the Brownian motion points.
//
{
  double *direction;
  int i;
  int j;
  double r;
  double *x;

  direction = new double[dim_num];
  x = new double[dim_num*n];
//
//  Initial point.
//
  j = 0;
  for ( i = 0; i < dim_num; i++ )
  {
    x[i+j*dim_num] = 0.0;
  }
//
//  Generate angles and steps.
//
  for ( j = 1; j < n; j++ )
  {
    r = r8_normal_01 ( seed );
    r = fabs ( r );

    direction_uniform_nd ( dim_num, seed, direction );

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = x[i+(j-1)*dim_num] + r * direction[i];
    }

  }

  delete [] direction;

  return x;
}
//****************************************************************************80

void daxpy ( int n, double da, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DAXPY computes constant times a vector plus a vector.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Lawson, Hanson, Kincaid, Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of elements in DX and DY.
//
//    Input, double DA, the multiplier of DX.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries of DX.
//
//    Input/output, double DY[*], the second vector.
//    On output, DY[*] has been replaced by DY[*] + DA * DX[*].
//
//    Input, int INCY, the increment between successive entries of DY.
//
{
  int i;
  int ix;
  int iy;
  int m;

  if ( n <= 0 )
  {
    return;
  }

  if ( da == 0.0 )
  {
    return;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dy[iy] = dy[iy] + da * dx[ix];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 4;

    for ( i = 0; i < m; i++ )
    {
      dy[i] = dy[i] + da * dx[i];
    }

    for ( i = m; i < n; i = i + 4 )
    {
      dy[i  ] = dy[i  ] + da * dx[i  ];
      dy[i+1] = dy[i+1] + da * dx[i+1];
      dy[i+2] = dy[i+2] + da * dx[i+2];
      dy[i+3] = dy[i+3] + da * dx[i+3];
    }

  }

  return;
}
//****************************************************************************80

double ddot ( int n, double dx[], int incx, double dy[], int incy )

//****************************************************************************80
//
//  Purpose:
//
//    DDOT forms the dot product of two vectors.
//
//  Discussion:
//
//    This routine uses unrolled loops for increments equal to one.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Jack Dongarra.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979.
//
//    Lawson, Hanson, Kincaid, Krogh,
//    Basic Linear Algebra Subprograms for Fortran Usage,
//    Algorithm 539,
//    ACM Transactions on Mathematical Software,
//    Volume 5, Number 3, September 1979, pages 308-323.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double DX[*], the first vector.
//
//    Input, int INCX, the increment between successive entries in DX.
//
//    Input, double DY[*], the second vector.
//
//    Input, int INCY, the increment between successive entries in DY.
//
//    Output, double DDOT, the sum of the product of the corresponding
//    entries of DX and DY.
//
{
  double dtemp;
  int i;
  int ix;
  int iy;
  int m;

  dtemp = 0.0;

  if ( n <= 0 )
  {
    return dtemp;
  }
//
//  Code for unequal increments or equal increments
//  not equal to 1.
//
  if ( incx != 1 || incy != 1 )
  {
    if ( 0 <= incx )
    {
      ix = 0;
    }
    else
    {
      ix = ( - n + 1 ) * incx;
    }

    if ( 0 <= incy )
    {
      iy = 0;
    }
    else
    {
      iy = ( - n + 1 ) * incy;
    }

    for ( i = 0; i < n; i++ )
    {
      dtemp = dtemp + dx[ix] * dy[iy];
      ix = ix + incx;
      iy = iy + incy;
    }
  }
//
//  Code for both increments equal to 1.
//
  else
  {
    m = n % 5;

    for ( i = 0; i < m; i++ )
    {
      dtemp = dtemp + dx[i] * dy[i];
    }

    for ( i = m; i < n; i = i + 5 )
    {
      dtemp = dtemp + dx[i  ] * dy[i  ]
                    + dx[i+1] * dy[i+1]
                    + dx[i+2] * dy[i+2]
                    + dx[i+3] * dy[i+3]
                    + dx[i+4] * dy[i+4];
    }

  }

  return dtemp;
}
//****************************************************************************80

double *dge_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGE_MXV multiplies a DGE matrix times a vector.
//
//  Discussion:
//
//    The DGE storage format is used for a general M by N matrix.  A physical storage
//    space is made for each logical entry.  The two dimensional logical
//    array is mapped to a vector, in which storage is by columns.
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
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the SGE matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DGE_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

void direction_uniform_nd ( int dim_num, int &seed, double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIRECTION_UNIFORM_ND generates a random direction vector in ND.
//
//  Discussion:
//
//    This is actually simply a random point on the unit sphere in ND.
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
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double W[DIM_NUM], a random direction vector, with unit norm.
//
{
  int i;
  double norm;
//
//  Sample the standard normal distribution.
//
  r8vec_normal_01 ( dim_num, seed, w );
//
//  Compute the length of the vector.
//
  norm = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    norm = norm + w[i] * w[i];
  }
  norm = sqrt ( norm );
//
//  Normalize the vector.
//
  for ( i = 0; i < dim_num; i++ )
  {
    w[i] = w[i] / norm;
  }

  return;
}
//****************************************************************************80

double *dut_mxv ( int m, int n, double a[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    DUT_MXV multiplies an DUT matrix times a vector.
//
//  Discussion:
//
//    The DUT storage format is used for an M by N upper triangular matrix,
//    and allocates space even for the zero entries.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 September 2003
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
//    Input, double A[M*N], the DUT matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double DUT_MXV[M], the product A * x.
//
{
  double *b;
  int i;
  int j;

  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    b[i] = 0.0;
    for ( j = i; j < n; j++ )
    {
      b[i] = b[i] + a[i+j*m] * x[j];
    }
  }

  return b;
}
//****************************************************************************80

unsigned long get_seed ( )

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
//    17 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, unsigned long GET_SEED, a random seed value.
//
{
# define UNSIGNED_LONG_MAX 4294967295UL
  time_t clock;
  int i;
  int hours;
  int minutes;
  int seconds;
  struct tm *lt;
  static unsigned long seed = 0;
  time_t tloc;
//
//  If the internal seed is 0, generate a value based on the time.
//
  if ( seed == 0 )
  {
    clock = time ( &tloc );
    lt = localtime ( &clock );
//
//  Extract HOURS.
//
    hours = lt->tm_hour;
//
//  In case of 24 hour clocks, shift so that HOURS is between 1 and 12.
//
    if ( 12 < hours )
    {
      hours = hours - 12;
    }
//
//  Move HOURS to 0, 1, ..., 11
//
    hours = hours - 1;

    minutes = lt->tm_min;

    seconds = lt->tm_sec;

    seed = seconds + 60 * ( minutes + 60 * hours );
//
//  We want values in [1,43200], not [0,43199].
//
    seed = seed + 1;
//
//  Remap SEED from [1,43200] to [1,UNSIGNED_LONG_MAX].
//
    seed = ( unsigned long )
      ( ( ( double ) seed )
      * ( ( double ) UNSIGNED_LONG_MAX ) / ( 60.0 * 60.0 * 12.0 ) );
  }
//
//  Never use a seed of 0.
//
  if ( seed == 0 )
  {
    seed = 1;
  }

  return seed;

# undef UNSIGNED_LONG_MAX
}
//****************************************************************************80

double *grid_in_cube01 ( int dim_num, int n, int center, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_IN_CUBE01 generates a grid dataset in the unit hypercube.
//
//  Discussion:
//
//    N points are needed in a DIM_NUM dimensional space.
//
//    The points are to lie on a uniform grid of side N_SIDE.
//
//    Unless the N = N_SIDE^DIM_NUM for some N_SIDE, we can't use all the
//    points on a grid.  What we do is find the smallest N_SIDE
//    that's big enough, and randomly omit some points.
//
//    If N_SIDE is 4, then the choices in 1D are:
//
//    A: 0,   1/3, 2/3, 1
//    B: 1/5, 2/5, 3/5, 4/5
//    C: 0,   1/4, 2/4, 3/4
//    D: 1/4, 2/4, 3/4, 1
//    E: 1/8, 3/8, 5/8, 7/8
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int CENTER, specifies the 1D grid centering:
//    1: first point is 0.0, last point is 1.0;
//    2: first point is 1/(N+1), last point is N/(N+1);
//    3: first point is 0, last point is (N-1)/N;
//    4: first point is 1/N, last point is 1;
//    5: first point is 1/(2*N), last point is (2*N-1)/(2*N);
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double GRID_IN_CUBE01[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  int n_grid;
  int n_side;
  double *r;
  int rank;
  int *rank_list;
  int *tuple;
//
//  Find the dimension of the smallest grid with N points.
//
  n_side = grid_side ( dim_num, n );
//
//  We need to select N points out of N_SIDE^DIM_NUM set.
//
  n_grid = ( int ) pow ( ( double ) n_side, dim_num );
//
//  Generate a random subset of N items from a set of size N_GRID.
//
  rank_list = new int[n];

  ksub_random2 ( n_grid, n, seed, rank_list );
//
//  Must make one dummy call to TUPLE_NEXT_FAST with RANK = -1.
//
  rank = -1;
  tuple = new int[dim_num];
  tuple_next_fast ( n_side, dim_num, rank, tuple );
//
//  Now generate the appropriate indices, and "center" them.
//
  r = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
    rank = rank_list[j] - 1;

    tuple_next_fast ( n_side, dim_num, rank, tuple );

    if ( center == 1 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] - 1 )
          / ( double ) ( n_side - 1 );
      }
    }
    else if ( center == 2 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] )
          / ( double ) ( n_side + 1 );
      }
    }
    else if ( center == 3 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] - 1 )
          / ( double ) ( n_side );
      }
    }
    else if ( center == 4 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( tuple[i] )
          / ( double ) ( n_side );
      }
    }
    else if ( center == 5 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        r[i+j*dim_num] = ( double ) ( 2 * tuple[i] - 1 )
          / ( double ) ( 2 * n_side );
      }
    }
  }

  delete [] rank_list;
  delete [] tuple;

  return r;
}
//****************************************************************************80

int grid_side ( int dim_num, int n )

//****************************************************************************80
//
//  Purpose:
//
//    GRID_SIDE finds the smallest DIM_NUM dimensional grid containing at least N points.
//
//  Discussion:
//
//    Each coordinate of the grid will have N_SIDE distinct values.
//    Thus the total number of points in the grid is N_SIDE**DIM_NUM.
//    This routine seeks the smallest N_SIDE such that N <= N_SIDE**DIM_NUM.
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
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Output, int GRID_SIDE, the length of one side of the smallest
//    grid in DIM_NUM dimensions that contains at least N points.
//
{
  double exponent;
  int n_grid;
  int n_side;

  if ( n <= 0 )
  {
    n_side = 0;
    return n_side;
  }

  if ( dim_num <= 0 )
  {
    n_side = -1;
    return n_side;
  }

  exponent = 1.0 / ( double ) dim_num;

  n_side = ( int ) pow ( n, exponent );

  if ( pow ( ( double ) n_side, ( double ) dim_num ) < n )
  {
    n_side = n_side + 1;
  }

  return n_side;
}
//****************************************************************************80

bool halham_dim_num_check ( int dim_num )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_DIM_NUM_CHECK checks DIM_NUM for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//    DIM_NUM must be positive.
//
//    Output, bool HALHAM_DIM_NUM_CHECK, is true if DIM_NUM is legal.
//
{
  bool value;

  if ( dim_num < 1 )
  {
    cerr << "\n";
    cerr << "HALHAM_DIM_NUM_CHECK - Fatal error!\n";
    cerr << "  DIM_NUM < 0.";
    cerr << "  DIM_NUM = " << dim_num << "\n";
    exit ( 1 );
  }

  value = true;

  return value;
}
//****************************************************************************80

bool halham_leap_check ( int dim_num, int leap[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_LEAP_CHECK checks LEAP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int LEAP[DIM_NUM], the successive jumps in the sequence.
//    Each entry must be greater than 0.
//
//    Output, bool HALHAM_LEAP_CHECK, is true if LEAP is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( leap[i] < 1 )
    {
      cerr << "\n";
      cerr << "HALHAM_LEAP_CHECK - Fatal error!\n";
      cerr << "  Leap entries must be greater than 0.\n";
      cerr << "  leap[" << i << "] = " << leap[i] << "\n";
      exit ( 1 );
    }
  }

  return value;
}
//****************************************************************************80

bool halham_n_check ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_N_CHECK checks N for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of elements of the subsequence.
//    N must be positive.
//
//    Output, bool HALHAM_N_CHECK, is true if N is legal.
//
{
  bool value;

  if ( n < 1 )
  {
    cerr << "\n";
    cerr << "HALHAM_N_CHECK - Fatal error!\n";
    cerr << "  N < 0.";
    cerr << "  N = " << n << "\n";
    exit ( 1 );
  }

  value = true;

  return value;
}
//****************************************************************************80

bool halham_seed_check ( int dim_num, int seed[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_SEED_CHECK checks SEED for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int SEED[DIM_NUM], the sequence index
//    corresponding to STEP = 0.  Each entry must be 0 or greater.
//
//    Output, bool HALHAM_SEED_CHECK, is true if SEED is legal.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( seed[i] < 0 )
    {
      cerr << "\n";
      cerr << "HALHAM_SEED_CHECK - Fatal error!\n";
      cerr << "  SEED entries must be nonnegative.\n";
      cerr << "  seed[" << i << "] = " << seed[i] << "\n";
      exit ( 1 );
    }
  }

  return value;
}
//****************************************************************************80

bool halham_step_check ( int step )

//****************************************************************************80
//
//  Purpose:
//
//    HALHAM_STEP_CHECK checks STEP for a Halton or Hammersley sequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int STEP, the index of the subsequence element.
//    STEP must be 1 or greater.
//
//    Output, bool HALHAM_STEP_CHECK, is true if STEP is legal.
//
{
  int i;
  bool value;

  if ( step < 0 )
  {
    cerr << "\n";
    cerr << "HALHAM_STEP_CHECK - Fatal error!\n";
    cerr << "  STEP < 0.";
    cerr << "  STEP = " << step << "\n";
    exit ( 1 );
  }

  value = true;

  return value;
}
//****************************************************************************80

bool halton_base_check ( int dim_num, int base[] )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_BASE_CHECK is TRUE if BASE is legal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the Halton bases.
//    Each base must be greater than 1.
//
//    Output, bool HALTON_BASE_CHECK.
//
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( base[i] <= 1 )
    {
      cerr << "\n";
      cerr << "HALTON_BASE_CHECK - Fatal error!\n";
      cerr << "  Bases must be greater than 1.\n";
      cerr << "  base[" << i << "] = " << base[i] << "\n";
      exit ( 1 );
    }
  }

  return value;
}
//****************************************************************************80

double *halton_in_circle01_accept ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_IN_CIRCLE01_ACCEPT accepts Halton points in a unit circle.
//
//  Discussion:
//
//    The acceptance/rejection method is used.
//
//    The unit circle is centered at the origin and has radius 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CIRCLE01_ACCEPT[DIM_NUM*N], the points.
//
{
  int *base;
  int have;
  int i;
  int j;
  int *leap;
  int *seed_vec;
  int step;
  double total;
  double *u;
  double *x;

  base = new int[dim_num];
  leap = new int[dim_num];
  seed_vec = new int[dim_num];
  u = new double[dim_num];
  x = new double[dim_num*n];

  have = 0;

  for ( i = 0; i < dim_num; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < dim_num; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < dim_num; i++ )
  {
    base[i] = prime ( i + 1 );
  }

  while ( have < n )
  {
    step = seed;

    i4_to_halton ( dim_num, step, seed_vec, leap, base, u );

    seed = seed + 1;

    total = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      u[i] = 2.0 * u[i] - 1.0;
      total = total + u[i] * u[i];
    }

    if ( total <= 1.0 )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        x[i+have*dim_num] = u[i];
      }
      have = have + 1;
    }
  }

  delete [] base;
  delete [] leap;
  delete [] seed_vec;
  delete [] u;

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *halton_in_circle01_map ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_IN_CIRCLE01_MAP maps Halton points into a unit circle.
//
//  Discussion:
//
//    The unit circle is centered at the origin and has radius 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CIRCLE01_MAP[DIM_NUM*N], the points.
//
{
# define PI 3.141592653589793

  int base[1];
  int j;
  int leap[1];
  double *r;
  double rval;
  int step;
  int seed_vec[1];
  double *t;
  double *x;

  r = new double[n];
  t = new double[n];
  x = new double[dim_num*n];

  step = 0;
  seed_vec[0] = seed;
  leap[0] = 1;
  base[0] = prime ( 1 );

  i4_to_halton_sequence ( 1, n, step, seed_vec, leap, base, r );

  for ( j = 0; j < n; j++ )
  {
    r[j] = sqrt ( r[j] );
  }

  step = 0;
  seed_vec[0] = seed;
  leap[0] = 1;
  base[0] = prime ( 2 );

  i4_to_halton_sequence ( 1, n, step, seed_vec, leap, base, t );

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * t[j];
  }

  for ( j = 0; j < n; j++ )
  {
    x[0+j*dim_num] = r[j] * cos ( t[j] );
    x[1+j*dim_num] = r[j] * sin ( t[j] );
  }

  seed = seed + n;

  delete [] r;
  delete [] t;

  return x;
# undef PI
}
//****************************************************************************80

double *halton_in_cube01 ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HALTON_IN_CUBE01 generates Halton points in the unit hypercube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double HALTON_IN_CUBE01[DIM_NUM*N], the points
//
{
  int *base;
  int i;
  int *leap;
  int *seed_vec;
  int step;
  double *x;

  base = new int[dim_num];
  leap = new int[dim_num];
  seed_vec = new int[dim_num];
  x = new double[dim_num*n];

  step = seed;
  for ( i = 0; i < dim_num; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < dim_num; i++ )
  {
    leap[i] = 1;
  }
  for ( i = 0; i < dim_num; i++ )
  {
    base[i] = prime ( i + 1 );
  }

  i4_to_halton_sequence ( dim_num, n, step, seed_vec, leap, base, x );

  seed = seed + n;

  delete [] base;
  delete [] leap;
  delete [] seed_vec;

  return x;
}
//****************************************************************************80

bool hammersley_base_check ( int dim_num, int base[] )

//****************************************************************************80
//
//  Purpose:
//
//    HAMMERSLEY_BASE_CHECK is TRUE if BASE is legal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int BASE[DIM_NUM], the bases.
//
//    Output, bool HAMMERSLEY_BASE_CHECK.
{
  int i;
  bool value;

  value = true;

  for ( i = 0; i < dim_num; i++ )
  {
    if ( base[i] == 0 || base[i] == 1 )
    {
      cerr << "\n";
      cerr << "HAMMERSLEY_BASE_CHECK - Fatal error!\n";
      cerr << "  Bases may not be 0 or 1.\n";
      i4vec_transpose_print ( dim_num, base, "BASE:  " );
      exit ( 1 );
    }
  }

  return value;
}
//****************************************************************************80

double *hammersley_in_cube01 ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HAMMERSLEY_IN_CUBE01 computes Hammersley points in the unit hypercube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double HAMMERSLEY_IN_CUBE01[DIM_NUM*N], the points.
//
{
  int *base;
  int i;
  int *leap;
  int *seed_vec;
  int step;
  double *x;

  base = new int[dim_num];
  leap = new int[dim_num];
  seed_vec = new int[dim_num];
  x = new double[dim_num*n];

  step = seed;
  for ( i = 0; i < dim_num; i++ )
  {
    seed_vec[i] = 0;
  }
  for ( i = 0; i < dim_num; i++ )
  {
    leap[i] = 1;
  }
  base[0] = -n;
  for ( i = 1; i < dim_num; i++ )
  {
    base[i] = prime ( i );
  }

  i4_to_hammersley_sequence ( dim_num, n, step, seed_vec, leap, base, x );

  seed = seed + n;

  delete [] base;
  delete [] leap;
  delete [] seed_vec;

  return x;
}
//****************************************************************************80

int i4_factorial ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    I4_FACTORIAL returns N!.
//
//  Discussion:
//
//    N! = Product ( 1 <= I <= N ) I
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the argument of the factorial function.
//    0 <= N.
//
//    Output, int I4_FACTORIAL, the factorial of N.
//
{
  int fact;
  int i;
//
//  Check.
//
  if ( n < 0 )
  {
    cerr << "\n";
    cerr << "I4_FACTORIAL - Fatal error!\n";
    cerr << "  N < 0.\n";
    exit ( 1 );
  }

  fact = 1;

  for ( i = 2; i <= n; i++ )
  {
    fact = fact * i;
  }

  return fact;
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
//    I4_MIN returns the smaller of two I4's.
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
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
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

void i4_to_halton ( int dim_num, int step, int seed[], int leap[], int base[],
  double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON computes one element of a leaped Halton subsequence.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM is required.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required.
//
//    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
//    to STEP = 0.
//    0 <= SEED(1:DIM_NUM) is required.
//
//    Input, int LEAP[DIM_NUM], the successive jumps in the Halton sequence.
//    1 <= LEAP(1:DIM_NUM) is required.
//
//    Input, int BASE[DIM_NUM], the Halton bases.
//    1 < BASE(1:DIM_NUM) is required.
//
//    Output, double R[DIM_NUM], the STEP-th element of the leaped
//    Halton subsequence.
//
{
  double base_inv;
  int digit;
  int i;
  int seed2;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }

  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  for ( i = 0; i < dim_num; i++ )
  {
    seed2 = seed[i] + step * leap[i];

    r[i] = 0.0;

    base_inv = 1.0 / ( ( double ) base[i] );

    while ( seed2 != 0 )
    {
      digit = seed2 % base[i];
      r[i] = r[i] + ( ( double ) digit ) * base_inv;
      base_inv = base_inv / ( ( double ) base[i] );
      seed2 = seed2 / base[i];
    }
  }

  return;
}
//****************************************************************************80

void i4_to_halton_sequence ( int dim_num, int n, int step, int seed[],
  int leap[], int base[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON_SEQUENCE computes N elements of a leaped Halton subsequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Halton sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Halton sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
//
//    of the original Halton sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is DIM_NUM * N, where DIM_NUM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first DIM_NUM entries of R, followed by the DIM_NUM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*DIM_NUM through (DIM_NUM-1) + (J-1)*DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 July 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J H Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, 1960, pages 84-90.
//
//    J H Halton and G B Smith,
//    Algorithm 247: Radical-Inverse Quasi-Random Point Sequence,
//    Communications of the ACM,
//    Volume 7, 1964, pages 701-702.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements of the sequence.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[DIM_NUM], the Halton sequence index corresponding
//    to STEP = 0.
//
//    Input, int LEAP[DIM_NUM], the succesive jumps in the Halton sequence.
//
//    Input, int BASE[DIM_NUM], the Halton bases.
//
//    Output, double R[DIM_NUM*N], the next N elements of the
//    leaped Halton subsequence, beginning with element STEP.
//
{
  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }

  if ( !halham_n_check ( n ) )
  {
    exit ( 1 );
  }

  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, leap ) )
  {
    exit ( 1 );
  }

  if ( !halton_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < dim_num; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      seed2[j] = seed[i] + ( step + j ) * leap[i];
    }

    for ( j = 0; j < n; j++ )
    {
      r[i+j*dim_num] = 0.0;
    }

    for ( j = 0; j < n; j++ )
    {
      base_inv = 1.0 / ( ( double ) base[i] );

      while ( seed2[j] != 0 )
      {
        digit = seed2[j] % base[i];
        r[i+j*dim_num] = r[i+j*dim_num] + ( ( double ) digit ) * base_inv;
        base_inv = base_inv / ( ( double ) base[i] );
        seed2[j] = seed2[j] / base[i];
      }
    }
  }

  delete [] seed2;

  return;
}
//****************************************************************************80

void i4_to_hammersley ( int dim_num, int step, int seed[], int leap[],
  int base[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HAMMERSLEY computes one element of a leaped Hammersley subsequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.  If the base is
//    greater than 1, a standard 1-dimensional
//    van der Corput sequence is generated.  But if the base is
//    negative, this is a signal that the much simpler sequence J/(-BASE)
//    is to be generated.  For the standard Hammersley sequence, the
//    first spatial coordinate uses a base of (-N), and subsequent
//    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
//    This program allows the user to specify any combination of bases,
//    included nonprimes and repeated values.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J M Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, 1960, pages 844-874.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//    1 <= DIM_NUM is required.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required.
//
//    Input, int SEED[DIM_NUM], the Hammersley sequence index corresponding
//    to STEP = 0.
//    0 <= SEED(1:DIM_NUM) is required.
//
//    Input, int LEAP[DIM_NUM], the successive jumps in the Hammersley sequence.
//    1 <= LEAP(1:DIM_NUM) is required.
//
//    Input, int BASE[DIM_NUM], the Hammersley bases.
//
//    Output, double R[DIM_NUM], the STEP-th element of the leaped
//    Hammersley subsequence.
//
{
# define FIDDLE 0.0

  double base_inv;
  int digit;
  int i;
  int seed2;
  int temp;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }

  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  for ( i = 0; i < dim_num; i++ )
  {
    if ( 1 < base[i] )
    {
      seed2 = seed[i] + step * leap[i];

      r[i] = 0.0;

      base_inv = 1.0 / ( ( double ) base[i] );

      while ( seed2 != 0 )
      {
        digit = seed2 % base[i];
        r[i] = r[i] + ( ( double ) digit ) * base_inv;
        base_inv = base_inv / ( ( double ) base[i] );
        seed2 = seed2 / base[i];
      }
    }
//
//  In the following computation, the value of FIDDLE can be:
//
//    0,   for the sequence 0/N, 1/N, ..., N-1/N
//    1,   for the sequence 1/N, 2/N, ..., N/N
//    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
//
    else
    {
      temp = ( seed[i] + step * leap[i] ) % ( -base[i] );
      r[i] = ( ( double ) ( temp ) + FIDDLE ) / ( double ) ( -base[i] );
    }
  }

  return;
# undef FIDDLE
}
//****************************************************************************80

void i4_to_hammersley_sequence ( int dim_num, int n, int step, int seed[],
  int leap[], int base[], double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HAMMERSLEY_SEQUENCE computes N elements of a leaped Hammersley subsequence.
//
//  Discussion:
//
//    The DIM_NUM-dimensional Hammersley sequence is really DIM_NUM separate
//    sequences, each generated by a particular base.  If the base is
//    greater than 1, a standard 1-dimensional
//    van der Corput sequence is generated.  But if the base is
//    negative, this is a signal that the much simpler sequence J/(-BASE)
//    is to be generated.  For the standard Hammersley sequence, the
//    first spatial coordinate uses a base of (-N), and subsequent
//    coordinates use bases of successive primes (2, 3, 5, 7, 11, ...).
//    This program allows the user to specify any combination of bases,
//    included nonprimes and repeated values.
//
//    This routine selects elements of a "leaped" subsequence of the
//    Hammersley sequence.  The subsequence elements are indexed by a
//    quantity called STEP, which starts at 0.  The STEP-th subsequence
//    element is simply element
//
//      SEED(1:DIM_NUM) + STEP * LEAP(1:DIM_NUM)
//
//    of the original Hammersley sequence.
//
//
//    The data to be computed has two dimensions.
//
//    The number of data items is DIM_NUM * N, where DIM_NUM is the spatial dimension
//    of each element of the sequence, and N is the number of elements of the sequence.
//
//    The data is stored in a one dimensional array R.  The first element of the
//    sequence is stored in the first DIM_NUM entries of R, followed by the DIM_NUM entries
//    of the second element, and so on.
//
//    In particular, the J-th element of the sequence is stored in entries
//    0+(J-1)*DIM_NUM through (DIM_NUM-1) + (J-1)*DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    J M Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, 1960, pages 844-874.
//
//    Ladislav Kocis and William Whiten,
//    Computational Investigations of Low-Discrepancy Sequences,
//    ACM Transactions on Mathematical Software,
//    Volume 23, Number 2, 1997, pages 266-294.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of elements of the sequence.
//
//    Input, int STEP, the index of the subsequence element.
//    0 <= STEP is required
//
//    Input, int SEED[DIM_NUM], the Hammersley sequence index corresponding
//    to STEP = 0.
//
//    Input, int LEAP[DIM_NUM], the succesive jumps in the Hammersley sequence.
//
//    Input, int BASE[DIM_NUM], the Hammersley bases.
//
//    Output, double R[DIM_NUM*N], the next N elements of the
//    leaped Hammersley subsequence, beginning with element STEP.
//
{
# define FIDDLE 0.0

  double base_inv;
  int digit;
  int i;
  int j;
  int *seed2;
  int temp;
//
//  Check the input.
//
  if ( !halham_dim_num_check ( dim_num ) )
  {
    exit ( 1 );
  }

  if ( !halham_n_check ( n ) )
  {
    exit ( 1 );
  }

  if ( !halham_step_check ( step ) )
  {
    exit ( 1 );
  }

  if ( !halham_seed_check ( dim_num, seed ) )
  {
    exit ( 1 );
  }

  if ( !halham_leap_check ( dim_num, leap ) )
  {
    exit ( 1 );
  }

  if ( !hammersley_base_check ( dim_num, base ) )
  {
    exit ( 1 );
  }
//
//  Calculate the data.
//
  seed2 = new int[n];

  for ( i = 0; i < dim_num; i++ )
  {
    if ( 1 < base[i] )
    {
      for ( j = 0; j < n; j++ )
      {
        seed2[j] = seed[i] + ( step + j ) * leap[i];
      }

      for ( j = 0; j < n; j++ )
      {
        r[i+j*dim_num] = 0.0;
      }

      for ( j = 0; j < n; j++ )
      {
        base_inv = 1.0 / ( ( double ) base[i] );

        while ( seed2[j] != 0 )
        {
          digit = seed2[j] % base[i];
          r[i+j*dim_num] = r[i+j*dim_num] + ( ( double ) digit ) * base_inv;
          base_inv = base_inv / ( ( double ) base[i] );
          seed2[j] = seed2[j] / base[i];
        }
      }
    }
//
//  In the following computation, the value of FIDDLE can be:
//
//    0,   for the sequence 0/N, 1/N, ..., N-1/N
//    1,   for the sequence 1/N, 2/N, ..., N/N
//    1/2, for the sequence 1/(2N), 3/(2N), ..., (2*N-1)/(2N)
//
    else
    {
      for ( j = 0; j < n; j++ )
      {
        temp = ( seed[i] + ( step + j ) * leap[i] ) % ( -base[i] );

        r[i+j*dim_num] = ( ( double ) ( temp ) + FIDDLE )
                    / ( double ) ( -base[i] );
      }
    }
  }

  delete [] seed2;

  return;
# undef FIDDLE
}
//****************************************************************************80

int i4_uniform_ab ( int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB returns a scaled pseudorandom I4 between A and B.
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
//    02 October 2012
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
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  int c;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4_UNIFORM_AB - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }

  r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
  r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
    +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
  value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
  if ( value < a )
  {
    value = a;
  }
  if ( b < value )
  {
    value = b;
  }

  return value;
}
//****************************************************************************80

void i4vec_transpose_print ( int n, int a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_TRANSPOSE_PRINT prints an I4VEC "transposed".
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//  Example:
//
//    A = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }
//    TITLE = "My vector:  "
//
//    My vector:      1    2    3    4    5
//                    6    7    8    9   10
//                   11
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2004
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
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank or NULL.
//
{
  int i;
  int ihi;
  int ilo;
  int title_len;

  if ( 0 < s_len_trim ( title ) )
  {
    title_len = title.length ( );

    for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
    {
      ihi = i4_min ( ilo + 5 - 1, n );
      if ( ilo == 1 )
      {
        cout << title;
      }
      else
      {
        for ( i = 1; i <= title_len; i++ )
        {
          cout << " ";
        }
      }
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(12) << a[i-1];
      }
      cout << "\n";
    }
  }
  else
  {
    for ( ilo = 1; ilo <= n; ilo = ilo + 5 )
    {
      ihi = i4_min ( ilo + 5 - 1, n );
      for ( i = ilo; i <= ihi; i++ )
      {
        cout << setw(12) << a[i-1];
      }
      cout << "\n";
    }
  }

  return;
}
//****************************************************************************80

void ksub_random2 ( int n, int k, int &seed, int a[] )

//****************************************************************************80
//
//  Purpose:
//
//    KSUB_RANDOM2 selects a random subset of size K from a set of size N.
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
//    FORTRAN77 original version by Albert Nijenhuis, Herbert Wilf.
//    C++ version by John Burkardt.
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
//    Input, int N, the size of the set from which subsets are drawn.
//
//    Input, int K, number of elements in desired subsets.  K must
//    be between 0 and N.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, int A[K].  A(I) is the I-th element of the
//    output set.  The elements of A are in order.
//
{
  int available;
  int candidate;
  int have;
  int need;
  double r;

  if ( k < 0 || n < k )
  {
    cerr << "\n";
    cerr << "KSUB_RANDOM2 - Fatal error!\n";
    cerr << "  N = " << n << "\n";
    cerr << "  K = " << k << "\n";
    cerr << "  but 0 <= K <= N is required!\n";
    exit ( 1 );
  }

  if ( k == 0 )
  {
    return;
  }

  need = k;
  have = 0;
  available = n;
  candidate = 0;

  for ( ; ; )
  {
    candidate = candidate + 1;

    r = r8_uniform_01 ( seed );

    if ( r * ( double ) available <= ( double ) need )
    {
      need = need - 1;
      a[have] = candidate;
      have = have + 1;

      if ( need <= 0 )
      {
        break;
      }

    }

    available = available - 1;

  }

  return;
}
//****************************************************************************80

double *normal ( int dim_num, int n, double r[], double mu[], int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL creates normally distributed points in DIM_NUM space.
//
//  Discussion:
//
//    The multivariate normal distribution for the DIM_NUM dimensional vector X
//    has the form:
//
//      pdf(X) = (2*pi*det(V))^(-DIM_NUM/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
//
//    where MU is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine requires that the user supply the upper triangular
//    Cholesky factor R, which has the property that
//
//      V = R' * R
//
//    This factorization always exists if V is actually symmetric and
//    positive definite.  This factorization can be computed by the
//    routine SPO_FA.
//
//    The user also supplies the mean vector MU.
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
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double R[DIM_NUM*DIM_NUM], the upper triangular Cholesky factor
//    of the variance-covariance matrix.
//
//    Input, double MU[DIM_NUM], the mean vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double NORMAL[DIM_NUM*N], the random points.
//
{
  int i;
  int j;
  int k;
  double *v;
  double *x;

  v = new double[dim_num];
  x = new double[dim_num*n];
//
//  Get a matrix V of normal data.
//  Compute X = MU + R' * V.
//  We actually carry out this computation in the equivalent form X' * R.
//
  for ( j = 0; j < n; j++ )
  {
    r8vec_normal_01 ( dim_num, seed, v );

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = mu[i];
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*dim_num] = x[i+j*dim_num] + v[k] * r[k+i*dim_num];
      }
    }
  }

  delete [] v;

  return x;
}
//****************************************************************************80

double *normal_circular ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CIRCULAR creates circularly normal points in 2 space.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, page 936.
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space, which must be 2.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double NORMAL_CIRULAR[DIM_NUM*N], the random points.
//
{
# define PI 3.141592653589793

  int j;
  double *r;
  double *t;
  double *x;

  r = new double[n];
  t = new double[n];
  x = new double[dim_num*n];
//
//  The angle varies uniformly from 0 to 2 pi.
//
  r8vec_uniform_01 ( n, seed, t );

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * t[j];
  }
//
//  The radius is normally distributed.
//
  r8vec_normal_01 ( n, seed, r );

  for ( j = 0; j < n; j++ )
  {
    x[0+j*dim_num] = r[j] * cos ( t[j] );
    x[1+j*dim_num] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;

  return x;
# undef PI
}
//****************************************************************************80

double *normal_multivariate ( int m, int n, double r[], double mu[],
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MULTIVARIATE samples a multivariate normal distribution.
//
//  Discussion:
//
//    The multivariate normal distribution for the M dimensional vector X
//    has the form:
//
//      pdf(X) = (2*pi*det(V))^(-M/2) * exp(-0.5*(X-MU)'*inverse(V)*(X-MU))
//
//    where MU is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine samples points associated with the M dimensional
//    normal distribution with mean MU and covariance matrix V.
//
//    This routine requires that the user supply the upper triangular
//    Cholesky factor R of V, which has the property that
//
//      V = R' * R
//
//    This factorization always exists if V is actually symmetric and
//    positive definite.  This factorization can be computed by the
//    routine R8PO_FA.
//
//    The user also supplies the mean vector MU.
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
//  Reference:
//
//    Russell Cheng,
//    Random Variate Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998, pages 167-168.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, double R[M*M], the upper triangular Cholesky factor
//    of the variance-covariance matrix.
//
//    Input, double MU[M], the mean vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double NORMAL_MULTIVARIATE[DIM_NUM*N], corresponding
//    points associated with the multivariate normal distribution.
//
{
  int i;
  int j;
  int k;
  double *v;
  double *x;

  v = new double[m];
  x = new double[m*n];
//
//  Compute X = MU + R' * V.
//  We actually carry out this computation in the equivalent form MU + V' * R.
//
  for ( j = 0; j < n; j++ )
  {
    r8vec_normal_01 ( m, seed, v );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = mu[i];
      for ( k = 0; k <= i; k++ )
      {
        x[i+j*m] = x[i+j*m] + v[k] * r[k+i*m];
      }
    }
  }

  delete [] v;

  return x;
}
//****************************************************************************80

double *normal_simple ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_SIMPLE creates normally distributed points in DIM_NUM space.
//
//  Discussion:
//
//    The multivariate normal distribution has the form:
//
//      f(x) = (2*pi*det(V))^(-DIM_NUM/2) * exp(-0.5*(x-mu)'*inverse(V)*(x-mu))
//
//    where mu is the mean vector, and V is a positive definite symmetric
//    matrix called the variance-covariance matrix.
//
//    This routine implements the simplest version of a multivariate
//    normal distribution.  The variance-covariance matrix is the identity,
//    and the mean vector is entirely zero.  Thus, a sample on N points
//    is simply DIM_NUM*N scalar values generated under the univariate
//    normal distribution with zero mean and unit variance.
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
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double NORMAL_SIMPLE[DIM_NUM*N], the random points.
//
{
  double *x;

  x = new double[dim_num*n];

  r8vec_normal_01 ( dim_num * n, seed, x );

  return x;
}
//****************************************************************************80

double *polygon_centroid_2d ( int n, double v[] )

//****************************************************************************80
//
//  Purpose:
//
//    POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
//
//  Formula:
//
//    Denoting the centroid coordinates by CENTROID, then
//
//      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
//      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
//
//    Green's theorem states that
//
//      Integral ( Polygon boundary ) ( M dx + N dy ) =
//      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
//
//    Using M = 0 and N = x * x / 2, we get:
//
//      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
//
//    which becomes
//
//      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
//        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
//
//    where, when I = N, the index "I+1" is replaced by 1.
//
//    A similar calculation gives us a formula for CENTROID(2).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerard Bashein and Paul Detmer,
//    Centroid of a Polygon,
//    Graphics Gems IV, edited by Paul Heckbert,
//    AP Professional, 1994.
//
//  Parameters:
//
//    Input, int N, the number of sides of the polygonal shape.
//
//    Input, double V[2*N], the coordinates of the vertices
//    of the shape.
//
//    Output, double POLYGON_CENTROID_2D[2], the coordinates of the
//    centroid of the shape.
//
{
  double area;
  double *centroid;
  int i;
  int ip1;
  double temp;
//
  area = 0.0;
  centroid = new double[2];
  centroid[0] = 0.0;
  centroid[1] = 0.0;

  for ( i = 0; i < n; i++ )
  {
    if ( i < n-1 )
    {
      ip1 = i + 1;
    }
    else
    {
      ip1 = 0;
    }

    temp = ( v[0+i*2] * v[1+ip1*2] - v[0+ip1*2] * v[1+i*2] );

    area = area + temp;

    centroid[0] = centroid[0] + ( v[0+ip1*2] + v[0+i*2] ) * temp;
    centroid[1] = centroid[1] + ( v[1+ip1*2] + v[1+i*2] ) * temp;

  }

  area = area / 2.0;

  centroid[0] = centroid[0] / ( 6.0 * area );
  centroid[1] = centroid[1] / ( 6.0 * area );

  return centroid;
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
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, pages 870-873.
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
    cerr << "\n";
    cerr << "PRIME - Fatal error!\n";
    cerr << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

float r4_abs ( float x )

//****************************************************************************80
//
//  Purpose:
//
//    R4_ABS returns the absolute value of an R4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, float X, the quantity whose absolute value is desired.
//
//    Output, float R4_ABS, the absolute value of X.
//
{
  float value;

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
//    14 November 2006
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
  int value;

  if ( x < 0.0 )
  {
    value = - ( int ) ( r4_abs ( x ) + 0.5 );
  }
  else
  {
    value =   ( int ) ( r4_abs ( x ) + 0.5 );
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
//    Output double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

double r8_min ( double x, double y )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
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
//    Output double R8_MIN, the minimum of X and Y.
//
{
  if ( x < y )
  {
    return x;
  }
  else
  {
    return y;
  }
}
//****************************************************************************80

int r8_nint ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NINT returns the nearest integer to a double precision real value.
//
//  Example:
//
//        X         R8_NINT
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
//    Input, double X, the real value.
//
//    Output, int R8_NINT, the nearest integer to X.
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

double r8_normal_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, 
//  return one, and save one.
//
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = r8_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
}
//****************************************************************************80

double r8_pi ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_PI returns the value of PI to 16 digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_PI, the value of PI.
//
{
  return 3.141592653589793;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

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
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
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
//    09 April 2012
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
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double *r8mat_normal_01_new ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORMAL_01_NEW returns a unit pseudonormal R8MAT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 April 2012
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
//    Input, int M, N, the number of rows and columns in the array.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, double R8MAT_NORMAL_01_NEW[M*N], the array of pseudonormal values.
//
{
  double *r;

  r = r8vec_normal_01_new ( m * n, seed );

  return r;
}
//****************************************************************************80

void r8mat_print ( int m, int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT, with an optional title.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
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
//    Input, string TITLE, a title to be printed.
//
{
  r8mat_print_some ( m, n, a, 1, 1, m, n, title );

  return;
}
//****************************************************************************80

void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,
  int jhi, string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values,  stored as a vector
//    in column-major order.
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
//    Input, string TITLE, a title for the matrix.
//
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

  return;
# undef INCX
}
//****************************************************************************80

double *r8mat_uniform_01_new ( int m, int n, int &seed )

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
//    Input/output, int &SEED, the "seed" value.  Normally, this
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
      k = seed / 127773;

      seed = 16807 * ( seed - k * 127773 ) - k * 2836;

      if ( seed < 0 )
      {
        seed = seed + 2147483647;
      }
      r[i+j*m] = ( double ) ( seed ) * 4.656612875E-10;
    }
  }

  return r;
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
      output << "  " << setw(10) << table[i+j*m];
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

int r8po_fa ( double a[], int lda, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_FA factors a real symmetric positive definite matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch and Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input/output, double A[LDA*N].  On input, the symmetric matrix
//    to be  factored.  Only the diagonal and upper triangle are used.
//    On output, an upper triangular matrix R so that A = R'*R
//    where R' is the transpose.  The strict lower triangle is unaltered.
//    If INFO /= 0, the factorization is not complete.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Output, int R8PO_FA, error flag.
//    0, for normal return.
//    K, signals an error condition.  The leading minor of order K is not
//    positive definite.
//
{
  int info;
  int j;
  int k;
  double s;
  double t;

  for ( j = 1; j <= n; j++ )
  {
    s = 0.0;

    for ( k = 1; k <= j-1; k++ )
    {
      t = a[k-1+(j-1)*lda] - ddot ( k-1, a+0+(k-1)*lda, 1, a+0+(j-1)*lda, 1 );
      t = t / a[k-1+(k-1)*lda];
      a[k-1+(j-1)*lda] = t;
      s = s + t * t;
    }

    s = a[j-1+(j-1)*lda] - s;

    if ( s <= 0.0 )
    {
      info = j;
      return info;
    }

    a[j-1+(j-1)*lda] = sqrt ( s );
  }

  info = 0;

  return info;
}
//****************************************************************************80

void r8po_sl ( double a[], int lda, int n, double b[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8PO_SL solves a linear system factored by R8PO_FA.
//
//  Discussion:
//
//    A division by zero will occur if the input factor contains
//    a zero on the diagonal.  Technically this indicates
//    singularity but it is usually caused by improper subroutine
//    arguments.  It will not occur if the subroutines are called
//    correctly and INFO == 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    FORTRAN77 original version by Dongarra, Moler, Bunch, Stewart.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Dongarra, Moler, Bunch and Stewart,
//    LINPACK User's Guide,
//    SIAM, (Society for Industrial and Applied Mathematics),
//    3600 University City Science Center,
//    Philadelphia, PA, 19104-2688.
//    ISBN 0-89871-172-X
//
//  Parameters:
//
//    Input, double A[LDA*N], the output from R8PO_FA.
//
//    Input, int LDA, the leading dimension of the array A.
//
//    Input, int N, the order of the matrix.
//
//    Input/output, double B[N].  On input, the right hand side.
//    On output, the solution.
//
{
  int k;
  double t;
//
//  Solve R' * Y = B.
//
  for ( k = 1; k <= n; k++ )
  {
    t = ddot ( k-1, a+0+(k-1)*lda, 1, b, 1 );
    b[k-1] = ( b[k-1] - t ) / a[k-1+(k-1)*lda];
  }
//
//  Solve R * X = Y.
//
  for ( k = n; 1 <= k; k-- )
  {
    b[k-1] = b[k-1] / a[k-1+(k-1)*lda];
    t = -b[k-1];
    daxpy ( k-1, t, a+0+(k-1)*lda, 1, b, 1 );
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
//      R8VEC_NORM = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
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

void r8vec_normal_01 ( int n, int &seed, double x[] )

//****************************************************************************80
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
//    Input/output, int &SEED, a seed for the random number generator.
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

double *r8vec_normal_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
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
//    02 February 2005
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, double Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  double pi = 3.141592653589793;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the internal data.
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
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r8vec_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT prints an R8VEC.
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
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title to be printed first.
//    TITLE may be blank.
//
{
  int i;

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(8)  << i
         << "  " << setw(14) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
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
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
}
//****************************************************************************80

void r8vec_uniform_01 ( int n, int &seed, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *r8vec_zero_new ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_ZERO_NEW creates and zeroes an R8VEC.
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
//    10 July 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Output, double R8VEC_ZERO_NEW[N], a vector of zeroes.
//
{
  double *a;
  int i;

  a = new double[n];

  for ( i = 0; i < n; i++ )
  {
    a[i] = 0.0;
  }
  return a;
}
//****************************************************************************80

unsigned long random_initialize ( unsigned long seed )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_INITIALIZE initializes the RANDOM random number generator.
//
//  Discussion:
//
//    If you don't initialize RANDOM, the random number generator,
//    it will behave as though it were seeded with value 1.
//    This routine will either take a user-specified seed, or
//    (if the user passes a 0) make up a "random" one.  In either
//    case, the seed is passed to SRANDOM (the appropriate routine
//    to call when setting the seed for RANDOM).  The seed is also
//    returned to the user as the value of the function.
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
//    Input, unsigned long SEED, is either 0, which means that the user
//    wants this routine to come up with a seed, or nonzero, in which
//    case the user has supplied the seed.
//
//    Output, unsigned long RANDOM_INITIALIZE, is the value of the seed
//    passed to SRANDOM, which is either the user's input value, or if
//    that was zero, the value selected by this routine.
//
{
# define DEBUG 0

  if ( seed != 0 )
  {
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with user SEED = " << seed << "\n";
    }
  }
  else
  {
    seed = get_seed ( );
    if ( DEBUG )
    {
      cout << "\n";
      cout << "RANDOM_INITIALIZE\n";
      cout << "  Initialize RANDOM with arbitrary SEED = " << seed << "\n";
    }
  }
//
//  Now set the seed.
//
  srandom ( seed );

  return seed;
# undef DEBUG
}
//****************************************************************************80

int s_len_trim ( string s )

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
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

void scale_from_simplex01 ( int dim_num, int n, double t[], double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SCALE_FROM_SIMPLEX01 rescales data from a unit to non-unit simplex.
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
//  Reference:
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
//    Input, double T[DIM_NUM*(DIM_NUM+1)], the coordinates of the DIM_NUM+1
//    points that define the simplex.  T[0:DIM_NUM-1,0] corresponds to the
//    origin, and T[0:DIM_NUM-1,J] will be the image of the J-th unit 
//    coordinate vector.
//
//    Input/output, double X[DIM_NUM*N], the data to be modified.
//
{
  double *a;
  int i;
  int j;
  double *v;

  a = new double[dim_num*dim_num];
  v = new double[dim_num];

  for ( j = 0; j < dim_num; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      a[i+j*dim_num] = t[i+j*(dim_num+1)] - t[i+0*(dim_num+1)];
    }
  }

  for ( j = 0; j < n; j++ )
  {

    for ( i = 0; i < dim_num; i++ )
    {
      v[i] = x[i+j*dim_num];
    }

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = t[i+0*(dim_num+1)];
      for ( j = 0; j < n; j++ )
      {
        x[i+j*dim_num] = x[i+j*dim_num] + a[i+j*dim_num] * v[j];
      }
    }

  }

  delete [] a;
  delete [] v;

  return;
}
//****************************************************************************80

void scale_to_ball01 ( int dim_num, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SCALE_TO_BALL01 translates and rescales data to fit within the unit ball.
//
//  Discussion:
//
//    Completely arbitrary input data is given.
//
//    The average of the data is computed, and taken as the coordinates
//    of the center C of a sphere.  The radius R of that sphere is the
//    distance from the center to the furthest point in the data set.
//
//    Then each point is transformed to the ball of center 0 and radius
//    1 by subtracting C and dividing by R:
//
//      X(1:DIM_NUM,J) -> ( X(1:DIM_NUM,J) - C(1:DIM_NUM) ) / R
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[DIM_NUM*N], the data to be modified.
//
{
  int i;
  int j;
  double r;
  double scale;
  double *xave;
//
//  Determine the center.
//
  xave = new double[dim_num];

  for ( i = 0; i < dim_num; i++ )
  {
    xave[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      xave[i] = xave[i] + x[i+j*dim_num];
    }
    xave[i] = xave[i] / ( double ) n;
  }
//
//  Determine the maximum distance of any point from the center.
//
  for ( j = 0; j < n; j++ )
  {
    r = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      r = r + pow ( x[i+j*dim_num] - xave[i], 2 );
    }
    if ( scale < r )
    {
      scale = r ;
    }
  }

  scale = sqrt ( scale );

  if ( 0.0 < scale )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < dim_num; i++)
      {
        x[i+j*dim_num] = ( x[i+j*dim_num] - xave[i] ) / scale;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < dim_num; i++)
      {
        x[i+j*dim_num] = 0.0;
      }
    }
  }

  delete [] xave;

  return;
}
//****************************************************************************80

void scale_to_block01 ( int dim_num, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SCALE_TO_BLOCK01 translates and rescales data to fit in the unit block.
//
//  Discussion:
//
//    The minimum and maximum coordinate values M1(I) and M2(I) are
//    determined, and the maximum of M2(I) - M1(I) is used to scale
//    all the coordinates by the same factor.
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
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[DIM_NUM*N], the data to be modified.
//
{
  int i;
  int j;
  double *xmax;
  double *xmin;
  double xrange;
  double xrange2;

  xmax = new double[dim_num];
  xmin = new double[dim_num];
//
//  Determine the extremes in each dimension.
//
  xrange = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    xmin[i] = x[i+0*dim_num];
    xmax[i] = x[i+0*dim_num];
    for ( j = 1; j < n; j++ )
    {
      xmin[i] = r8_min ( xmin[i], x[i+j*dim_num] );
      xmax[i] = r8_max ( xmax[i], x[i+j*dim_num] );
    }
    xrange = r8_max ( xrange, xmax[i] - xmin[i] );
  }
//
//  Extend all the extremes so that the range is the same in each dimension.
//
  for ( i = 0; i < dim_num; i++ )
  {
    xrange2 = xrange - ( xmax[i] - xmin[i] );
    xmax[i] = xmax[i] + 0.5 * xrange2;
    xmin[i] = xmin[i] - 0.5 * xrange2;
  }
//
//  Now map the data to [0,1], using a single dilation factor
//  for all dimensions.
//
  if ( 0.0 == xrange )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        x[i+j*dim_num] = 0.5;
      }
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        x[i+j*dim_num] = ( x[i+j*dim_num] - xmin[i] ) / xrange;
      }
    }
  }

  delete [] xmax;
  delete [] xmin;

  return;
}
//****************************************************************************80

void scale_to_cube01 ( int dim_num, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    SCALE_TO_CUBE01 translates and rescales data to the unit hypercube.
//
//  Discussion:
//
//    In each coordinate dimension I, the minimum and maximum coordinate
//    values M1(I) and M2(I) are determined.
//
//    Then, in each coordinate, the points are rescaled as
//
//      X(I) -> ( X(I) - M1(I) ) / ( M2(I) - M1(I) ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, double X[DIM_NUM*N], the data to be modified.
//
{
  int i;
  int j;
  double xmax;
  double xmin;

  for ( i = 0; i < dim_num; i++ )
  {
    xmin = x[i+0*dim_num];
    xmax = x[i+0*dim_num];
    for ( j = 1; j < n; j++ )
    {
      if ( x[i+j*dim_num] < xmin )
      {
        xmin = x[i+j*dim_num];
      }
      if ( xmax < x[i+j*dim_num] )
      {
        xmax = x[i+j*dim_num];
      }
    }

    if ( 0.0 < xmax - xmin )
    {
      for ( j = 0; j < n; j++ )
      {
        x[i+j*dim_num] = ( x[i+j*dim_num] - xmin ) / ( xmax - xmin );
      }
    }
    else
    {
      for ( j = 0; j < n; j++ )
      {
        x[i+j*dim_num] = 0.0;
      }
    }
  }

  return;
}
//****************************************************************************80

double stri_angles_to_area ( double r, double a, double b, double c )

//****************************************************************************80
//
//  Purpose:
//
//    STRI_ANGLES_TO_AREA computes the area of a spherical triangle.
//
//  Discussion:
//
//    A sphere centered at 0 in 3D satisfies the equation:
//
//      X**2 + Y**2 + Z**2 = R**2
//
//    A spherical triangle is specified by three points on the surface
//    of the sphere.
//
//    The area formula is known as Girard's formula.
//
//    The area of a spherical triangle is:
//
//      AREA = ( A + B + C - PI ) * R**2
//
//    where A, B and C are the (surface) angles of the triangle.
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
//    Input, double R, the radius of the sphere.
//
//    Input, double A, B, C, the angles of the triangle.
//
//    Output, double STRI_ANGLES_TO_AREA_3D, the area of the spherical triangle.
//
{
  double area;
  double pi = 3.141592653589793;

  area = r * r * ( a + b + c - pi );

  return area;
}
//****************************************************************************80

void stri_sides_to_angles ( double r, double as, double bs, double cs,
  double *a, double *b, double *c )

//****************************************************************************80
//
//  Purpose:
//
//    STRI_SIDES_TO_ANGLES computes spherical triangle angles.
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
//    Input, double R, the radius of the sphere.
//
//    Input, double AS, BS, CS, the (geodesic) length of the sides of the
//    triangle.
//
//    Output, double *A, *B, *C, the spherical angles of the triangle.
//    Angle A is opposite the side of length AS, and so on.
//
{
  double asu;
  double bsu;
  double csu;
  double ssu;
  double tan_a2;
  double tan_b2;
  double tan_c2;

  asu = as / r;
  bsu = bs / r;
  csu = cs / r;
  ssu = ( asu + bsu + csu ) / 2.0;

  tan_a2 = sqrt ( ( sin ( ssu - bsu ) * sin ( ssu - csu ) ) /
                  ( sin ( ssu ) * sin ( ssu - asu )     ) );

  *a = 2.0 * atan ( tan_a2 );

  tan_b2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - csu ) ) /
                  ( sin ( ssu ) * sin ( ssu - bsu )     ) );

  *b = 2.0 * atan ( tan_b2 );

  tan_c2 = sqrt ( ( sin ( ssu - asu ) * sin ( ssu - bsu ) ) /
                  ( sin ( ssu ) * sin ( ssu - csu )     ) );

  *c = 2.0 * atan ( tan_c2 );

  return;
}
//****************************************************************************80

void stri_vertices_to_sides ( double r, double v1[3], double v2[3],
  double v3[3], double *as, double *bs, double *cs )

//****************************************************************************80
//
//  Purpose:
//
//    STRI_VERTICES_TO_SIDES_3D computes spherical triangle sides.
//
//  Discussion:
//
//    We can use the ACOS system call here, but the ARC_COSINE routine
//    will automatically take care of cases where the input argument is
//    (usually slightly) out of bounds.
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
//    Input, double R, the radius of the sphere.
//
//    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
//    triangle.
//
//    Output, double *AS, *BS, *CS, the (geodesic) length of the sides of the
//    triangle.
//
{
  *as = r * arc_cosine ( r8vec_dot_product ( 3, v2, v3 ) / ( r * r ) );
  *bs = r * arc_cosine ( r8vec_dot_product ( 3, v3, v1 ) / ( r * r ) );
  *cs = r * arc_cosine ( r8vec_dot_product ( 3, v1, v2 ) / ( r * r ) );

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
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
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
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.  AREA will
//    be nonnegative.
//
{
  double area;

  area = fabs ( 0.5 * (
    t[0+0*2] * ( t[1+2*2] - t[1+1*2] ) +
    t[0+1*2] * ( t[1+0*2] - t[1+2*2] ) +
    t[0+2*2] * ( t[1+1*2] - t[1+0*2] ) ) );

  return area;
}
//****************************************************************************80

void tuple_next_fast ( int m, int n, int rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank          X
//    ----          ----
//   -1            -1 -1
//
//    0             1  1
//    1             1  2
//    2             1  3
//    3             2  1
//    4             2  2
//    5             2  3
//    6             3  1
//    7             3  2
//    8             3  3
//    9             1  1
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
//    Input, int M, the maximum entry in each component.
//    M must be greater than 0.
//
//    Input, int N, the number of components.
//    N must be greater than 0.
//
//    Input, integer RANK, indicates the rank of the tuples.
//    Typically, 0 <= RANK < N**M; values larger than this are legal
//    and meaningful, and are equivalent to the corresponding value
//    MOD N**M.  If RANK < 0, this indicates that this is the first
//    call for the given values of (M,N).  Initialization is done,
//    and X is set to a dummy value.
//
//    Output, int X[N], the next tuple, or a dummy value if initialization
//    is being done.
//
{
  static int *base = NULL;
  int i;
//
  if ( rank < 0 )
  {
    if ( m <= 0 )
    {
      cerr << "\n";
      cerr << "TUPLE_NEXT_FAST - Fatal error!\n";
      cerr << "  The value M <= 0 is not legal.\n";
      cerr << "  M = " << m << "\n";
      exit ( 1 );
    }
    if ( n <= 0 )
    {
      cerr << "\n";
      cerr << "TUPLE_NEXT_FAST - Fatal error!\n";
      cerr << "  The value N <= 0 is not legal.\n";
      cerr << "  N = " << n << "\n";
      exit ( 1 );
    }

    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
    for ( i = 0; i < n; i++ )
    {
      x[i] = -1;
    }
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      x[i] = ( ( rank / base[i] ) % m ) + 1;
    }
  }
  return;
}
//****************************************************************************80

double *uniform_in_annulus ( double pc[], double r1, double r2, int n,
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_ANNULUS samples a circular annulus.
//
//  Discussion:
//
//    A circular annulus with center PC, inner radius R1 and
//    outer radius R2, is the set of points P so that
//
//      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, double PC[2], the center.
//
//    Input, double R1, R2, the inner and outer radii.
//
//    Input, int N, the number of points to generate.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ANNULUS[2*N], sample points.
//
{
# define DIM_NUM 2
# define PI 3.141592653589793

  int j;
  double *p;
  double r;
  double theta;
  double u;
  double v;

  p = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    u = r8_uniform_01 ( seed );
    theta = u * 2.0 * PI;
    v = r8_uniform_01 ( seed );
    r = sqrt ( ( 1.0 - v ) * r1 * r1
           +           v   * r2 * r2 );

    p[0+j*2] = pc[0] + r * cos ( theta );
    p[1+j*2] = pc[1] + r * sin ( theta );
  }

  return p;
# undef DIM_NUM
# undef PI
}
//****************************************************************************80

double *uniform_in_annulus_accept ( double pc[], double r1, double r2, int n,
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_ANNULUS_ACCEPT accepts points in an annulus.
//
//  Discussion:
//
//    A circular annulus with center PC, inner radius R1 and
//    outer radius R2, is the set of points P so that
//
//      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
//
//    The acceptance/rejection method is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double PC[2], the center.
//
//    Input, double R1, R2, the inner and outer radii.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ANNULUS_ACCEPT[2*N], the points.
//
{
# define DIM_NUM 2

  int i;
  int j;
  double *p;
  double r_squared;
  double u[DIM_NUM];

  if ( r2 <= r1 )
  {
    cerr << "\n";
    cerr << "UNIFORM_IN_ANNULUS_ACCEPT - Fatal error!\n";
    cerr << "  R2 <= R1.\n";
    exit ( 1 );
  }

  p = new double[DIM_NUM*n];
//
//  Generate points in a square of "radius" R2.
//  Accept those points which lie inside the circle of radius R2, and outside
//  the circle of radius R1.
//
  for ( j = 0; j < n; j++ )
  {
    for ( ; ; )
    {
      r8vec_uniform_01 ( DIM_NUM, seed, u );

      for ( i = 0; i < DIM_NUM; i++ )
      {
        u[i] = ( 2.0 * u[i] - 1.0 ) * r2;
      }

      r_squared = 0.0;
      for ( i = 0; i < DIM_NUM; i++ )
      {
        r_squared = r_squared + u[i] * u[i];
      }

      if ( r1 * r1 <= r_squared && r_squared <= r2 * r2 )
      {
        break;
      }
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      p[i+j*DIM_NUM] = pc[i] + u[i];
    }

  }

  return p;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_annulus_sector ( double pc[], double r1, double r2,
  double theta1, double theta2, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_ANNULUS_SECTOR samples an annular sector in 2D.
//
//  Discussion:
//
//    An annular sector with center PC, inner radius R1 and
//    outer radius R2, and angles THETA1, THETA2, is the set of points
//    P so that
//
//      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
//
//    and
//
//      THETA1 <= THETA ( P - PC ) <= THETA2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, double PC[2], the center.
//
//    Input, double R1, R2, the inner and outer radii.
//
//    Input, double THETA1, THETA2, the angles.
//
//    Input, int N, the number of points to generate.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ANNULUS_SECTOR[2*N], sample points.
//
{
# define DIM_NUM 2

  int j;
  double *p;
  double r;
  double theta;
  double u;
  double v;

  p = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    u = r8_uniform_01 ( seed );

    theta = ( 1.0 - u ) * theta1
          +         u   * theta2;

    v = r8_uniform_01 ( seed );

    r = sqrt ( ( 1.0 - v ) * r1 * r1
           +           v   * r2 * r2 );

    p[0+j*2] = pc[0] + r * cos ( theta );
    p[1+j*2] = pc[1] + r * sin ( theta );
  }

  return p;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_circle01_map ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_CIRCLE01_MAP maps uniform points into the unit circle.
//
//  Discussion:
//
//    The unit circle is centered at the origin and has radius 1.
//
//    This routine is valid for spatial dimension DIM_NUM = 2.
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
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_CIRCLE01_MAP[DIM_NUM*N], the points.
//
{
# define DIM_NUM 2
# define PI 3.141592653589793

  int j;
  double *r;
  double *t;
  double *x;

  r = new double[n];
  t = new double[n];
  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    r[j] = r8_uniform_01 ( seed );
    r[j] = sqrt ( r[j] );
  }

  for ( j = 0; j < n; j++ )
  {
    t[j] = 2.0 * PI * r8_uniform_01 ( seed );
  }

  for ( j = 0; j < n; j++ )
  {
    x[0+j*DIM_NUM] = r[j] * cos ( t[j] );
    x[1+j*DIM_NUM] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;

  return x;
# undef DIM_NUM
# undef PI
}
//****************************************************************************80

double *uniform_in_cube01 ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_CUBE01 creates uniform points in the unit hypercube.
//
//  Discussion:
//
//    The unit hypercube is defined as points all of whose components are between
//    0 and 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_CUBE01[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  double *x;

  x = new double[dim_num*n];

  r8vec_uniform_01 ( dim_num*n, seed, x );

  return x;
}
//****************************************************************************80

double *uniform_in_ellipsoid_map ( int dim_num, int n, double a[], double r,
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_ELLIPSOID_MAP maps uniform points into an ellipsoid.
//
//  Discussion:
//
//    The points X in the ellipsoid are described by a DIM_NUM by DIM_NUM positive
//    definite symmetric matrix A, and a "radius" R, such that
//
//      X' * A * X <= R * R
//
//    The algorithm computes the Cholesky factorization of A:
//
//      A = U' * U.
//
//    A set of uniformly random points Y is generated, satisfying:
//
//      Y' * Y <= R * R.
//
//    The appropriate points in the ellipsoid are found by solving
//
//      U * X = Y
//
//    Thanks to Dr Karl-Heinz Keil for pointing out that the original
//    coding was actually correct only if A was replaced by its inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Input, double A[DIM_NUM*DIM_NUM], the matrix that describes the ellipsoid.
//
//    Input, double R, the right hand side of the ellipsoid equation.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_ELLIPSOID_MAP[DIM_NUM*N], the points.
//
{
  int i;
  int info;
  int j;
  int k;
  double *u;
  double *x;
//
//  Get the upper triangular Cholesky factor U of A.
//
  u = new double[dim_num*dim_num];

  for ( j = 0; j < dim_num; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      u[i+j*dim_num] = a[i+j*dim_num];
    }
  }

  info = r8po_fa ( u, dim_num, dim_num );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "UNIFORM_IN_ELLIPSOID_MAP - Fatal error!\n";
    cerr << "  R8PO_FA reports that the matrix A\n";
    cerr << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }
//
//  Get the points Y that satisfy Y' * Y <= R * R.
//
  x = uniform_in_sphere01_map ( dim_num, n, seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = r * x[i+j*dim_num];
    }
  }
//
//  Solve U * X = Y.
//
  for ( j = 0; j < n; j++ )
  {
    r8po_sl ( u, dim_num, dim_num, x+j*dim_num );
  }

  delete [] u;

  return x;
}
//****************************************************************************80

double *uniform_in_parallelogram_map ( double v1[2], double v2[2],
  double v3[2], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_PARALLELOGRAM_MAP maps uniform points into a parallelogram.
//
//  Discussion:
//
//    The parallelogram is defined by three vertices, V1, V2 and V3.
//    The missing vertex V4 is equal to V2+V3-V1.
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
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_PARALLELOGRAM_MAP[2*N], the points.
//
{
# define DIM_NUM 2

  int i;
  int j;
  double r;
  double s;
  double *x;

  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    r = r8_uniform_01 ( seed );
    s = r8_uniform_01 ( seed );

    for ( i = 0; i < DIM_NUM; i++ )
    {
      x[i+j*DIM_NUM] = ( 1.0 - r - s ) * v1[i]
                             + r       * v2[i]
                                 + s   * v3[i];
    }
  }

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_polygon_map ( int nv, double v[], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_POLYGON_MAP maps uniform points into a polygon.
//
//  Discussion:
//
//    If the polygon is regular, or convex, or at least star-shaped,
//    this routine will work.
//
//    This routine assumes that all points between the centroid and
//    any point on the boundary lie within the polygon.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NV, the number of vertices.
//
//    Input, double V[2*NV], the vertices of the polygon, listed in
//    clockwise or counterclockwise order.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_POLYGON_MAP[2*N], the points.
//
{
# define DIM_NUM 2

  double *area;
  double *centroid;
  int i;
  int ip1;
  int j;
  int k;
  double r[DIM_NUM];
  double t[2*3];
  int triangle;
  int triangle_p1;
  double total;
  double u;
  double *x;

  area = new double[nv];
  x = new double[DIM_NUM*n];
//
//  Find the centroid.
//
  centroid = polygon_centroid_2d ( nv, v );
//
//  Determine the areas of each triangle.
//
  total = 0.0;
  for ( i = 0; i < nv; i++ )
  {
    if ( i < nv-1 )
    {
      ip1 = i + 1;
    }
    else
    {
      ip1 = 0;
    }

    t[0+0*2] = v[0+i*2];
    t[1+0*2] = v[1+i*2];

    t[0+1*2] = v[0+ip1*2];
    t[1+1*2] = v[1+ip1*2];

    t[0+2*2] = centroid[0];
    t[1+2*2] = centroid[1];

    area[i] = triangle_area_2d ( t );

    total = total + area[i];
  }
//
//  Normalize the areas.
//
  for ( i = 0; i < nv; i++ )
  {
    area[i] = area[i] / total;
  }
//
//  Replace each area by the sum of itself and all previous ones.
//
  for ( i = 1; i < nv; i++ )
  {
    area[i] = area[i] + area[i-1];
  }

  for ( j = 0; j < n; j++ )
  {
//
//  Choose a triangle T at random, based on areas.
//
    u = r8_uniform_01 ( seed );

    for ( k = 0; k < nv; k++ )
    {
      if ( u <= area[k] )
      {
        triangle = k;
        break;
      }
    }
//
//  Now choose a point at random in the triangle.
//
    if ( triangle < nv-1 )
    {
      triangle_p1 = triangle + 1;
    }
    else
    {
      triangle_p1 = 0;
    }

    r8vec_uniform_01 ( DIM_NUM, seed, r );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      x[i+j*DIM_NUM] = ( 1.0 - r[0] - r[1] ) * v[i+DIM_NUM*triangle]
                       +       r[0]          * v[i+DIM_NUM*triangle_p1]
                       +              r[1]   * centroid[i];
    }

  }

  delete [] area;
  delete [] centroid;

  return x;
# undef M
}
//****************************************************************************80

double *uniform_in_sector_map ( double r1, double r2, double t1,
  double t2, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_SECTOR_MAP maps uniform points into a circular sector.
//
//  Discussion:
//
//    The sector lies between circles with center at 0 and radius R1 and R2,
//    and between rays from the center at the angles T1 and T2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, double R1, R2, the two radii.
//
//    Input, double T1, T2, the two angles, which should
//    be measured in radians, with T1 < T2.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_SECTOR_MAP[2*N], the points.
//
{
# define DIM_NUM 2

  int j;
  double *r;
  double *t;
  double *u;
  double *v;
  double *x;

  x = new double[DIM_NUM*n];
  r = new double[n];
  t = new double[n];
  u = new double[n];
  v = new double[n];

  r8vec_uniform_01 ( n, seed, u );
  r8vec_uniform_01 ( n, seed, v );

  for ( j = 0; j < n; j++ )
  {
    t[j] = ( 1.0 - u[j] ) * t1 + u[j] * t2;
    r[j] = sqrt ( ( 1.0 - v[j] ) * r1 * r1 + v[j] * r2 * r2 );

    x[0+j*DIM_NUM] = r[j] * cos ( t[j] );
    x[1+j*DIM_NUM] = r[j] * sin ( t[j] );
  }

  delete [] r;
  delete [] t;
  delete [] u;
  delete [] v;

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_simplex01_map ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_SIMPLEX01 maps uniform points into the unit simplex.
//
//  Discussion:
//
//    The interior of the unit DIM_NUM dimensional simplex is the set of points X(1:DIM_NUM)
//    such that each X(I) is nonnegative, and sum(X(1:DIM_NUM)) <= 1.
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
//  Reference:
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_SIMPLEX01_MAP[DIM_NUM*N], the points.
//
{
  double *e;
  int i;
  int j;
  double total;
  double *x;
//
//  The construction begins by sampling DIM_NUM+1 points from the
//  exponential distribution with parameter 1.
//
  e = new double[dim_num+1];
  x = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( dim_num+1, seed, e );

    for ( i = 0; i <= dim_num; i++ )
    {
      e[i] = -log ( e[i] );
    }

    total = 0.0;
    for ( i = 0; i <= dim_num; i++ )
    {
      total = total + e[i];
    }

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+dim_num*j] = e[i] / total;
    }

  }

  delete [] e;

  return x;
}
//****************************************************************************80

double *uniform_in_sphere01_map ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_SPHERE01_MAP maps uniform points into the unit sphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    We first generate a point ON the sphere, and then distribute it
//    IN the sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2014
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
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double X[M*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *v;
  double *x;

  exponent = 1.0 / ( double ) ( m );

  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    v = r8vec_normal_01_new ( m, seed );
//
//  Compute the length of the vector.
//
    norm = r8vec_norm ( m, v );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      v[i] = v[i] / norm;
    }
//
//  Now compute a value to map the point ON the sphere INTO the sphere.
//
    r = r8_uniform_01 ( seed );
    r = pow ( r, exponent );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * v[i];
    }

    delete [] v;
  }

  return x;
}
//****************************************************************************80

double *uniform_in_tetrahedron ( double v[], int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_TETRAHEDRON returns uniform points in a tetrahedron.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Claudio Rocchini, Paolo Cignoni,
//    Generating Random Points in a Tetrahedron,
//    Journal of Graphics Tools,
//    Volume 5, Number 5, 2000, pages 9-12.
//
//  Parameters:
//
//    Input, double V[3*4], the vertices of the tetrahedron.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double UNIFORM_IN_TETRAHEDRON[3*N], the points.
//
{
  double c[4];
  int i;
  int j;
  int k;
  double t;
  double *x;

  x = new double[3*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( 3, seed, c );

    if ( 1.0 < c[0] + c[1] )
    {
      c[0] = 1.0 - c[0];
      c[1] = 1.0 - c[1];
    }

    if ( 1.0 < c[1] + c[2] )
    {
      t = c[2];
      c[2] = 1.0 - c[0] - c[1];
      c[1] = 1.0 - t;
    }
    else if ( 1.0 < c[0] + c[1] + c[2] )
    {
       t = c[2];
       c[2] = c[0] + c[1] + c[2] - 1.0;
       c[0] = 1.0 - c[1] - t;
    }
    c[3] = 1.0 - c[0] - c[1] - c[2];

    for ( i = 0; i < 3; i++ )
    {
      x[i+j*3] = 0.0;
      for ( k = 0; k < 4; k++ )
      {
        x[i+j*3] = x[i+j*3] + v[i+k*3] * c[k];
      }
    }
  }

  return x;
}
//****************************************************************************80

double *uniform_in_triangle_map1 ( double v1[2], double v2[2], double v3[2],
  int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE_MAP1 maps uniform points into a triangle.
//
//  Discussion:
//
//    This routine uses Turk's rule 1.
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
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_TRIANGLE_MAP1[2*N], the points.
//
{
# define DIM_NUM 2

  double a;
  double b;
  double c;
  int i;
  int j;
  double r[DIM_NUM];
  double total;
  double *x;

  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( DIM_NUM, seed, r );

    r[1] = sqrt ( r[1] );

    a = 1.0 - r[1];
    b = ( 1.0 - r[0] ) * r[1];
    c = r[0] * r[1];

    for ( i = 0; i < DIM_NUM; i++ )
    {
      x[i+j*DIM_NUM] = a * v1[i] + b * v2[i] + c * v3[i];
    }
  }

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_triangle_map2 ( double v1[2], double v2[2], double v3[2],
  int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE_MAP2 maps uniform points into a triangle.
//
//  Discussion:
//
//    The triangle is defined by three vertices.
//
//    This routine uses Turk's rule 2.
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
//  Reference:
//
//    Greg Turk,
//    Generating Random Points in a Triangle,
//    in Graphics Gems,
//    edited by Andrew Glassner,
//    AP Professional, 1990, pages 24-28.
//
//  Parameters:
//
//    Input, double V1[2], V2[2], V3[2], the vertices.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_IN_TRIANGLE_MAP2[2*N], the points.
//
{
# define DIM_NUM 2

  int i;
  int j;
  double r[DIM_NUM];
  double total;
  double *x;

  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( DIM_NUM, seed, r );

    if ( 1.0 < r[0] + r[1] )
    {
      r[0] = 1.0 - r[0];
      r[1] = 1.0 - r[1];
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      x[i+j*DIM_NUM] = ( 1.0 - r[0] - r[1] ) * v1[i]
                     +             r[0]          * v2[i]
                     +                    r[1]   * v3[i];
    }
  }

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_in_triangle01_map ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_IN_TRIANGLE01_MAP maps uniform points into the unit triangle.
//
//  Discussion:
//
//    The triangle is defined by the three vertices (1,0), (0,1) and (0,0).
//    Because this is a right triangle, it is easy to generate sample points.
//    In the case of a general triangle, more care must be taken.
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
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number geneator.
//
//    Output, double UNIFORM_IN_TRIANGLE01_MAP[2*N], the points.
//
{
# define DIM_NUM 2

  int i;
  int j;
  double r[DIM_NUM];
  double total;
  double *x;

  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( DIM_NUM, seed, r );

    total = 0.0;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      total = total + r[i];
    }
    if ( 1.0 < total )
    {
      for ( i = 0; i < DIM_NUM; i++ )
      {
        r[i] = 1.0 - r[i];
      }
    }

    for ( i = 0; i < DIM_NUM; i++ )
    {
      x[i+j*DIM_NUM] = r[i];
    }
  }

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_on_cube ( int m, int n, double c[], double r, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_CUBE returns random points on the surface of a cube.
//
//  Discussion:
//
//    The cube is assumed to be aligned with the coordinate axes.
//
//    The cube has center C and radius R.  Any point on the surface of
//    the cube is described by
//
//      X = C + R * PM
//
//    where PM is an M-dimensional vector whose entries are between
//    -1 and +1, and for which at least one value has norm 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= M.
//
//    Input, int N, the number of points.
//
//    Input, double C[M], the coordinates of the center.
//
//    Input, double R, the radius.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double UNIFORM_ON_CUBE[M*N], the coordinates of N points, chosen
//    uniformly at random from the surface of the M-cube of center C and 
//    radius R.
//
{
  int i;
  int j;
  int k;
  double *x;
//
//  Choose random points within the cube of radius 1.
//
  x = r8mat_uniform_01_new ( m, n, seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = 2.0 * x[i+j*m] - 1.0;
    }
  }
//
//  For each point, select a coordinate at random, and set it to +1 or -1.
//
  for ( j = 0; j < n; j++ )
  {
    i = i4_uniform_ab ( 0, m - 1, seed );
    k = i4_uniform_ab ( 0, 1, seed );
    if ( k == 0 )
    {
      x[i+j*m] = 1.0;
    }
    else
    {
      x[i+j*m] = -1.0;
    }
  }
//
//  Shift by C and scale by R.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = c[i] + r * x[i+j*m];
    }
  }

  return x;
}
//****************************************************************************80

double *uniform_on_cube01 ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_CUBE01 returns random points on the surface of the unit cube.
//
//  Discussion:
//
//    The cube is assumed to be aligned with the coordinate axes.
//
//    The cube has center at the origin and radius 1. Any point on the surface
//    of the cube is described by
//
//      X = PM
//
//    where PM is an M-dimensional vector whose entries are between
//    -1 and +1, and for which at least one value has norm 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 April 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    1 <= M.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double UNIFORM_ON_CUBE[M*N], the coordinates of N points, chosen
//    uniformly at random from the surface of the unit M-cube.
//
{
  int i;
  int j;
  int k;
  double *x;
//
//  Choose random points within the cube of radius 1.
//
  x = r8mat_uniform_01_new ( m, n, seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = 2.0 * x[i+j*m] - 1.0;
    }
  }
//
//  For each point, select a coordinate at random, and set it to +1 or -1.
//
  for ( j = 0; j < n; j++ )
  {
    i = i4_uniform_ab ( 0, m - 1, seed );
    k = i4_uniform_ab ( 0, 1, seed );
    if ( k == 0 )
    {
      x[i+j*m] = 1.0;
    }
    else
    {
      x[i+j*m] = -1.0;
    }
  }

  return x;
}
//****************************************************************************80

double *uniform_on_ellipsoid_map ( int dim_num, int n, double a[],
  double r, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_ELLIPSOID_MAP maps uniform points onto an ellipsoid.
//
//  Discussion:
//
//    The points X on the ellipsoid are described by a DIM_NUM by DIM_NUM
//    positive definite symmetric matrix A, and a "radius" R, such that
//
//      X' * A * X = R * R
//
//    The algorithm computes the Cholesky factorization of A:
//
//      A = U' * U.
//
//    A set of uniformly random points Y is generated, satisfying:
//
//      Y' * Y = R * R.
//
//    The appropriate points in the ellipsoid are found by solving
//
//      U * X = Y
//
//    Thanks to Dr Karl-Heinz Keil for pointing out that the original
//    coding was actually correct only if A was replaced by its inverse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
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
//    Input, double A[DIM_NUM*DIM_NUM], the matrix that describes the ellipsoid.
//
//    Input, double R, the right hand side of the ellipsoid equation.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_ELLIPSOID_MAP[DIM_NUM*N], the points.
//
{
  int i;
  int info;
  int j;
  int k;
  double *u;
  double *x;
//
//  Get the factor U.
//
  u = new double[dim_num*dim_num];

  for ( j = 0; j < dim_num; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      u[i+j*dim_num] = a[i+j*dim_num];
    }
  }

  info = r8po_fa ( a, dim_num, dim_num );

  if ( info != 0 )
  {
    cerr << "\n";
    cerr << "UNIFORM_ON_ELLIPSOID_MAP - Fatal error!\n";
    cerr << "  R8PO_FA reports that the matrix A \n";
    cerr << "  is not positive definite symmetric.\n";
    exit ( 1 );
  }
//
//  Get the points Y that satisfy Y' * Y = R * R.
//
  x = uniform_on_sphere01_map ( dim_num, n, seed );

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = r * x[i+j*dim_num];
    }
  }
//
//  Solve U * X = Y.
//
  for ( j = 0; j < n; j++ )
  {
    r8po_sl ( u, dim_num, dim_num, x+j*dim_num );
  }

  delete [] u;

  return x;
}
//****************************************************************************80

double *uniform_on_hemisphere01_phong ( int n, int m, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_HEMISPHERE01_PHONG maps uniform points onto the unit hemisphere.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    The Phong density is used, with exponent M:
//
//    rho ( theta, phi; m ) = ( m + 1 ) * cos ( phi )**M / ( 2 * pi )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, int M, the Phong exponent.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_HEMISPHERE01_PHONG[3*N], the points.
//
{
# define DIM_NUM 3
# define PI 3.141592653589793

  int j;
  double phi;
  double power;
  double theta;
  double *x;

  x = new double[DIM_NUM*n];

  power = 1.0 / ( double ) ( m + 1 );

  for ( j = 0; j < n; j++ )
  {
    phi = r8_uniform_01 ( seed );

    phi = acos ( pow ( 1.0 - phi, power ) );

    theta = r8_uniform_01 ( seed );

    theta = 2.0 * PI * theta;

    x[0+j*DIM_NUM] = cos ( theta ) * sin ( phi );
    x[1+j*DIM_NUM] = sin ( theta ) * sin ( phi );
    x[2+j*DIM_NUM] = cos ( phi );
  }

  return x;
# undef DIM_NUM
# undef PI
}
//****************************************************************************80

double *uniform_on_simplex01_map ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SIMPLEX01_MAP maps uniform points onto the unit simplex.
//
//  Discussion:
//
//    The surface of the unit DIM_NUM-dimensional simplex is the set of points
//    X(1:DIM_NUM) such that each X(I) is nonnegative,
//    every X(I) is no greater than 1, and
//
//    ( X(I) = 0 for some I, or sum ( X(1:DIM_NUM) ) = 1. )
//
//    In DIM_NUM dimensions, there are DIM_NUM sides, and one main face.
//    This code picks a point uniformly with respect to "area".
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
//  Reference:
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SIMPLEX01_MAP[DIM_NUM*N], the points.
//
{
  double area1;
  double area2;
  double *e;
  int i;
  int j;
  double r;
  double total;
  double u;
  double *x;
//
//  The construction begins by sampling DIM_NUM points from the
//  exponential distribution with parameter 1.
//
  e = new double[dim_num];
  x = new double[dim_num*n];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( dim_num, seed, e );

    for ( i = 0; i < dim_num; i++ )
    {
      e[i] = -log ( e[i] );
    }

    total = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      total = total + e[i];
    }
//
//  Based on their relative areas, choose a side of the simplex,
//  or the main face.
//
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = e[i] / total;
    }

    area1 = sqrt ( ( double ) dim_num );
    area2 = ( double ) dim_num;

    r = r8_uniform_01 ( seed );

    if ( area1 / ( area1 + area2 ) < r )
    {
      i = i4_uniform_ab ( 0, dim_num-1, seed );
      x[i+j*dim_num] = 0.0;
    }

  }
  delete [] e;

  return x;
}
//****************************************************************************80

double *uniform_on_sphere01_map ( int dim_num, int n, int &seed )

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
//    12 November 2010
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
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_MAP[DIM_NUM*N], the points.
//
{
  int i;
  int j;
  double norm;
  double *x;

  x = r8mat_normal_01_new ( dim_num, n, seed );

  for ( j = 0; j < n; j++ )
  {
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      norm = norm + x[i+j*dim_num] * x[i+j*dim_num];
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = x[i+j*dim_num] / norm;
    }

  }

  return x;
}
//****************************************************************************80

double *uniform_on_sphere01_patch_tp ( int n, double phi1, double phi2,
  double theta1, double theta2, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_PATCH_TP maps uniform points to a spherical TP patch.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    A sphere TP patch on the surface of the unit sphere contains those
//    points with radius R = 1 and angles (THETA,PHI) such that
//
//      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
//      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
//
//    mapped to
//
//      X = cos ( THETA ) * sin ( PHI )
//      Y = sin ( THETA ) * sin ( PHI )
//      Z =                 cos ( PHI )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double PHI1, PHI2, the latitudinal angle range.
//
//    Input, double THETA1, THETA2, the longitudinal angle range.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_PATCH_TP[2*N], the TP points.
//
{
  int j;
  double *tp;

  tp = r8mat_uniform_01_new ( 2, n, seed );

  for ( j = 0; j < n; j++ )
  {
    tp[0+j*2] = ( 1.0 - tp[0+j*2] ) * theta1
              +         tp[0+j*2]   * theta2;

    tp[1+j*2] = acos ( ( 1.0 - tp[1+j*2] ) * cos ( phi1 )
                     +         tp[1+j*2]   * cos ( phi2 ) );
  }

  return tp;
}
//****************************************************************************80

double *uniform_on_sphere01_patch_xyz ( int n, double phi1, double phi2,
  double theta1, double theta2, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform points to a spherical XYZ patch.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    A sphere XYZ patch on the surface of the unit sphere contains those
//    points with radius R = 1 and angles (THETA,PHI) such that
//
//      0.0 <= THETA1 <= THETA <= THETA2 <= 2 * PI
//      0.0 <= PHI1   <= PHI   <= PHI2   <=     PI
//
//    mapped to
//
//      X = cos ( THETA ) * sin ( PHI )
//      Y = sin ( THETA ) * sin ( PHI )
//      Z =                 cos ( PHI )
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
//  Reference:
//
//    Peter Shirley,
//    Nonuniform Random Point Sets Via Warping,
//    Graphics Gems, Volume III,
//    edited by David Kirk,
//    AP Professional, 1992,
//    ISBN: 0122861663,
//    LC: T385.G6973.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double PHI1, PHI2, the latitudinal angle range.
//
//    Input, double THETA1, THETA2, the longitudinal angle range.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_PATCH_XYZ[3*N], the points.
//
{
# define DIM_NUM 3

  int j;
  double phi;
  double theta;
  double *x;

  x = new double[DIM_NUM*n];

  for ( j = 0; j < n; j++ )
  {
    phi = r8_uniform_01 ( seed );

    phi = acos ( ( 1.0 - phi ) * cos ( phi1 )
               +         phi   * cos ( phi2 ) );

    theta = r8_uniform_01 ( seed );

    theta = ( 1.0 - theta ) * theta1
          +         theta   * theta2;

    x[0+j*DIM_NUM] = cos ( theta ) * sin ( phi );
    x[1+j*DIM_NUM] = sin ( theta ) * sin ( phi );
    x[2+j*DIM_NUM] = cos ( phi );
  }

  return x;
# undef DIM_NUM
}
//****************************************************************************80

double *uniform_on_sphere01_triangle_xyz ( int n, double v1[], double v2[],
  double v3[], int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_SPHERE01_TRIANGLE_XYZ: sample spherical triangle, XYZ coordinates.
//
//  Discussion:
//
//    The sphere has center 0 and radius 1.
//
//    A spherical triangle on the surface of the unit sphere contains those
//    points with radius R = 1, bounded by the vertices V1, V2, V3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Arvo,
//    Stratified sampling of spherical triangles,
//    Computer Graphics Proceedings, Annual Conference Series,
//    ACM SIGGRAPH '95, pages 437-438, 1995.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double V1[3], V2[3], V3[3], the XYZ coordinates of
//    the vertices of the spherical triangle.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double UNIFORM_ON_SPHERE01_TRIANGLE_XYZ[3*N], the XYZ
//    coordinates of the sample points.
//
{
  double a;
  double alpha;
  double area;
  double area_hat;
  double b;
  double beta;
  double c;
  double gamma;
  int i;
  int j;
  double q;
  double r;
  double s;
  double t;
  double temp;
  double u;
  double v;
  double *v31;
  double *v4;
  double *v42;
  double w;
  double *x;
  double xsi1;
  double xsi2;
  double z;
//
//  Compute the sides, angles, and area of the spherical triangle;
//  for now, we assume R = 1.
//
  r = 1.0;

  stri_vertices_to_sides ( r, v1, v2, v3, &a, &b, &c );

  stri_sides_to_angles ( r, a, b, c, &alpha, &beta, &gamma );

  area = stri_angles_to_area ( r, alpha, beta, gamma );

  x = new double[3*n];
  v31 = new double[3];
  v4 = new double[3];
  v42 = new double[3];

  for ( j = 0; j < n; j++ )
  {
//
//  Select the new area.
//
    xsi1 = r8_uniform_01 ( seed );
    area_hat = xsi1 * area;
//
//  Compute the sine and cosine of the angle phi.
//
    s = sin ( area_hat - alpha );
    t = cos ( area_hat - alpha );
//
//  Compute the pair that determines beta_hat.
//
    u = t - cos ( alpha );
    v = s + sin ( alpha ) * cos ( c );
//
//  Q is the cosine of the new edge length b_hat.
//
    q = ( ( v * t - u * s ) * cos ( alpha ) - v )
      / ( ( v * s + u * t ) * sin ( alpha ) );
//
//  V31 = normalized ( V3 - ( V3 dot V1 ) * V1 )
//
    w = r8vec_dot_product ( 3, v3, v1 );

    for ( i = 0; i < 3; i++ )
    {
      v31[i] = v3[i] - w * v1[i];
    }
    temp = r8vec_norm ( 3, v31 );
    for ( i = 0; i < 3; i++ )
    {
      v31[i] = v31[i] / temp;
    }
//
//  V4 is the third vertex of the subtriangle V1, V2, V4.
//
    for ( i = 0; i < 3; i++ )
    {
      v4[i] = q * v1[i] + sqrt ( 1.0 - q * q ) * v31[i];
    }
//
//  Select cos theta, which will sample along the edge from V2 to V4.
//
    xsi2 = r8_uniform_01 ( seed );
    z = 1.0 - xsi2 * ( 1.0 - r8vec_dot_product ( 3, v4, v2 ) );
//
//  V42 = normalized ( V4 - ( V4 dot V2 ) * V2 )
//
    w = r8vec_dot_product ( 3, v4, v2 );
    for ( i = 0; i < 3; i++ )
    {
      v42[i] = v4[i] - w * v2[i];
    }
    temp = r8vec_norm ( 3, v42 );
    for ( i = 0; i < 3; i++ )
    {
      v42[i] = v42[i] / temp;
    }
//
//  Construct the point.
//
    for ( i = 0; i < 3; i++ )
    {
      x[i+j*3] = z * v2[i] + sqrt ( 1.0 - z * z ) * v42[i];
    }
  }

  delete [] v31;
  delete [] v4;
  delete [] v42;

  return x;
}
//****************************************************************************80

double *uniform_on_triangle ( int n, double v[], int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_TRIANGLE maps uniform points onto the boundary of a triangle.
//
//  Discussion:
//
//    The triangle is defined by the three vertices V1, V2, V3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int  N, the number of points.
//
//    Input, double V[2*3], the vertices of the triangle.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double UNIFORM_ON_TRIANGLE[2*N], the points.
//
{
  int j;
  double l1;
  double l2;
  double l3;
  int m = 2;
  double r;
  double s;
  double t;
  double *x;

  l1 = sqrt ( pow ( v[0+1*2] - v[0+0*2], 2 )
            + pow ( v[1+1*2] - v[1+0*2], 2 ) );

  l2 = sqrt ( pow ( v[0+2*2] - v[0+1*2], 2 ) 
            + pow ( v[1+2*2] - v[1+1*2], 2 ) );

  l3 = sqrt ( pow ( v[0+0*2] - v[0+2*2], 2 ) 
            + pow ( v[1+0*2] - v[1+2*2], 2 ) );
  
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  R can be regarded as the distance of the point on the perimeter,
//  as measured from the origin, along the perimeter.
//
    r = ( l1 + l2 + l3 ) * r8_uniform_01 ( seed );
//
//  Case 1: between V1 and V2.
//
    if ( r <= l1 )
    {
      s = ( l1 - r ) / l1;
      t =        r   / l1;
      x[0+j*2] = s * v[0+0*2] + t * v[0+1*2];
      x[1+j*2] = s * v[1+0*2] + t * v[1+1*2];
    }
//
//  Case 2: between V2 and V3.
//
    else if ( r <= l1 + l2 )
    {
      s = ( l2 - r + l1 ) / l2;
      t = (      r - l1 ) / l2;
      x[0+j*2] = s * v[0+1*2] + t * v[0+2*2];
      x[1+j*2] = s * v[1+1*2] + t * v[1+2*2];
    }
//
//  Case 3: between V3 and V1.
//
    else
    {
      s = ( l3 - r + l1 + l2 ) / l3;
      t = (      r - l1 - l2 ) / l3;
      x[0+j*2] = s * v[0+2*2] + t * v[0+0*2];
      x[1+j*2] = s * v[1+2*2] + t * v[1+0*2];
    }
  }

  return x;
}
//****************************************************************************80

double *uniform_on_triangle01 ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_ON_TRIANGLE01 maps uniform points onto the unit triangle.
//
//  Discussion:
//
//    The unit triangle is defined by the three vertices (1,0), (0,1) and (0,0).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 May 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input/output, int SEED, a seed for the random
//    number generator.
//
//    Output, double X[2*N], the points.
//
{
  double a;
  double b;
  int j;
  int m = 2;
  double r;
  double s;
  double *x;

  s = sqrt ( 2.0 );

  a =   1.0       / ( 2.0 + s );
  b = ( 1.0 + s ) / ( 2.0 + s );
  
  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
//
//  R can be regarded as the distance of the point on the perimeter,
//  as measured from the origin, along the perimeter.
//
    r = ( 2.0 + s ) * r8_uniform_01 ( seed );
//
//  Case 1: between (0,0) and (1,0).
//
    if ( r <= a )
    {
      x[0+j*2] = 0.0;
      x[1+j*2] = r;
    }
//
//  Case 2: between (1,0) and (0,1).
//
    else if ( r <= b )
    {
      x[0+j*2] = 1.0 - ( r - a ) * sqrt ( 2.0 ) / 2.0;
      x[1+j*2] = 0.0 + ( r - a ) * sqrt ( 2.0 ) / 2.0;
    }
//
//  Case 3: between (0,1) and (0,0).
//
    else
    {
      x[0+j*2] = 0.0;
      x[1+j*2] = 1.0 - ( r - b );
    }
  }

  return x;
}
//****************************************************************************80

double *uniform_walk ( int dim_num, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    UNIFORM_WALK generates points on a uniform random walk.
//
//  Discussion:
//
//    The first point is at the origin.  Uniform random numbers are
//    generated to determine the direction of the next step, which
//    is always of length 1, and in coordinate direction.
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
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double UNIFORM_WALK[DIM_NUM*N], the points.
//
{
  double arg;
  double dir;
  int i;
  int j;
  double *x;

  x = new double[dim_num*n];

  j = 0;
  for ( i = 0; i < dim_num; i++ )
  {
    x[i+j*dim_num] = 0.0;
  }

  for ( j = 1; j < n; j++ )
  {
    dir = r8_uniform_01 ( seed );
    dir = ( double ) ( 2 * dim_num ) * ( dir - 0.5 );

    for ( i = 0; i < dim_num; i++ )
    {
      x[i+j*dim_num] = x[i+(j-1)*dim_num];
    }
    arg = fabs ( dir ) + 0.5;
    i = r8_nint ( arg );
    i = i4_min ( i, dim_num );
    i = i4_max ( i, 1 );
    i = i - 1;

    if ( dir < 0.0 )
    {
      x[i+j*dim_num] = x[i+j*dim_num] - 1.0;
    }
    else
    {
      x[i+j*dim_num] = x[i+j*dim_num] + 1.0;
    }

  }

  return x;
}
