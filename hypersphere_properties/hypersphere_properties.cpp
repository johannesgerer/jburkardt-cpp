# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "hypersphere_properties.hpp"

//****************************************************************************80

void cartesian_to_hypersphere ( int m, int n, double c[], double x[], 
  double r[], double theta[] )

//****************************************************************************80
//
//  Purpose:
//
//    CARTESIAN_TO_HYPERSPHERE: Cartesian to hypersphere coordinate transform.
//
//  Discussion:
//
//    We allow the trivial case M = 1; in that case alone, the value R
//    must be assumed to have a sign.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
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
//    Input, int N, the number of points to transform.
//
//    Input, double C[M], the center of the hypersphere.
//
//    Input, double X[M*N], the Cartesian coordinates of the points.
//
//    Output, double R[N], the radius of the points on the 
//    hypersphere.  Except for the trivial case M = 1, R is assumed nonnegative.
//
//    Output, double THETA[(M-1)*N], the coordinate angles of the 
//    points, measured in radians.
//
{
  int i;
  int i1;
  int j;
  double t;
  double *top;
  double *x2;
//
//  Handle special case of M = 1.
//
  if ( m == 1 )
  {
    for ( j = 0; j < n; j++ )
    {
      r[j] = x[0+j*m] - c[0];
    }
    return;
  }
//
//  Subtract the center.
//
  x2 = new double[ m * n ];
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x2[i+j*m] = x[i+j*m] - c[i];
    }
  }
//
//  Compute R.
//
  for ( j = 0; j < n; j++ )
  {
    t = 0.0;
    for ( i = 0; i < m; i++ )
    {
      t = t + pow ( x2[i+m*j], 2 );
    }
    r[j] = sqrt ( t );
  }
//
//  Compute M-2 components of THETA.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      theta[i+j*(m-1)] = 0.0;
    }
  }

  for ( i = 1; i < m - 1; i++ )
  {
    for ( i1 = 0; i1 <= i - 1; i1++ )
    {
      for ( j = 0; j < n; j++ )
      {
        theta[i1+j*(m-1)] = theta[i1+j*(m-1)] + pow ( x2[i+j*m], 2 );
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 2; i++ )
    {
      theta[i+j*(m-1)] = theta[i+j*(m-1)] + pow ( x2[m-1+j*m], 2 );
    }
  }

  for ( i = 0; i < m - 2; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      theta[i+j*(m-1)] = atan2 ( sqrt ( theta[i+j*(m-1)] ), x2[i+j*m] );
    }
  }
//
//  Compute last component of THETA.
//
  top = new double[n];

  for ( j = 0; j < n; j++ )
  {
    top[j] = sqrt ( pow ( x2[m-1+j*m], 2 ) + pow ( x2[m-2+j*m], 2 ) ) + x2[m-2+j*m];
  }

  for ( j = 0; j < n; j++ )
  {
    theta[m-2+j*(m-1)] = 2.0 * atan2 ( x2[m-1+j*m], top[j] );
  }

  delete [] top;
  delete [] x2;

  return;
}
//****************************************************************************80

double hypersphere_01_area ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_AREA computes the surface area of a unit hypersphere.
//
//  Discussion:
//
//    The unit hypersphere satisfies the equation:
//
//      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
//
//    M   Area
//
//     2    2        * PI
//     3    4        * PI
//     4  ( 2 /   1) * PI^2
//     5  ( 8 /   3) * PI^2
//     6  ( 1 /   1) * PI^3
//     7  (16 /  15) * PI^3
//     8  ( 1 /   3) * PI^4
//     9  (32 / 105) * PI^4
//    10  ( 1 /  12) * PI^5
//
//    For the unit hypersphere, Area(M) = M * Volume(M)
//
//    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Output, double HYPERSPHERE_01_AREA, the area.
//
{
  double area;
  int i;
  int m2;
  double pi = 3.141592653589793;

  if ( ( m % 2 ) == 0 )
  {
    m2 = m / 2;
    area = 2.0 * pow ( pi, m2 );
    for ( i = 1; i <= m2 - 1; i++ )
    {
      area = area / ( ( double ) i );
    }
  }
  else
  {
    m2 = ( m - 1 ) / 2;
    area = pow ( 2.0, m ) * pow ( pi, m2 );
    for ( i = m2 + 1; i <= 2 * m2; i++ )
    {
      area = area / ( ( double ) i );
    }
  }

  return area;
}
//****************************************************************************80

void hypersphere_01_area_values ( int &n_data, int &m, double &area )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_AREA_VALUES returns some areas of the unit hypersphere.
//
//  Discussion:
//
//    The unit hypersphere satisfies the equation:
//
//      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
//
//     M         Area
//
//     2    2        * PI
//     3  ( 4 /    ) * PI
//     4  ( 2 /   1) * PI^2
//     5  ( 8 /   3) * PI^2
//     6  ( 1 /   1) * PI^3
//     7  (16 /  15) * PI^3
//     8  ( 1 /   3) * PI^4
//     9  (32 / 105) * PI^4
//    10  ( 1 /  12) * PI^5
//
//    For the unit hypersphere, Area(M) = M * Volume(M)
//
//    Sphere_Unit_Area ( M ) = 2 * PI^(M/2) / Gamma ( M / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &M, the spatial dimension.
//
//    Output, double &AREA, the area.
//
{
# define NMAX 20

  double area_vec[NMAX] = {
    0.2000000000000000E+01, 
    0.6283185307179586E+01, 
    0.1256637061435917E+02, 
    0.1973920880217872E+02, 
    0.2631894506957162E+02, 
    0.3100627668029982E+02, 
    0.3307336179231981E+02, 
    0.3246969701133415E+02, 
    0.2968658012464836E+02, 
    0.2550164039877345E+02, 
    0.2072514267328890E+02, 
    0.1602315322625507E+02, 
    0.1183817381218268E+02, 
    0.8389703410491089E+01, 
    0.5721649212349567E+01, 
    0.3765290085742291E+01, 
    0.2396678817591364E+01, 
    0.1478625959000308E+01, 
    0.8858104195716824E+00, 
    0.5161378278002812E+00 };
  int m_vec[NMAX] = {
     1,
     2,
     3, 
     4,
     5, 
     6,
     7, 
     8,
     9, 
    10,
    11,
    12,
    13, 
    14,
    15, 
    16,
    17, 
    18,
    19, 
    20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( NMAX <= n_data )
  {
    n_data = 0;
    m = 0;
    area = 0.0;
  }
  else
  {
    m = m_vec[n_data];
    area = area_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef NMAX
}
//****************************************************************************80

double *hypersphere_01_interior_uniform ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_INTERIOR_UNIFORM: uniform points inside unit hypersphere.
//
//  Discussion:
//
//    The hypersphere has center 0 and radius 1.
//
//    This routine is valid for any spatial dimension.
//
//    We first generate a point ON the hypersphere, and then distribute it
//    IN the hypersphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
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
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double HYPERSPHERE_01_INTERIOR_UNIFORM[M*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *v;
  double *x;

  x = new double[m*n];

  exponent = 1.0 / ( double ) ( m );

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    v = r8vec_normal_01_new ( m, seed );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < m; i++ )
    {
      norm = norm + pow ( v[i], 2 );
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = v[i] / norm;
    }
//
//  Now compute a value to map the point ON the hypersphere INTO the hypersphere.
//
    r = r8_uniform_01 ( seed );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = pow ( r, exponent ) * x[i+j*m];
    }
    delete [] v;
  }

  return x;
}
//****************************************************************************80

double *hypersphere_01_surface_uniform ( int m, int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_SURFACE_UNIFORM: uniform points on unit hypersphere surface.
//
//  Discussion:
//
//    The hypersphere has center 0 and radius 1.
//
//    This procedure is valid for any spatial dimension DIM_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
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
//    George Marsaglia,
//    Choosing a point from the surface of a sphere,
//    Annals of Mathematical Statistics,
//    Volume 43, Number 2, April 1972, pages 645-646.
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
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double HYPERSPHERE_01_UNIFORM_SURFACE[M*N], the points.
//
{
  int i;
  int j;
  double norm;
  double *x;
//
//  Fill a matrix with normally distributed values.
//
  x = r8mat_normal_01_new ( m, n, seed );
//
//  Normalize each column.
//
  for ( j = 0; j < n; j++ )
  {
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < m; i++ )
    {
      norm = norm + pow ( x[i+j*m], 2 );
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+j*m] / norm;
    }
  }
  return x;
}
//****************************************************************************80

double hypersphere_01_volume ( int m )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_VOLUME computes the volume of a unit hypersphere.
//
//  Discussion:
//
//    The unit hypersphere satisfies the equation:
//
//      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
//
//     M  Volume
//
//     1    2
//     2    1        * PI
//     3  ( 4 /   3) * PI
//     4  ( 1 /   2) * PI^2
//     5  ( 8 /  15) * PI^2
//     6  ( 1 /   6) * PI^3
//     7  (16 / 105) * PI^3
//     8  ( 1 /  24) * PI^4
//     9  (32 / 945) * PI^4
//    10  ( 1 / 120) * PI^5
//
//    For the unit hypersphere, Volume(M) = 2 * PI * Volume(M-2)/ M
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Output, double HYPERSPHERE_01_VOLUME, the volume.
//
{
  int i;
  int m2;
  double pi = 3.141592653589793;
  double volume;

  if ( m % 2== 0 )
  {
    m2 = m / 2;
    volume = 1.0;
    for ( i = 1; i <= m2; i++ )
    {
      volume = volume * pi / ( ( double ) i );
    }
  }
  else
  {
    m2 = ( m - 1 ) / 2;
    volume = pow ( pi, m2 ) * pow ( 2.0, m );
    for ( i = m2 + 1; i <= 2 * m2 + 1; i++ )
    {
      volume = volume / ( ( double ) i );
    }
  }

  return volume;
}
//****************************************************************************80

void hypersphere_01_volume_values ( int &n_data, int &m, double &volume )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_01_VOLUME_VALUES returns some volumes of the unit hypersphere.
//
//  Discussion:
//
//    The unit hypersphere satisfies the equation:
//
//      Sum ( 1 <= I <= M ) X(I) * X(I) = 1
//
//     M  Volume
//
//     1    1
//     2    1        * PI
//     3  ( 4 /   3) * PI
//     4  ( 1 /   2) * PI^2
//     5  ( 8 /  15) * PI^2
//     6  ( 1 /   6) * PI^3
//     7  (16 / 105) * PI^3
//     8  ( 1 /  24) * PI^4
//     9  (32 / 945) * PI^4
//    10  ( 1 / 120) * PI^5
//
//    For the unit hypersphere, Volume(M) = 2 * PI * Volume(M-2) / M
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &M, the spatial dimension.
//
//    Output, double &VOLUME, the volume.
//
{
# define N_MAX 20

  int m_vec[N_MAX] = {
    1,
    2, 
    3,
    4, 
    5,
    6, 
    7,
    8, 
    9,
   10,
   11,
   12, 
   13,
   14, 
   15,
   16, 
   17,
   18, 
   19,
   20 };
  double volume_vec[N_MAX] = {
    0.2000000000000000E+01, 
    0.3141592653589793E+01, 
    0.4188790204786391E+01, 
    0.4934802200544679E+01, 
    0.5263789013914325E+01, 
    0.5167712780049970E+01, 
    0.4724765970331401E+01, 
    0.4058712126416768E+01, 
    0.3298508902738707E+01, 
    0.2550164039877345E+01, 
    0.1884103879389900E+01, 
    0.1335262768854589E+01, 
    0.9106287547832831E+00, 
    0.5992645293207921E+00, 
    0.3814432808233045E+00, 
    0.2353306303588932E+00, 
    0.1409811069171390E+00, 
    0.8214588661112823E-01, 
    0.4662160103008855E-01, 
    0.2580689139001406E-01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    m = 0;
    volume = 0.0;
  }
  else
  {
    m = m_vec[n_data];
    volume = volume_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

double hypersphere_area ( int m, double r )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_AREA computes the surface area of a hypersphere.
//
//  Discussion:
//
//    A hypersphere satisfies the equation:
//
//      sum ( ( P(1:M) - C(1:M) )^2 ) = R^2
//
//    M   Area
//
//    2      2       * PI   * R
//    3      4       * PI   * R^2
//    4      2       * PI^2 * R^3
//    5      (8/3)   * PI^2 * R^4
//    6                PI^3 * R^5
//    7      (16/15) * PI^3 * R^6
//
//    Sphere_Area ( M, R ) = 2 * PI^(M/2) * R^(M-1) / Gamma ( M / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, double R, the radius.
//
//    Output, double HYPERSPHERE_AREA, the area.
//
{
  double value;

  value = pow ( r, m - 1  ) * hypersphere_01_area ( m );

  return value;
}
//****************************************************************************80

double *hypersphere_stereograph ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_STEREOGRAPH: stereographic mapping of points on a hypersphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//   16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    M must be at least 2.
//
//    Input, int N, the number of points.
//
//    Input, double X[M*N], the points to be mapped.
//
//    Output, double HYPERSPHERE_STEREOGRAPH[M-1)*N], the stereographically 
//    mapped points.
//
{
  int i;
  int j;
  double *x2;

  x2 = new double[(m-1)*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      x2[i+j*(m-1)] = x[i+j*m];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      x2[i+j*(m-1)] = x2[i+j*(m-1)] / ( 1.0 - x[m-1+j*m] );
    }
  }
  
  return x2;
}
//****************************************************************************80

double *hypersphere_stereograph_inverse ( int m, int n, double x2[] )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_STEREOGRAPH_INVERSE inverts a stereographic map.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    M must be at least 2.
//
//    Input, int N, the number of points.
//
//    Input, double X2[(M-1)*N], points in the plane.
//
//    Input, double HYPERSPHERE_STEREOGRAPH_INVERSE[M*N], points mapped 
//    back to the hypersphere.
//
{
  double *d;
  int i;
  int j;
  double *x;

  x = new double[m*n];

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m - 1; i++ )
    {
      x[i+j*m] = 2.0 * x2[i+j*(m-1)];
    }
  }

  d = new double[n];

  for ( j = 0; j < n; j++ )
  {
    d[j] = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      d[j] = d[j] + pow ( x2[i+j*(m-1)], 2 );
    }
  }

  for ( j = 0; j < n; j++ )
  {
    x[m-1+j*m] = d[j] - 1.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+j*m] / ( d[j] + 1.0 );
    }
  }

  delete [] d;

  return x;
}
//****************************************************************************80

double *hypersphere_surface_uniform ( int m, int n, double r, double c[], 
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_SURFACE_UNIFORM: uniform hypersphere surface samples
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
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
//    George Marsaglia,
//    Choosing a point from the surface of a sphere,
//    Annals of Mathematical Statistics,
//    Volume 43, Number 2, April 1972, pages 645-646.
//
//    Reuven Rubinstein,
//    Monte Carlo Optimization, Simulation, and Sensitivity
//    of Queueing Networks,
//    Wiley, 1986, page 234.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input, real R, the radius.
//
//    Input, real C[M], the center.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, real HYPERSPHERE_SURFACE_UNIFORM[M*N], the points.
//
{
  int i;
  int j;
  double *x;

  x = hypersphere_01_surface_uniform ( m, n, seed );
//
//  Scale by the radius.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * x[i+j*m];
    }
  }
//
//  Shift to the center.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+j*m] + c[i];
    }
  }

  return x;
}
//****************************************************************************80

double *hypersphere_to_cartesian ( int m, int n, double c[], double r[], 
  double theta[] )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_TO_CARTESIAN: hypersphere to Cartesian coordinate transform.
//
//  Discussion:
//
//    We allow the trivial case M = 1; in that case alone, the value R
//    must be assumed to have a sign.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer M, the spatial dimension.
//    1 <= M.
//
//    Input, integer N, the number of points to transform.
//
//    Input, real C[M], the center of the hypersphere.
//
//    Input, real R[N], the radius of the points on the hypersphere.
//    Except for the trivial case M = 1, R is assumed nonnegative.
//
//    Input, real THETA[(M-1)*N], the coordinate angles of the points,
//    measured in radians.
//
//    Output, real HYPERSPHERE_TO_CARTESIAN[M*N], the Cartesian 
//    coordinates of the points.
//
{
  int i;
  int i1;
  int i2;
  int j;
  double *x;

  x = new double[m*n];

  if ( m == 1 )
  {
    for ( j = 0; j < n; j++ )
    {
      x[0+j*m] = r[j];
    }
  }
  else
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < m; i++ )
      {
        x[i+j*m] = r[j];
      }
    }

    for ( j = 0; j < n; j++ )
    {
      for ( i1 = 0; i1 < m - 1; i1++ )
      {
        x[i1+j*m] = x[i1+j*m] * cos ( theta[i1+j*(m-1)] );
        for ( i2 = i1 + 1; i2 < m; i2++ )
        {
          x[i2+j*m] = x[i2+j*m] * sin ( theta[i1+j*(m-1)] );
        }
      }
    }

  }
//
//  Add the center.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+j*m] + c[i];
    }
  }

  return x;
}
//****************************************************************************80

double hypersphere_volume ( int m, double r )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERSPHERE_VOLUME computes the volume of a hypersphere.
//
//  Discussion:
//
//    A hypersphere satisfies the equation:
//
//      sum ( ( X(1:M) - PC(1:M) )^2 ) = R^2
//
//    where R is the radius and PC is the center.
//
//    Results include:
//
//    M     Volume
//    -     -----------------------
//    2                PI   * R^2
//    3     (4/3)    * PI   * R^3
//    4     (1/2)    * PI^2 * R^4
//    5     (8/15)   * PI^2 * R^5
//    6     (1/6)    * PI^3 * R^6
//    7     (16/105) * PI^3 * R^7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    16 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, double R, the radius.
//
//    Output, double HYPERSPHERE_VOLUME, the volume.
//
{
  double

  value = pow ( r, m ) * hypersphere_01_volume ( m );

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

double r8mat_norm_fro_affine ( int m, int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_NORM_FRO_AFFINE returns the Frobenius norm of an R8MAT difference.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)^2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows.
//
//    Input, int N, the number of columns.
//
//    Input, double A1[M*N], A2[M,N], the matrice for whose difference the 
//    Frobenius norm is desired.
//
//    Output, double R8MAT_NORM_FRO_AFFINE, the Frobenius norm of A1 - A2.
//
{
  int i;
  int j;
  double value;

  value = 0.0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      value = value + pow ( a1[i+j*m] - a2[i+j*m], 2 );
    }
  }
  value = sqrt ( value );

  return value;
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.
//
{
  int i;
  int m;
  const double pi = 3.141592653589793;
  double *r;
  double *x;
  int x_hi;
  int x_lo;

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  If we need just one new value, do that here to avoid null arrays.
//
  if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );

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

    delete [] r;
  }

  return x;
}
//****************************************************************************80

void r8vec_transpose_print ( int n, double a[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_TRANSPOSE_PRINT prints an R8VEC "transposed".
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Example:
//
//    A = (/ 1.0, 2.1, 3.2, 4.3, 5.4, 6.5, 7.6, 8.7, 9.8, 10.9, 11.0 /)
//    TITLE = 'My vector:  '
//
//    My vector:
//        1.0    2.1    3.2    4.3    5.4
//        6.5    7.6    8.7    9.8   10.9
//       11.0
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
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A[N], the vector to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;
  int ihi;
  int ilo;

  cout << "\n";
  cout << title << "\n";

  if ( n <= 0 )
  {
    cout << "  (Empty)\n";
    return;
  }

  for ( ilo = 0; ilo < n; ilo = ilo + 5 )
  {
    ihi = i4_min ( ilo + 5, n );
    for ( i = ilo; i < ihi; i++ )
    {
      cout << "  " << setw(12) << a[i];
    }
    cout << "\n";
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
//    16 December 2013
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
//    16 December 2013
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

  p = new double[ m * n ];

  for ( j = 0; j < n; j++ )
  {
    qn = 0.0;
    for ( i = 0; i < m - 1; i++ )
    {
      qn = qn + pow ( q[i+j*(m-1)], 2 );
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
