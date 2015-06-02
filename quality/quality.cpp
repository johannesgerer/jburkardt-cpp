# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>

using namespace std;

# include "quality.hpp"

//****************************************************************************80

double alpha_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    ALPHA_MEASURE determines the triangulated pointset quality measure ALPHA.
//
//  Discusion:
//
//    The ALPHA measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    We compute the minimum angle among all the triangles in the triangulated
//    dataset and divide by the maximum possible value (which, in degrees,
//    is 60).  The best possible value is 1, and the worst 0.  A good
//    triangulation should have an ALPHA score close to 1.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, real ( kind = 8 ) Z(2,N), the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
//    the triangulation.
//
//    Output, double ALPHA_MEASURE, the ALPHA quality measure.
//
{
  double a_angle;
  int a_index;
  double a_x;
  double a_y;
  double ab_len;
  double alpha;
  double b_angle;
  int b_index;
  double b_x;
  double b_y;
  double bc_len;
  double c_angle;
  int c_index;
  double c_x;
  double c_y;
  double ca_len;
  double pi = 3.141592653589793;
  int triangle;
  double value;

  alpha = r8_huge ( );

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*3];
    b_index = triangle_node[1+triangle*3];
    c_index = triangle_node[2+triangle*3];

    a_x = z[0+(a_index-1)*2];
    a_y = z[1+(a_index-1)*2];
    b_x = z[0+(b_index-1)*2];
    b_y = z[1+(b_index-1)*2];
    c_x = z[0+(c_index-1)*2];
    c_y = z[1+(c_index-1)*2];

    ab_len = sqrt ( pow ( a_x - b_x, 2 ) + pow ( a_y - b_y, 2 ) );
    bc_len = sqrt ( pow ( b_x - c_x, 2 ) + pow ( b_y - c_y, 2 ) );
    ca_len = sqrt ( pow ( c_x - a_x, 2 ) + pow ( c_y - a_y, 2 ) );
//
//  Take care of a ridiculous special case.
//
    if ( ab_len == 0.0 && bc_len == 0.0 && ca_len == 0.0 )
    {
      a_angle = 2.0 * pi / 3.0;
      b_angle = 2.0 * pi / 3.0;
      c_angle = 2.0 * pi / 3.0;
    }
    else
    {
      if ( ca_len == 0.0 || ab_len == 0.0 )
      {
        a_angle = pi;
      }
      else
      {
        a_angle = arc_cosine (
          ( ca_len * ca_len + ab_len * ab_len - bc_len * bc_len )
          / ( 2.0 * ca_len * ab_len ) );
      }

      if ( ab_len == 0.0 || bc_len == 0.0 )
      {
        b_angle = pi;
      }
      else
      {
        b_angle = arc_cosine (
          ( ab_len * ab_len + bc_len * bc_len - ca_len * ca_len )
          / ( 2.0 * ab_len * bc_len ) );
      }

      if ( bc_len == 0.0 || ca_len == 0.0 )
      {
        c_angle = pi;
      }
      else
      {
        c_angle = arc_cosine (
          ( bc_len * bc_len + ca_len * ca_len - ab_len * ab_len )
          / ( 2.0 * bc_len * ca_len ) );
      }
    }
    alpha = r8_min ( alpha, a_angle );
    alpha = r8_min ( alpha, b_angle );
    alpha = r8_min ( alpha, c_angle );
  }
//
//  Normalize angles from [0,60] degrees into qualities in [0,1].
//
  value = alpha * 3.0 / pi;

  return value;
}
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
  const double r8_pi = 3.141592653589793;
  double value;

  if ( c <= -1.0 )
  {
    value = r8_pi;
  }
  else if ( 1.0 <= c )
  {
    value = 0.0;
  }
  else
  {
    value = acos ( c );
  }

  return value;
}
//****************************************************************************80

double area_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    AREA_MEASURE determines the area ratio quality measure.
//
//  Discusion:
//
//    This measure computes the area of every triangle, and returns
//    the ratio of the minimum to the maximum triangle.  A value of
//    1 is "perfect", indicating that all triangles have the same area.
//    A value of 0 is the worst possible result.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double AREA_MEASURE, the AREA quality measure.
//
{
  double area;
  double area_max;
  double area_min;
  int triangle;
  double value;
  double x1;
  double x2;
  double x3;
  double y1;
  double y2;
  double y3;

  area_max = 0.0;
  area_min = r8_huge ( );

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    x1 = z[0+(triangle_node[0+triangle*3]-1)*2];
    y1 = z[1+(triangle_node[0+triangle*3]-1)*2];
    x2 = z[0+(triangle_node[1+triangle*3]-1)*2];
    y2 = z[1+(triangle_node[1+triangle*3]-1)*2];
    x3 = z[0+(triangle_node[2+triangle*3]-1)*2];
    y3 = z[1+(triangle_node[2+triangle*3]-1)*2];

    area = 0.5 * fabs ( x1 * ( y2 - y3 )
                      + x2 * ( y3 - y1 )
                      + x3 * ( y1 - y2 ) );

    area_min = r8_min ( area_min, area );
    area_max = r8_max ( area_max, area );
  }
  if ( 0.0 < area_max )
  {
    value = area_min / area_max;
  }
  else
  {
    value = 0.0;
  }

  return value;
}
//****************************************************************************80

void bandwidth_mesh ( int element_order, int element_num, int element_node[],
  int *ml, int *mu, int *m )

//****************************************************************************80
//
//  Purpose:
//
//    BANDWIDTH_MESH determines the bandwidth of the coefficient matrix.
//
//  Discussion:
//
//    The quantity computed here is the "geometric" bandwidth determined
//    by the finite element mesh alone.
//
//    If a single finite element variable is associated with each node
//    of the mesh, and if the nodes and variables are numbered in the
//    same way, then the geometric bandwidth is the same as the bandwidth
//    of a typical finite element matrix.
//
//    The bandwidth M is defined in terms of the lower and upper bandwidths:
//
//      M = ML + 1 + MU
//
//    where
//
//      ML = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but earlier column,
//
//      MU = maximum distance from any diagonal entry to a nonzero
//      entry in the same row, but later column.
//
//    Because the finite element node adjacency relationship is symmetric,
//    we are guaranteed that ML = MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input,  ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM];
//    ELEMENT_NODE(I,J) is the global index of local node I in element J.
//
//    Output, int *ML, *MU, the lower and upper bandwidths of the matrix.
//
//    Output, int *M, the bandwidth of the matrix.
//
{
  int element;
  int global_i;
  int global_j;
  int local_i;
  int local_j;

  *ml = 0;
  *mu = 0;

  for ( element = 0; element < element_num; element++ )
  {
    for ( local_i = 0; local_i < element_order; local_i++ )
    {
      global_i = element_node[local_i+element*element_order];

      for ( local_j = 0; local_j < element_order; local_j++ )
      {
        global_j = element_node[local_j+element*element_order];

        *mu = i4_max ( *mu, global_j - global_i );
        *ml = i4_max ( *ml, global_i - global_j );
      }
    }
  }

  *m = *ml + 1 + *mu;

  return;
}
//****************************************************************************80

double beta_measure ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_MEASURE determines the pointset quality measure BETA.
//
//  Discussion:
//
//    The BETA measure of point distribution quality for a set Z of
//    N points in an DIM_NUM dimensional region is defined as follows:
//
//    For each point Z(I), determine the nearest distinct element of
//    the pointset by
//
//      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
//
//    Let GAMMA_AVE be the average of GAMMA(1:N).
//
//    Let GAMMA_STD be the standard deviation of the GAMMA's:
//
//      GAMMA_STD = sqrt ( 1 / ( N - 1 )
//        * sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 ) )
//
//    Then BETA is the standard deviation normalized by the average:
//
//      BETA = GAMMA_STD / GAMMA_AVE.
//
//    For an ideally regular mesh, the GAMMA(I)'s will be equal and
//    BETA will be zero.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Output, double BETA_MEASURE, the BETA quality measure.
//
{
  double *gamma;
  double gamma_ave;
  double gamma_std;
  int i;
  double value;

  gamma = pointset_spacing ( dim_num, n, z );

  gamma_ave = 0.0;
  for ( i = 0; i < n; i++ )
  {
    gamma_ave = gamma_ave + gamma[i];
  }
  gamma_ave = gamma_ave / ( double ) ( n );

  if ( 1 < n )
  {
    gamma_std = 0.0;
    for ( i = 0; i < n; i++ )
    {
      gamma_std = gamma_std + pow ( gamma[i] - gamma_ave, 2 );
    }
    gamma_std = sqrt ( gamma_std / ( double ) ( n - 1 ) );
  }
  else
  {
    gamma_std = 0.0;
  }

  if ( 0.0 < gamma_ave )
  {
    value = gamma_std / gamma_ave;
  }
  else
  {
    value = 0.0;
  }

  delete [] gamma;

  return value;
}
//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 )
  {
    c = c - 32;
  }

  return c;
}
//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 )
  {
    c1 = c1 - 32;
  }
  if ( 97 <= c2 && c2 <= 122 )
  {
    c2 = c2 - 32;
  }

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

double chi_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_MEASURE determines the pointset quality measure CHI.
//
//  Discussion:
//
//    The CHI measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//    Assign every point X in the region to the nearest element
//    Z(I) of the point set.  For each Z(I), let H(I) be the maximum
//    distance between Z(I) and any point assigned to it by this process.
//
//    For each point Z(I), we determine the nearest distinct element of
//    the pointset by
//
//      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
//
//    Then
//
//      CHI(I) = 2 * H(I) / GAMMA(I)
//
//    and
//
//      CHI = maximum ( 1 <= I <= N ) CHI(I)
//
//    This quantity can be estimated by using sampling to pick a large
//    number of points in the region, rather than all of them.
//
//    For an ideally regular mesh, all the CHI(I)'s will be equal.
//    Any deviation from regularity increases the value of some entries
//    of CHI; thus, given two meshes, the one with a lower value of
//    CHI is to be recommended.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double CHI_MEASURE, the CHI quality measure.
//
{
  double chi;
  double *chi_vec;
  int closest[1];
  double dist;
  double *gamma;
  double *h;
  int i;
  int j;
  int k;
  int seed;
  double *x;

  seed = seed_init;

  chi_vec = new double[n];
  h = new double[n];

  for ( j = 0; j < n; j++ )
  {
    h[j] = 0.0;
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*(dim_num)], 2 );
    }
    h[closest[0]] = r8_max ( h[closest[0]], dist );

    delete [] x;
  }

  gamma = pointset_spacing ( dim_num, n, z );

  chi = 0.0;

  for ( j = 0; j < n; j++ )
  {
    chi_vec[j] = 2.0 * sqrt ( h[j] ) / gamma[j];
    chi = r8_max ( chi, chi_vec[j] );
  }

  delete [] chi_vec;
  delete [] gamma;
  delete [] h;

  return chi;
}
//****************************************************************************80

double d_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    D_MEASURE determines the pointset quality measure D.
//
//  Discussion:
//
//    The D measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//    For each point Z(I) in the pointset, let V(I) be the subregion
//    defined by the intersection of the region with the Voronoi
//    region associated with Z(I).
//
//    Let D(I) be the determinant of the deviatoric tensor associated with
//    the region V(I).
//
//    Then D = maximum ( 1 <= I <= N ) D(I).
//
//    This quantity can be estimated using sampling.  A given number of
//    sample points are generated in the region, assigned to the nearest
//    element of the pointset, and used to approximate the Voronoi regions
//    and the deviatoric tensors.
//
//    In an ideally regular mesh, each deviatoric tensor would have a
//    zero determinant, and hence D would be zero.  In general, the smaller
//    D, the better.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double D_MEASURE, the D quality measure.
//
{
  double *a;
  double *centroid;
  int closest[1];
  double d;
  double di;
  int *hit;
  int i;
  int i1;
  int i2;
  int j;
  int k;
  double *moment;
  int seed;
  double *tri;
  double *x;

  a = new double[dim_num*dim_num];
  centroid = new double[dim_num*n];
  hit = new int[n];
  moment = new double[dim_num*dim_num*n];
  tri = new double[n];

  seed = seed_init;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      centroid[i+j*dim_num] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    hit[j] = 0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i2 = 0; i2 < dim_num; i2++ )
    {
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        moment[i1+i2*dim_num+j*dim_num*dim_num] = 0.0;
      }
    }
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    hit[closest[0]] = hit[closest[0]] + 1;

    for ( i = 0; i < dim_num; i++ )
    {
      centroid[i+closest[0]*dim_num] = centroid[i+closest[0]*dim_num] + x[i];
    }

    for ( i1 = 0; i1 < dim_num; i1++ )
    {
      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        moment[i1+i2*dim_num+closest[0]*dim_num*dim_num]
        = moment[i1+i2*dim_num+closest[0]*dim_num*dim_num] + x[i1] * x[i2];
      }
    }
    delete [] x;
  }

  for ( j = 0; j < n; j++ )
  {
    if ( 0 < hit[j] )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        centroid[i+j*dim_num] = centroid[i+j*dim_num] / ( double ) ( hit[j] );
      }
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        for ( i2 = 0; i2 < dim_num; i2++ )
        {
          moment[i1+i2*dim_num+j*dim_num*dim_num]
            = moment[i1+i2*dim_num+j*dim_num*dim_num] / ( double ) ( hit[j] );
        }
      }
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        for ( i2 = 0; i2 < dim_num; i2++ )
        {
          moment[i1+i2*dim_num+j*dim_num*dim_num] =
          moment[i1+i2*dim_num+j*dim_num*dim_num]
            - centroid[i1+j*dim_num] * centroid[i2+j*dim_num];
        }
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    tri[j] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      tri[j] = tri[j] + moment[i+i*dim_num+j*dim_num*dim_num];
    }
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      moment[i+i*dim_num+j*dim_num*dim_num] =
      moment[i+i*dim_num+j*dim_num*dim_num]
        - tri[j] / ( double ) ( dim_num );
    }
  }

  d = 0.0;

  for ( j = 0; j < n; j++ )
  {
    for ( i2 = 0; i2 < dim_num; i2++ )
    {
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        a[i1+i2*dim_num] = moment[i1+i2*dim_num+j*dim_num*dim_num];
      }
    }
    di = dge_det ( dim_num, a );

    d = r8_max ( d, di );

  }

  delete [] a;
  delete [] centroid;
  delete [] hit;
  delete [] moment;
  delete [] tri;

  return d;
}
//****************************************************************************80

double dge_det ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    DGE_DET computes the determinant of a square matrix in DGE storage.
//
//  Discussion:
//
//    The DGE storage format is used for a general M by N matrix.  A storage
//    space is made for each logical entry.  The two dimensional logical
//    array is mapped to a vector, in which storage is by columns.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Dongarra, Bunch, Moler, Stewart,
//    LINPACK User's Guide,
//    SIAM, 1979
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//    N must be positive.
//
//    Input/output, double A[N*N], the matrix to be analyzed.
//    On output, the matrix has been overwritten by factorization information.
//
//    Output, double DGE_DET, the determinant of the matrix.
//
{
  double det;
  int i;
  int j;
  int k;
  int l;
  double t;

  det = 1.0;

  for ( k = 1; k <= n-1; k++ )
  {
//
//  Find L, the index of the pivot row.
//
    l = k;
    for ( i = k+1; i <= n; i++ )
    {
      if ( fabs ( a[(l-1)+(k-1)*n] ) < fabs ( a[(i-1)+(k-1)*n] ) )
      {
        l = i;
      }
    }

    det = det * a[(l-1)+(k-1)*n];

    if ( a[(l-1)+(k-1)*n] == 0.0 )
    {
      return det;
    }
//
//  Interchange rows L and K if necessary.
//
    if ( l != k )
    {
      t                = a[(l-1)+(k-1)*n];
      a[(l-1)+(k-1)*n] = a[(k-1)+(k-1)*n];
      a[(k-1)+(k-1)*n] = t;
    }
//
//  Normalize the values that lie below the pivot entry A(K,K).
//
    for ( i = k+1; i <= n; i++ )
    {
      a[(i-1)+(k-1)*n] = -a[(i-1)+(k-1)*n] / a[(k-1)+(k-1)*n];
    }
//
//  Row elimination with column indexing.
//
    for ( j = k+1; j <= n; j++ )
    {
      if ( l != k )
      {
        t                = a[(l-1)+(j-1)*n];
        a[(l-1)+(j-1)*n] = a[(k-1)+(j-1)*n];
        a[(k-1)+(j-1)*n] = t;
      }

      for ( i = k+1; i <= n; i++ )
      {
        a[(i-1)+(j-1)*n] = a[(i-1)+(j-1)*n]
          + a[(i-1)+(k-1)*n] * a[(k-1)+(j-1)*n];
      }
    }
  }

  det = det * a[(n-1)+(n-1)*n];

  return det;
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
//    Original FORTRAN77 version by Barry Joe,
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

double *dtable_data_read ( char *input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_DATA_READ reads the data from a real TABLE file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double DTABLE_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "DTABLE_DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80

void dtable_header_read ( char *input_filename, int *m, int *n )

//****************************************************************************80
//
//  Purpose:
//
//    DTABLE_HEADER_READ reads the header from a real TABLE file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int *M, the number of spatial dimensions.
//
//    Output, int *N, the number of points
//
{
  *m = file_column_count ( input_filename );

  if ( *m <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_COLUMN_COUNT failed.\n";
    *n = -1;
    return;
  }

  *n = file_row_count ( input_filename );

  if ( *n <= 0 )
  {
    cout << "\n";
    cout << "DTABLE_HEADER_READ - Fatal error!\n";
    cout << "  FILE_ROW_COUNT failed.\n";
    return;
  }

  return;
}
//****************************************************************************80

int dtris2 ( int point_num, double point_xy[], int *tri_num,
  int tri_vert[], int tri_nabe[] )

//****************************************************************************80
//
//  Purpose:
//
//    DTRIS2 constructs a Delaunay triangulation of 2D vertices.
//
//  Discussion:
//
//    The routine constructs the Delaunay triangulation of a set of 2D vertices
//    using an incremental approach and diagonal edge swaps.  Vertices are
//    first sorted in lexicographically increasing (X,Y) order, and
//    then are inserted one at a time from outside the convex hull.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 January 2004
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe,
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
//    Input, int POINT_NUM, the number of vertices.
//
//    Input/output, double POINT_XY[POINT_NUM*2], the coordinates of the vertices.
//    On output, the vertices have been sorted into dictionary order.
//
//    Output, int *TRI_NUM, the number of triangles in the triangulation;
//    TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the number
//    of boundary vertices.
//
//    Output, int TRI_VERT[TRI_NUM*3], the nodes that make up each triangle.
//    The elements are indices of POINT_XY.  The vertices of the triangles are
//    in counter clockwise order.
//
//    Output, int TRI_NABE[TRI_NUM*3], the triangle neighbor list.
//    Positive elements are indices of TIL; negative elements are used for links
//    of a counter clockwise linked list of boundary edges; LINK = -(3*I + J-1)
//    where I, J = triangle, edge index; TRI_NABE[I,J] refers to
//    the neighbor along edge from vertex J to J+1 (mod 3).
//
//    Output, int RTRIS, is 0 for no error.
{
  double cmax;
  int e;
  int error;
  int i;
  int *indx;
  int j;
  int k;
  int l;
  int ledg;
  int lr;
  int ltri;
  int m;
  int m1;
  int m2;
  int n;
  int redg;
  int rtri;
  int *stack;
  int t;
  double tol;
  int top;
//
  stack = new int[point_num];

  tol = 100.0 * r8_epsilon ( );
//
//  Sort the vertices by increasing (x,y).
//
  indx = r82vec_sort_heap_index_a ( point_num, point_xy );

  r82vec_permute ( point_num, point_xy, indx );
//
//  Make sure that the data points are "reasonably" distinct.
//
  m1 = 1;

  for ( i = 2; i <= point_num; i++ )
  {
    m = m1;
    m1 = i;

    k = -1;

    for ( j = 0; j <= 1; j++ )
    {
      cmax = r8_max ( fabs ( point_xy[2*(m-1)+j] ),
                     fabs ( point_xy[2*(m1-1)+j] ) );

      if ( tol * ( cmax + 1.0 )
           < fabs ( point_xy[2*(m-1)+j] - point_xy[2*(m1-1)+j] ) )
      {
        k = j;
        break;
      }

    }

    if ( k == -1 )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      cout << "  Fails for point number I = " << i << "\n";
      cout << "  M =  " << m  << "\n";
      cout << "  M1 = " << m1 << "\n";
      cout << "  X,Y(M)  = " << point_xy[2*(m-1)+0] << "  "
                             << point_xy[2*(m-1)+1] << "\n";
      cout << "  X,Y(M1) = " << point_xy[2*(m1-1)+0] << "  "
                             << point_xy[2*(m1-1)+1] << "\n";
      delete [] stack;
      return 224;
    }
  }
//
//  Starting from points M1 and M2, search for a third point M that
//  makes a "healthy" triangle (M1,M2,M)
//
  m1 = 1;
  m2 = 2;
  j = 3;

  for ( ; ; )
  {
    if ( point_num < j )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      delete [] stack;
      return 225;
    }

    m = j;

    lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1],
      point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1],
      point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

    if ( lr != 0 )
    {
      break;
    }

    j = j + 1;

  }
//
//  Set up the triangle information for (M1,M2,M), and for any other
//  triangles you created because points were collinear with M1, M2.
//
  *tri_num = j - 2;

  if ( lr == -1 )
  {
    tri_vert[3*0+0] = m1;
    tri_vert[3*0+1] = m2;
    tri_vert[3*0+2] = m;
    tri_nabe[3*0+2] = -3;

    for ( i = 2; i <= *tri_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      tri_vert[3*(i-1)+0] = m1;
      tri_vert[3*(i-1)+1] = m2;
      tri_vert[3*(i-1)+2] = m;
      tri_nabe[3*(i-1)+0] = -3 * i;
      tri_nabe[3*(i-1)+1] = i;
      tri_nabe[3*(i-1)+2] = i - 1;

    }

    tri_nabe[3*(*tri_num-1)+0] = -3 * (*tri_num) - 1;
    tri_nabe[3*(*tri_num-1)+1] = -5;
    ledg = 2;
    ltri = *tri_num;
  }
  else
  {
    tri_vert[3*0+0] = m2;
    tri_vert[3*0+1] = m1;
    tri_vert[3*0+2] = m;
    tri_nabe[3*0+0] = -4;

    for ( i = 2; i <= *tri_num; i++ )
    {
      m1 = m2;
      m2 = i+1;
      tri_vert[3*(i-1)+0] = m2;
      tri_vert[3*(i-1)+1] = m1;
      tri_vert[3*(i-1)+2] = m;
      tri_nabe[3*(i-2)+2] = i;
      tri_nabe[3*(i-1)+0] = -3 * i - 3;
      tri_nabe[3*(i-1)+1] = i - 1;
    }

    tri_nabe[3*(*tri_num-1)+2] = -3 * (*tri_num);
    tri_nabe[3*0+1] = -3 * (*tri_num) - 2;
    ledg = 2;
    ltri = 1;
  }
//
//  Insert the vertices one at a time from outside the convex hull,
//  determine visible boundary edges, and apply diagonal edge swaps until
//  Delaunay triangulation of vertices (so far) is obtained.
//
  top = 0;

  for ( i = j+1; i <= point_num; i++ )
  {
    m = i;
    m1 = tri_vert[3*(ltri-1)+ledg-1];

    if ( ledg <= 2 )
    {
      m2 = tri_vert[3*(ltri-1)+ledg];
    }
    else
    {
      m2 = tri_vert[3*(ltri-1)+0];
    }

    lr = lrline ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1],
      point_xy[2*(m1-1)+0], point_xy[2*(m1-1)+1],
      point_xy[2*(m2-1)+0], point_xy[2*(m2-1)+1], 0.0 );

    if ( 0 < lr )
    {
      rtri = ltri;
      redg = ledg;
      ltri = 0;
    }
    else
    {
      l = -tri_nabe[3*(ltri-1)+ledg-1];
      rtri = l / 3;
      redg = (l % 3) + 1;
    }

    vbedg ( point_xy[2*(m-1)+0], point_xy[2*(m-1)+1], point_num,
      point_xy, *tri_num, tri_vert, tri_nabe, &ltri, &ledg, &rtri, &redg );

    n = *tri_num + 1;
    l = -tri_nabe[3*(ltri-1)+ledg-1];

    for ( ; ; )
    {
      t = l / 3;
      e = ( l % 3 ) + 1;
      l = -tri_nabe[3*(t-1)+e-1];
      m2 = tri_vert[3*(t-1)+e-1];

      if ( e <= 2 )
      {
        m1 = tri_vert[3*(t-1)+e];
      }
      else
      {
        m1 = tri_vert[3*(t-1)+0];
      }

      *tri_num = *tri_num + 1;
      tri_nabe[3*(t-1)+e-1] = *tri_num;
      tri_vert[3*(*tri_num-1)+0] = m1;
      tri_vert[3*(*tri_num-1)+1] = m2;
      tri_vert[3*(*tri_num-1)+2] = m;
      tri_nabe[3*(*tri_num-1)+0] = t;
      tri_nabe[3*(*tri_num-1)+1] = *tri_num - 1;
      tri_nabe[3*(*tri_num-1)+2] = *tri_num + 1;
      top = top + 1;

      if ( point_num < top )
      {
        cout << "\n";
        cout << "DTRIS2 - Fatal error!\n";
        cout << "  Stack overflow.\n";
        delete [] stack;
        return 8;
      }

      stack[top-1] = *tri_num;

      if ( t == rtri && e == redg )
      {
        break;
      }

    }

    tri_nabe[3*(ltri-1)+ledg-1] = -3 * n - 1;
    tri_nabe[3*(n-1)+1] = -3 * (*tri_num) - 2;
    tri_nabe[3*(*tri_num-1)+2] = -l;
    ltri = n;
    ledg = 2;

    error = swapec ( m, &top, &ltri, &ledg, point_num, point_xy, *tri_num,
      tri_vert, tri_nabe, stack );

    if ( error != 0 )
    {
      cout << "\n";
      cout << "DTRIS2 - Fatal error!\n";
      cout << "  Error return from SWAPEC.\n";
      delete [] stack;
      return error;
    }

  }
//
//  Now account for the sorting that we did.
//
  for ( i = 0; i < 3; i++ )
  {
    for ( j = 0; j < *tri_num; j++ )
    {
      tri_vert[i+j*3] = indx [ tri_vert[i+j*3] - 1 ];
    }
  }
  perm_inv ( point_num, indx );

  r82vec_permute ( point_num, point_xy, indx );

  delete [] indx;
  delete [] stack;

  return 0;
}
//****************************************************************************80

double e_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    E_MEASURE determines the pointset quality measure E.
//
//  Discussion:
//
//    The E measure of point distribution quality for a set Z of
//    N points in an DIM_NUM dimensional region is defined as follows:
//
//    Assign every point X in the region to the nearest element
//    Z(I) of the point set.  For each point Z(I), let E_VEC(I) be the
//    integral of the distance between Z(I) and all the points assigned to
//    it:
//
//      E_VEC(I) = Integral ( all X nearest to Z(I) ) distance ( X, Z(I) )
//
//    If we let VOLUME be the volume of the region, then we define E by:
//
//      E = sum ( 1 <= I <= N ) E_VEC(I) / VOLUME
//
//    This quantity can be estimated by using sampling to pick a large
//    number of points in the region, rather than all of them.
//
//    The E measure is minimized by a centroidal Voronoi tessellation.
//
//    Given two meshes, the one with a lower value of E is to be recommended.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double E_MEASURE, the E quality measure.
//
{
  int closest[1];
  double dist;
  double e;
  double *e_vec;
  int i;
  int j;
  int k;
  int seed;
  double *x;

  seed = seed_init;

  e_vec = new double[n];

  for ( j = 0; j < n; j++ )
  {
    e_vec[j] = 0.0;
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*(dim_num)], 2 );
    }
    e_vec[closest[0]] = e_vec[closest[0]] + dist;

    delete [] x;
  }

  e = 0.0;
  for ( j = 0; j < n; j++ )
  {
    e = e + e_vec[j];
  }
  e = e / double ( ns );

  delete [] e_vec;

  return e;
}
//****************************************************************************80

int file_column_count ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
//
//  Discussion:
//
//    The file is assumed to be a simple text file.
//
//    Most lines of the file is presumed to consist of COLUMN_NUM words, separated
//    by spaces.  There may also be some blank lines, and some comment lines,
//    which have a "#" in column 1.
//
//    The routine tries to find the first non-comment non-blank line and
//    counts the number of words in that line.
//
//    If all lines are blanks or comments, it goes back and tries to analyze
//    a comment line.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the file.
//
//    Output, int FILE_COLUMN_COUNT, the number of columns assumed
//    to be in the file.
//
{
  int column_num;
  ifstream input;
  bool got_one;
  char line[256];
//
//  Open the file.
//
  input.open ( input_filename );

  if ( !input )
  {
    column_num = -1;
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Fatal error!\n";
    cout << "  Could not open the file:\n";
    cout << "  \"" << input_filename << "\"\n";
    return column_num;
  }
//
//  Read one line, but skip blank lines and comment lines.
//
  got_one = false;

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      continue;
    }

    if ( line[0] == '#' )
    {
      continue;
    }

    got_one = true;
    break;

  }

  if ( !got_one )
  {
    input.close ( );

    input.open ( input_filename );

    for ( ; ; )
    {
      input.getline ( line, sizeof ( line ) );

      if ( input.eof ( ) )
      {
        break;
      }

      if ( s_len_trim ( line ) == 0 )
      {
        continue;
      }

      got_one = true;
      break;

    }

  }

  input.close ( );

  if ( !got_one )
  {
    cout << "\n";
    cout << "FILE_COLUMN_COUNT - Warning!\n";
    cout << "  The file does not seem to contain any data.\n";
    return -1;
  }

  column_num = s_word_count ( line );

  return column_num;
}
//****************************************************************************80

int file_row_count ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    FILE_ROW_COUNT counts the number of row records in a file.
//
//  Discussion:
//
//    It does not count lines that are blank, or that begin with a
//    comment symbol '#'.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Output, int FILE_ROW_COUNT, the number of rows found.
//
{
  int bad_num;
  int comment_num;
  ifstream input;
  int i;
  char line[100];
  int record_num;
  int row_num;

  row_num = 0;
  comment_num = 0;
  record_num = 0;
  bad_num = 0;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "FILE_ROW_COUNT - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    exit ( 1 );
  }

  for ( ; ; )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    record_num = record_num + 1;

    if ( line[0] == '#' )
    {
      comment_num = comment_num + 1;
      continue;
    }

    if ( s_len_trim ( line ) == 0 )
    {
      comment_num = comment_num + 1;
      continue;
    }

    row_num = row_num + 1;

  }

  input.close ( );

  return row_num;
}
//****************************************************************************80

void find_closest ( int dim_num, int n, int sample_num, double s[], double r[],
  int nearest[] )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_CLOSEST finds the nearest R point to each S point.
//
//  Discussion:
//
//    This routine finds the closest Voronoi cell generator by checking every
//    one.  For problems with many cells, this process can take the bulk
//    of the CPU time.  Other approaches, which group the cell generators into
//    bins, can run faster by a large factor.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of cell generators.
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, double S[DIM_NUM*SAMPLE_NUM], the points to be checked.
//
//    Input, double R[DIM_NUM*N], the cell generators.
//
//    Output, int NEAREST[SAMPLE_NUM], the (0-based) index of the nearest
//    cell generator.
//
{
  double dist_sq_min;
  double dist_sq;
  int i;
  int jr;
  int js;

  for ( js = 0; js < sample_num; js++ )
  {
    dist_sq_min = r8_huge ( );
    nearest[js] = -1;

    for ( jr = 0; jr < n; jr++ )
    {
      dist_sq = 0.0;
      for ( i = 0; i < dim_num; i++ )
      {
        dist_sq = dist_sq + ( s[i+js*dim_num] - r[i+jr*dim_num] )
                          * ( s[i+js*dim_num] - r[i+jr*dim_num] );
      }

      if ( jr == 0 || dist_sq < dist_sq_min )
      {
        dist_sq_min = dist_sq;
        nearest[js] = jr;
      }
    }
  }

  return;
}
//****************************************************************************80

double gamma_measure ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_MEASURE determines the pointset quality measure GAMMA.
//
//  Discussion:
//
//    The GAMMA measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//      GAMMA = ( GAMMA_MAX / GAMMA_MIN ),
//
//    where
//
//      GAMMA_MAX = maximum ( 1 <= I <= N ) DIST_MIN(I)
//      GAMMA_MIN = minimum ( 1 <= I <= N ) DIST_MIN(I)
//
//    and
//
//      DIST_MIN(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
//
//
//    Note that, in this code, the variable DIST_SQ_MIN is actually the square
//    of the minimum point distance, and so when we compute GAMMA, we must
//    take the square root of the given ratio.
//
//    GAMMA must be at least 1.  For an ideally regular mesh, GAMMA would
//    be equal to one.  Given two meshes, this measure recommends the one
//    with the smaller value of GAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Output, double GAMMA_MEASURE, the GAMMA quality measure.
//
//  Local parameters:
//
//    Local, double GAMMA_SQ_MAX, the maximum, over all points,
//    of the minimum squared distance to a distinct point.
//
//    Local, double GAMMA_SQ_MIN, the minimum, over all points,
//    of the minimum squared distance to a distinct point.
//
{
  int i;
  int j1;
  int j2;
  double dist_sq;
  double dist_sq_min;
  double gamma;
  double gamma_sq_max;
  double gamma_sq_min;
//
//  Take care of ridiculous cases.
//
  if ( n <= 1 )
  {
    gamma = 0.0;
    return gamma;
  }

  gamma_sq_max = 0.0;
  gamma_sq_min = r8_huge ( );

  for ( j1 = 0; j1 < n; j1++ )
  {
    dist_sq_min = r8_huge ( );

    for ( j2 = 0; j2 < n; j2++ )
    {
      if ( j2 != j1 )
      {

        dist_sq = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          dist_sq = dist_sq + pow ( z[i+j1*dim_num] - z[i+j2*dim_num], 2 );
        }

        if ( dist_sq < dist_sq_min )
        {
          dist_sq_min = dist_sq;
        }
      }

    }

    gamma_sq_max = r8_max ( gamma_sq_max, dist_sq_min );
    gamma_sq_min = r8_min ( gamma_sq_min, dist_sq_min );

  }

  if ( gamma_sq_min <= 0.0 )
  {
    gamma = r8_huge ( );
  }
  else
  {
    gamma = sqrt ( gamma_sq_max / gamma_sq_min );
  }

  return gamma;
}
//****************************************************************************80

double h_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    H_MEASURE determines the pointset quality measure H.
//
//  Discussion:
//
//    The H measure of dispersion for a set of N points in an DIM_NUM-dimensional
//    region is the maximum distance between a point in the region and some
//    point in the set.
//
//    To compute this quantity exactly, for every point X in the region,
//    find the nearest element Z of the point set and compute the distance.
//    H is then the maximum of all these distances.
//
//    To ESTIMATE this quantity, carry out the same process, but only for
//    NS sample points in the region.
//
//    Under this measure, a mesh with a smaller value of H is preferable.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double H_MEASURE, the H quality measure.
//
{
  int closest[1];
  double dist;
  double h;
  int i;
  int k;
  int seed;
  double *x;

  seed = seed_init;
  h = 0.0;

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*(dim_num)], 2 );
    }

    h = r8_max ( h, dist );

    delete [] x;
  }

  h = sqrt ( h );

  return h;
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
//  Discussion:
//
//    The sign of 0 and all positive integers is taken to be +1.
//    The sign of all negative integers is -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 May 2003
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

int i4_wrap ( int ival, int ilo, int ihi )

//****************************************************************************80
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
//    13 January 2004
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
//    14 November 2003
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

  if ( 0 < s_len_trim ( title ) )
  {
    cout << "\n";
    cout << title << "\n";
  }

  cout << "\n";
  for ( i = 0; i <= n-1; i++ )
  {
    cout << setw(6) << i + 1 << "  "
         << setw(8) << a[i]  << "\n";
  }

  return;
}
//****************************************************************************80

double lambda_measure ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    LAMBDA_MEASURE determines the pointset quality measure LAMBDA.
//
//  Discussion:
//
//    The LAMBDA measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//    Let
//
//      GAMMA(I) = minimum ( 1 <= J <= N, I /= J ) distance ( Z(I), Z(J) )
//
//    and let
//
//      GAMMA_AVE = sum ( 1 <= I <= N ) GAMMA(I) / N
//
//    then
//
//      LAMBDA = sqrt ( sum ( 1 <= I <= N ) ( GAMMA(I) - GAMMA_AVE )**2 / N )
//        / GAMMA_AVE
//
//    An ideally regular mesh would have GAMMA(I) = GAMMA_AVE for all I,
//    so that LAMBDA would be 0.  Under this measure, the mesh with the
//    smaller value of LAMBDA is to be preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Output, double LAMBDA_MEASURE, the LAMBDA quality measure.
//
//  Local parameters:
//
//    Local, double GAMMA_MAX, the maximum, over all points,
//    of the minimum distance to a distinct point.
//
//    Local, double GAMMA_MIN, the minimum, over all points,
//    of the minimum distance to a distinct point.
//
{
  int i;
  int j;
  double dist;
  double *gamma;
  double gamma_ave;
  double lambda;
//
//  Take care of ridiculous cases.
//
  if ( n <= 1 )
  {
    lambda = 0.0;
    return lambda;
  }
//
//  Compute the minimum spacing between distinct points of the set.
//
  gamma = pointset_spacing ( dim_num, n, z );
//
//  Average the minimum spacing.
//
  gamma_ave = 0.0;
  for ( j = 0; j < n; j++ )
  {
    gamma_ave = gamma_ave + gamma[j];
  }
  gamma_ave = gamma_ave / ( double ) ( n );
//
//  Compute a weighted variance.
//
  if ( gamma_ave <= 0.0 )
  {
    lambda = r8_huge ( );
  }
  else
  {
    lambda = 0.0;
    for ( j = 0; j < n; j++ )
    {
      lambda = lambda + pow ( gamma[j] - gamma_ave, 2 );
    }
    lambda = sqrt ( lambda / ( double ) n );
    lambda = lambda / gamma_ave;
  }

  delete [] gamma;

  return lambda;
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
//    Original FORTRAN77 version by Barry Joe,
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
  double tol = 0.0000001;
  double tolabs;
  int value;
//
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

double mu_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    MU_MEASURE determines the pointset quality measure MU.
//
//  Discussion:
//
//    The MU measure of dispersion for a set of N points in an DIM_NUM-dimensional
//    region takes the ratio of the largest and smallest half-diameters
//    of the Voronoi cells defined by a pointset.
//
//    To compute this quantity exactly, for every point X in the region,
//    find the nearest element Z of the point set and compute the distance.
//
//    Then, for each element Z(I) of the point set, define H(I) to be the
//    maximum of these distances.
//
//    MU is then the ratio of the maximum and minimum values of H.
//
//    To ESTIMATE this quantity, carry out the same process, but only for
//    NS sample points in the region.
//
//    In an ideally regular mesh, MU would be 1.  MU must be at least 1.
//    Under this measure, the mesh with the smaller value of MU is to be
//    preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double MU_MEASURE, the MU quality measure.
//
{
  int k;
  int closest[1];
  double dist;
  double *h;
  double h_max;
  double h_min;
  int i;
  int j;
  double mu;
  int seed;
  double *x;

  h = new double[n];

  seed = seed_init;

  for ( j = 0; j < n; j++ )
  {
    h[j] = 0.0;
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*dim_num], 2 );
    }

    h[closest[0]] = r8_max ( h[closest[0]], dist );

    delete [] x;
  }

  h_max = h[0];
  for ( j = 1; j < n; j++ )
  {
    h_max = r8_max ( h_max, h[j] );
  }
  h_max = sqrt ( h_max );

  h_min = h[0];
  for ( j = 1; j < n; j++ )
  {
    h_min = r8_min ( h_min, h[j] );
  }
  h_min = sqrt ( h_min );

  if ( h_min == 0.0 )
  {
    mu = r8_huge ( );
  }
  else
  {
    mu = h_max / h_min;
  }

  delete [] h;

  return mu;
}
//****************************************************************************80

double nu_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    NU_MEASURE determines the pointset quality measure NU.
//
//  Discussion:
//
//    The NU measure of dispersion for a set of N points in an DIM_NUM-dimensional
//    region is defined as follows:
//
//    For each element Z(I) of the pointset, let VOLUME(I) be the volume
//    of the corresponding Voronoi subregion, restricted to the region.
//
//    Then
//
//      NU = max ( 1 <= I <= N ) VOLUME(I) / min ( 1 <= I <= N ) VOLUME(I)
//
//    This quantity can be estimated by using a large number of sampling
//    points to estimate the Voronoi volumes.
//
//    For an ideally uniform pointset, the Voronoi volumes would be equal,
//    so that NU would be 1.  In any case, NU must be 1 or greater.  In
//    comparing two meshes, the one with the lower value of NU would be
//    preferred.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double NU_MEASURE, the NU quality measure.
//
{
  int closest[1];
  int *hit;
  int j;
  int k;
  double nu;
  int seed;
  double *volume;
  double volume_max;
  double volume_min;
  double *x;

  hit = new int[n];
  volume = new double[n];

  seed = seed_init;

  for ( j = 0; j < n; j++ )
  {
    hit[j] = 0;
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    hit[closest[0]] = hit[closest[0]] + 1;

    delete [] x;
  }

  for ( j = 0; j < n; j++ )
  {
    volume[j] = ( double ) ( hit[j] ) / ( double ) ( ns );
  }

  volume_max = 0.0;
  for ( j = 0; j < n; j++ )
  {
    volume_max = r8_max ( volume_max, volume[j] );
  }

  volume_min = r8_huge ( );
  for ( j = 0; j < n; j++ )
  {
    volume_min = r8_min ( volume_min, volume[j] );
  }

  if ( volume_min == 0.0 )
  {
    nu = r8_huge ( );
  }
  else
  {
    nu = volume_max / volume_min;
  }

  delete [] hit;
  delete [] volume;

  return nu;
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

double *pointset_spacing ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINTSET_SPACING determines the minimum spacing between points in the set.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the point distribution.
//
//    Output, double POINTSET_SPACING(N), the minimum distance between each
//    point and a distinct point in the set.
//
{
  double dist;
  int i;
  int j1;
  int j2;
  double *gamma;

  gamma = new double[n];

  for ( j1 = 0; j1 < n; j1++ )
  {
    gamma[j1] = r8_huge ( );

    for ( j2 = 0; j2 < n; j2++ )
    {
      if ( j2 != j1 )
      {
        dist = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          dist = dist + pow ( z[i+j1*dim_num] - z[i+j2*dim_num], 2 );
        }
        gamma[j1] = r8_min ( gamma[j1], dist );
      }
    }
  }

  for ( j1 = 0; j1 < n; j1++ )
  {
    gamma[j1] = sqrt ( gamma[j1] );
  }

  return gamma;
}
//****************************************************************************80

double q_measure ( int n, double z[], int triangle_order, int triangle_num,
  int triangle_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    Q_MEASURE determines the triangulated pointset quality measure Q.
//
//  Discussion:
//
//    The Q measure evaluates the uniformity of the shapes of the triangles
//    defined by a triangulated pointset.
//
//    For a single triangle T, the value of Q(T) is defined as follows:
//
//      TAU_IN = radius of the inscribed circle,
//      TAU_OUT = radius of the circumscribed circle,
//
//      Q(T) = 2 * TAU_IN / TAU_OUT
//        = ( B + C - A ) * ( C + A - B ) * ( A + B - C ) / ( A * B * C )
//
//    where A, B and C are the lengths of the sides of the triangle T.
//
//    The Q measure computes the value of Q(T) for every triangle T in the
//    triangulation, and then computes the minimum of this
//    set of values:
//
//      Q_MEASURE = min ( all T in triangulation ) Q(T)
//
//    In an ideally regular mesh, all triangles would have the same
//    equilateral shape, for which Q = 1.  A good mesh would have
//    0.5 < Q.
//
//    Given the 2D coordinates of a set of N nodes, stored as Z(1:2,1:N),
//    a triangulation is a list of TRIANGLE_NUM triples of node indices that form
//    triangles.  Generally, a maximal triangulation is expected, namely,
//    a triangulation whose image is a planar graph, but for which the
//    addition of any new triangle would mean the graph was no longer planar.
//    A Delaunay triangulation is a maximal triangulation which maximizes
//    the minimum angle that occurs in any triangle.
//
//    The code has been modified to 'allow' 6-node triangulations.
//    However, no effort is made to actually process the midside nodes.
//    Only information from the vertices is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//    Per-Olof Persson and Gilbert Strang,
//    A Simple Mesh Generator in MATLAB,
//    SIAM Review,
//    Volume 46, Number 2, pages 329-345, June 2004.
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, double Z[2*N], the points.
//
//    Input, int TRIANGLE_ORDER, the order of the triangles.
//
//    Input, int TRIANGLE_NUM, the number of triangles.
//
//    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
//    the triangulation.
//
//    Output, double Q_MEASURE, the Q quality measure.
//
{
  int a_index;
  double ab_length;
  int b_index;
  double bc_length;
  int c_index;
  double ca_length;
  double q;
  double q_min;
  int triangle;
  double value;

  if ( triangle_num < 1 )
  {
    value = -1.0;
    return value;
  }

  q_min = r8_huge ( );

  for ( triangle = 0; triangle < triangle_num; triangle++ )
  {
    a_index = triangle_node[0+triangle*3];
    b_index = triangle_node[1+triangle*3];
    c_index = triangle_node[2+triangle*3];

    ab_length = sqrt (
        pow ( z[0+(a_index-1)*2] - z[0+(b_index-1)*2], 2 )
      + pow ( z[1+(a_index-1)*2] - z[1+(b_index-1)*2], 2 ) );

    bc_length = sqrt (
        pow ( z[0+(b_index-1)*2] - z[0+(c_index-1)*2], 2 )
      + pow ( z[1+(b_index-1)*2] - z[1+(c_index-1)*2], 2 ) );

    ca_length = sqrt (
        pow ( z[0+(c_index-1)*2] - z[0+(a_index-1)*2], 2 )
      + pow ( z[1+(c_index-1)*2] - z[1+(a_index-1)*2], 2 ) );

    q = ( bc_length + ca_length - ab_length )
      * ( ca_length + ab_length - bc_length )
      * ( ab_length + bc_length - ca_length )
      / ( ab_length * bc_length * ca_length );

    q_min = r8_min ( q_min, q );
  }

  value = q_min;

  return value;
}
//****************************************************************************80

double r0_measure ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    R0_MEASURE determines the pointset quality measure R0.
//
//  Discussion:
//
//    The R0 measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//      R0 = sum ( 1 <= I /= J <= N ) log ( 1 / distance ( Z(I), Z(J) ) )
//         / ( N * ( N - 1 ) )
//
//    The divisor of ( N * ( N - 1 ) ) means that R0 is essentially an
//
//    R0 is undefined if N < 2 or if any two points are equal.
//
//    R0 is known as the Riesz S-energy for S = 0.
//
//    Given two meshes, this measure recommends the one with the smaller
//    value of R0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    D P Hardin and E B Saff,
//    Discretizing Manifolds via Minimum Energy Points,
//    Notices of the AMS,
//    Volume 51, Number 10, November 2004, pages 1186-1194.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Output, double R0_MEASURE, the R0 quality measure.
//
{
  double dist;
  int i;
  int j1;
  int j2;
  double value;
//
//  Take care of ridiculous cases.
//
  if ( n <= 1 )
  {
    value = r8_huge ( );
    return value;
  }

  value = 0.0;

  for ( j1 = 0; j1 < n; j1++ )
  {
    for ( j2 = 0; j2 < n; j2++ )
    {
      if ( j2 != j1 )
      {
        dist = 0.0;
        for ( i = 0; i < dim_num; i++ )
        {
          dist = dist + pow ( z[i+j1*dim_num] - z[i+j2*dim_num], 2 );
        }
        dist = sqrt ( dist );

        if ( dist == 0.0 )
        {
          value = r8_huge ( );
          return value;
        }

        value = value + log ( 1.0 / dist );
      }
    }
  }

  value = value / ( double ) ( n * ( n - 1 ) );

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

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

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
//    Output, double R8_MAX, the maximum of X and Y.
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
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  }
  else
  {
    return x;
  }
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a portable pseudorandom number generator.
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
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
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

  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r82vec_permute ( int n, double a[], int p[] )

//****************************************************************************80
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
          cout << "  IGET = " << iget << "\n";
          cout << "  N =    " << n << "\n";
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

int *r82vec_sort_heap_index_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R82VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R28VEC.
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

bool r8mat_in_01 ( int m, int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_IN_01 is TRUE if the entries of an R8MAT are in the range [0,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2004
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
//    Input, double A[M*N], the matrix.
//
//    Output, bool R8MAT_IN_01, is TRUE if every entry of A is
//    between 0 and 1.
//
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      if ( a[i+j*m] < 0.0 || 1.0 < a[i+j*m] )
      {
        return false;
      }
    }
  }

  return true;
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

  return;
# undef INCX
}
//****************************************************************************80

double *r8mat_uniform_01 ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_UNIFORM_01 fills an R8MAT with pseudorandom values.
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
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double R8MAT_UNIFORM_01[M*N], a matrix of pseudorandom values.
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

      r[i+j*m] = ( double ) ( *seed ) * 4.656612875E-10;
    }
  }

  return r;
}
//****************************************************************************80

double r8vec_max ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MAX returns the value of the maximum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], a pointer to the first entry of the array.
//
//    Output, double R8VEC_MAX, the value of the maximum element.
//
{
  int i;
  double *r8vec_pointer;
  double value;

  value = - r8_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( value < r8vec[i] )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double r8vec_min ( int n, double r8vec[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_MIN returns the value of the minimum element in an R8VEC.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double R8VEC[N], the array to be checked.
//
//    Output, double R8VEC_MIN, the value of the minimum element.
//
{
  int i;
  double *r8vec_pointer;
  double value;

  value = r8_huge ( );

  if ( n <= 0 )
  {
    return value;
  }

  for ( i = 0; i < n; i++ )
  {
    if ( r8vec[i] < value )
    {
      value = r8vec[i];
    }
  }
  return value;
}
//****************************************************************************80

double *r8vec_normal_01 ( int n, int *seed )

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
//    18 October 2004
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
//    Output, double X(N), a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, real R(N+1), is used to store some uniform random values.
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
//    Local, real Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
# define PI 3.141592653589793

  int i;
  int m;
  static int made = 0;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;

  x = new double[n];
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
    r = r8vec_uniform_01 ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * PI * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * PI * r[1] );

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

    r = r8vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
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

    r = r8vec_uniform_01 ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * PI * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * PI * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
# undef PI
}
//****************************************************************************80

double *r8vec_uniform_01 ( int n, int *seed )

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

    r[i] = ( double ) ( *seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

double *radius_maximus ( int dim_num, int n, double z[], bool walls )

//****************************************************************************80
//
//  Purpose:
//
//    RADIUS_MAXIMUS finds the biggest possible nonintersecting sphere.
//
//  Discussion:
//
//    We are given a set of N points in DIM_NUM space.  We imagine that
//    at each point simultaneously, a sphere begins to expand.
//    Each sphere stops expanding as soon as it touches another sphere.
//    The radius of these spheres is to be computed.
//
//    If WALLS is true, then the spheres must not extend outside the
//    "walls" of the unit hypersquare.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the point coordinates.
//    If WALLS is TRUE, these values must be between 0 and 1.
//
//    Input, logical WALLS, is TRUE if the spheres must not extend
//    outside the unit hypercube.  If WALLS is FALSE, then this
//    restriction is not imposed.
//
//    Output, double RADIUS(N), the radius of the
//    maximal nonintersecting sphere around each point.
//
{
  double distance_j;
  double distance_min;
  bool done;
  int FIXED = 0;
  int FREE = 1;
  int i;
  int j;
  int j1;
  int j2;
  int next;
  double *radius;
  double radius_i;
  double radius_min;
  int *status;

  radius = new double[n];
  status = new int[n];

  if ( walls )
  {
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        if ( z[i+j*dim_num] < 0.0 )
        {
          cout << "\n";
          cout << "RADIUS_MAXIMUS - Fatal error!\n";
          cout << "  Some coordinate is less than 0.\n";
          exit ( 1 );
        }
        else if ( 1.0 < z[i+j*dim_num] )
        {
          cout << "\n";
          cout << "RADIUS_MAXIMUS - Fatal error!\n";
          cout << "  Some coordinate is greater than 1.\n";
          exit ( 1 );
        }
      }
    }

  }
//
//  Initially, all points are "free".
//
  for ( j = 0; j < n; j++ )
  {
    radius[j] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    status[j] = FREE;
  }

  for ( ; ; )
  {
//
//  If all points are fixed, we're done.
//
    done = true;

    for ( j = 0; j < n; j++ )
    {
      if ( status[j] != FIXED )
      {
        done = false;
        break;
      }
    }

    if ( done )
    {
      break;
    }
//
//  Look at all the free points.
//  Imagine an expanding sphere at each free point, and determine
//  which such sphere will first have to stop expanding.
//
    next = -1;
    radius_min = r8_huge ( );

    for ( j1 = 0; j1 < n; j1++ )
    {
      if ( status[j1] == FREE )
      {
        if ( walls )
        {
          radius_i = r8_huge ( );
          for ( i = 0; i < dim_num; i++ )
          {
            radius_i = r8_min ( radius_i, z[i+j1*dim_num] );
          }
          for ( i = 0; i < dim_num; i++ )
          {
            radius_i = r8_min ( radius_i, 1.0 - z[i+j1*dim_num] );
          }
        }
        else
        {
          radius_i = r8_huge ( );
        }

        for ( j2 = 0; j2 < n; j2++ )
        {
          if ( j2 != j1 )
          {
            distance_j = 0.0;
            for ( i = 0; i < dim_num; i++ )
            {
              distance_j = distance_j
                + pow ( z[i+j1*dim_num] - z[i+j2*dim_num], 2 );
            }
            distance_j = sqrt ( distance_j );

            if ( status[j2] == FREE )
            {
              radius_i = r8_min ( radius_i, distance_j / 2.0 );
            }
            else
            {
              radius_i = r8_min ( radius_i, distance_j - radius[j2] );
            }
          }
        }

        if ( radius_i < radius_min )
        {
          next = j1;
          radius_min = radius_i;
        }

      }

    }

    if ( next == -1 )
    {
      cout << "\n";
      cout << "RADIUS_MAXIMUS - Fatal error!\n";
      cout << "  There were points left to handle, but could\n";
      cout << "  not choose the next one to work on.\n";
      exit ( 1 );
    }

    radius[next] = radius_min;
    status[next] = FIXED;
  }

  return radius;
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
  char *t;

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

double s_to_r8 ( char *s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( char *s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;
  double x;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

int s_word_count ( char *s )

//****************************************************************************80
//
//  Purpose:
//
//    S_WORD_COUNT counts the number of "words" in a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be examined.
//
//    Output, int S_WORD_COUNT, the number of "words" in the string.
//    Words are presumed to be separated by one or more blanks.
//
{
  bool blank;
  int i;
  int nword;

  nword = 0;
  blank = true;

  while ( *s )
  {
    if ( *s == ' ' )
    {
      blank = true;
    }
    else if ( blank )
    {
      nword = nword + 1;
      blank = false;
    }
    *s++;
  }

  return nword;
}
//****************************************************************************80

double *sample_hypercube_uniform ( int dim_num, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_HYPERCUBE_UNIFORM returns sample points in the unit hypercube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points to compute.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SAMPLE_HYPERCUBE_UNIFORM[DIM_NUM*N], the sample points.
//
{
  double *x;

  x = r8mat_uniform_01 ( dim_num, n, seed );

  return x;
}
//****************************************************************************80

double *sample_sphere_uniform ( int m, int n, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    SAMPLE_SPHERE_UNIFORM samples points inside the unit sphere.
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
//    31 October 2004
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
//    Wiley, 1986, page 232.
//
//  Parameters:
//
//    Input, int M, the dimension of the space.
//
//    Input, int N, the number of points.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double SAMPLE_SPHERE_UNIFORM[M*N], the points.
//
{
  double exponent;
  int i;
  int j;
  double norm;
  double r;
  double *x;
  double *y;

  x = new double[m*n];
  y = new double[m];

  exponent = 1.0 / ( double ) ( m );

  for ( j = 0; j < n; j++ )
  {
//
//  Fill a vector with normally distributed values.
//
    y = r8vec_normal_01 ( m, seed );
//
//  Compute the length of the vector.
//
    norm = 0.0;
    for ( i = 0; i < m; i++ )
    {
      norm = norm + y[i] * y[i];
    }
    norm = sqrt ( norm );
//
//  Normalize the vector.
//
    for ( i = 0; i < m; i++ )
    {
      y[i] = y[i] / norm;
    }
//
//  Now compute a value to map the point ON the sphere INTO the sphere.
//
    r = r8_uniform_01 ( seed );
    r = pow ( r, exponent );

    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = r * y[i];
    }
  }

  return x;
}
//****************************************************************************80

double sphere_measure ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_MEASURE determines the pointset quality measure S.
//
//  Discussion:
//
//    This routine computes a measure of even spacing for a set of N points
//    in the DIM_NUM-dimensional unit hypercube.  We will discuss the program
//    as though the space is 2-dimensional and the spheres are circles, but
//    the program may be used for general DIM_NUM-dimensional data.
//
//    The points are assumed to lie in the unit square.
//
//    The program makes a circle-packing measurement on the points
//    by assuming that, at each point, a circle is centered; all
//    the circles start out with zero radius, and then expand
//    together at the same rate.  A circle stops expanding as soon
//    as it touches any other circle.
//
//    The amount of area covered by the circles is compared to the
//    area of the unit square.  This measurement has a certain amount
//    of boundary effect: some circles will naturally extend outside
//    the unit hypercube.  If this is a concern, is possible to restrict
//    the circles to remain inside the unit hypercube.  In any case,
//    this problem generally goes away as the number of points increases.
//
//    Since we are interested in the coverage of the unit hypercube,
//    it is probably best if the circles are restricted.  This way,
//    computing the area of the circles gives a measure of the even
//    coverage of the region, relative to the presumably best possible
//    covering, by the same number of circles, but of equal radius.
//
//    In the limit, the maximum relative packing density of a 2D
//    region with equal-sized circles is 0.9069.  In 3D, a density
//    of at least 0.74 can be achieved, and it is known that no
//    greater than 0.7796 is possible.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of points.
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double Z[DIM_NUM*N], the points.
//
//    Output, double SPHERE_MEASURE, the amount of volume taken up
//    by the nonintersecting spheres of maximum radius around each
//    point.  Ignoring boundary effects, the "ideal" value would be
//    1 (achievable only in 1 dimension), and the maximum value
//    possible is the sphere packing density in the given spatial
//    dimension.  If boundary effects can be ignored, the value of
//    SPHERE_VOLUME reports how closely the given set of points
//    behaves like a set of close-packed spheres.
//
//  Local Parameters:
//
//    Local, logical WALLS, is TRUE if the spheres are restricted
//    to lie within the unit hypercube.
//
{
  int i;
  int j;
  double *radius;
  double radius_ave;
  double radius_max;
  double radius_min;
  double sphere;
  bool verbose = false;
  double volume;
  bool walls = true;

  if ( !r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "\n";
    cout << "SPHERE_MEASURE - Fatal error!\n";
    cout << "  Some of the data is not inside the unit hypercube.\n";
    return r8_huge ( );
  }

  radius = radius_maximus ( dim_num, n, z, walls );

  sphere = 0.0;
  for ( i = 0; i < n; i++ )
  {
    volume = sphere_volume_nd ( dim_num, radius[i] );
    sphere = sphere + volume;
  }

  if ( verbose )
  {
    radius_ave = 0.0;
    radius_min = r8_huge ( );
    radius_max = 0.0;
    for ( j = 0; j < n; j++ )
    {
      radius_ave = radius_ave + radius[j];
      radius_min = r8_min ( radius_min, radius[j] );
      radius_max = r8_max ( radius_max, radius[j] );
    }

    cout << "\n";
    cout << "  Number of dimensions is " << dim_num << "\n";
    cout << "  Number of points is " << n << "\n";
    if ( walls )
    {
      cout << "  Spheres are required to stay in the unit hypercube.\n";
    }
    else
    {
      cout << "  Spheres are NOT required to stay in the unit hypercube.\n";
    }
    cout << "\n";
    cout << "  Average radius = " << radius_ave << "\n";
    cout << "  Minimum radius = " << radius_min << "\n";
    cout << "  Maximum radius = " << radius_max << "\n";
    cout << "  Sphere volume =  " << sphere     << "\n";
  }

  delete [] radius;

  return sphere;
}
//****************************************************************************80

double sphere_volume_nd ( int dim_num, double r )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_VOLUME_ND computes the volume of a sphere in ND.
//
//  Discussion:
//
//    A sphere in ND satisfies the equation:
//
//      sum ( ( X(1:N) - XC(1:N) )**2 ) = R**2
//
//    where R is the radius and XC is the center.
//
//    N  Volume
//
//    2             PI    * R**2
//    3  (4/3)    * PI    * R**3
//    4  (1/2)    * PI**2 * R**4
//    5  (8/15)   * PI**2 * R**5
//    6  (1/6)    * PI**3 * R**6
//    7  (16/105) * PI**3 * R**7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the dimension of the space.
//
//    Input, double R, the radius of the sphere.
//
//    Output, double SPHERE_VOLUME, the volume of the sphere.
//
{
# define PI 3.141592653589793

  int i;
  int m;
  double volume;

  if ( ( dim_num % 2 ) == 0 )
  {
    m = dim_num / 2;
    volume = pow ( PI, m );
    for ( i = 1; i <= m; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }
  else
  {
    m = ( dim_num - 1 ) / 2;
    volume = pow ( PI, m ) * pow ( 2.0, dim_num );
    for ( i = m+1; i <= 2*m+1; i++ )
    {
      volume = volume / ( double ) ( i );
    }
  }

  volume = volume * pow ( r, dim_num );

  return volume;
# undef PI
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe,
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

double tau_measure ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ),
  int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TAU_MEASURE determines the pointset quality measure TAU.
//
//  Discussion:
//
//    The TAU measure of point distribution quality for a set Z of
//    N points in an DIM_NUM-dimensional region is defined as follows:
//
//    For each point Z(I) in the pointset, let V(I) be the subregion
//    defined by the intersection of the region with the Voronoi
//    region associated with Z(I).
//
//    Let T(I) be the trace of the second moment tensor about the point
//    Z(I), associated with the subregion V(I).  Let T_BAR be the average
//    of the values of T(1:N).
//
//    Then TAU = maximum ( 1 <= I <= N ) abs ( T(I) - TBAR ).
//
//    This quantity can be estimated using sampling.  A given number of
//    sample points are generated in the region, assigned to the nearest
//    element of the pointset, and used to approximate the Voronoi regions
//    and the second moment tensors.
//
//    In an ideally regular mesh, the values of T would be equal, and so
//    TAU would be zero.  In general, the smaller TAU, the better.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 November 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Gunzburger and John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double Z[DIM_NUM*N], the point distribution.
//
//    Input, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region,
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double TAU_MEASURE, a quality measure.
//
{
  double *centroid;
  int closest[1];
  int *hit;
  int i;
  int i1;
  int i2;
  int j;
  int k;
  double *moment;
  int seed;
  double *t;
  double t_bar;
  double tau;
  double *x;

  centroid = new double[dim_num*n];
  hit = new int[n];
  moment = new double[dim_num*dim_num*n];
  t = new double[n];

  seed = seed_init;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      centroid[i+j*dim_num] = 0.0;
    }
  }

  for ( j = 0; j < n; j++ )
  {
    hit[j] = 0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i2 = 0; i2 < dim_num; i2++ )
    {
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        moment[i1+i2*dim_num+j*dim_num*dim_num] = 0.0;
      }
    }
  }

  for ( k = 1; k <= ns; k++ )
  {
    x = sample_routine ( dim_num, 1, &seed );

    find_closest ( dim_num, n, 1, x, z, closest );

    hit[closest[0]] = hit[closest[0]] + 1;

    for ( i = 0; i < dim_num; i++ )
    {
      centroid[i+closest[0]*dim_num] = centroid[i+closest[0]*dim_num] + x[i];
    }

    for ( i1 = 0; i1 < dim_num; i1++ )
    {
      for ( i2 = 0; i2 < dim_num; i2++ )
      {
        moment[i1+i2*dim_num+closest[0]*dim_num*dim_num]
        = moment[i1+i2*dim_num+closest[0]*dim_num*dim_num] + x[i1] * x[i2];
      }
    }
    delete [] x;
  }

  for ( j = 0; j < n; j++ )
  {
    if ( 0 < hit[j] )
    {
      for ( i = 0; i < dim_num; i++ )
      {
        centroid[i+j*dim_num] = centroid[i+j*dim_num] / ( double ) ( hit[j] );
      }
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        for ( i2 = 0; i2 < dim_num; i2++ )
        {
          moment[i1+i2*dim_num+j*dim_num*dim_num]
            = moment[i1+i2*dim_num+j*dim_num*dim_num] / ( double ) ( hit[j] );
        }
      }
      for ( i1 = 0; i1 < dim_num; i1++ )
      {
        for ( i2 = 0; i2 < dim_num; i2++ )
        {
            moment[i1+i2*dim_num+j*dim_num*dim_num] =
            moment[i1+i2*dim_num+j*dim_num*dim_num]
              - centroid[i1+j*dim_num] * centroid[i2+j*dim_num];
        }
      }
    }
  }

  for ( j = 0; j < n; j++ )
  {
    t[j] = 0.0;
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      t[j] = t[j] + moment[i+i*dim_num+j*dim_num*dim_num];
    }
  }

  t_bar = 0.0;

  for ( j = 0; j < n; j++ )
  {
    t_bar = t_bar + t[j];
  }
  t_bar = t_bar / ( double ) ( n );

  tau = 0.0;
  for ( j = 0; j < n; j++ )
  {
    tau = r8_max ( tau, fabs ( t[j] - t_bar ) );
  }

  delete [] centroid;
  delete [] hit;
  delete [] moment;
  delete [] t;

  return tau;
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
//    02 October 2003
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

char *timestring ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTRING returns the current YMDHMS date as a string.
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
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, char *TIMESTRING, a string containing the current YMDHMS date.
//
{
# define TIME_SIZE 40

  const struct tm *tm;
  size_t len;
  time_t now;
  char *s;

  now = time ( NULL );
  tm = localtime ( &now );

  s = new char[TIME_SIZE];

  len = strftime ( s, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  return s;
# undef TIME_SIZE
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
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2003
//
//  Author:
//
//    Original FORTRAN77 version by Barry Joe,
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
