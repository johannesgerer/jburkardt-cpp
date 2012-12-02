# ifdef ANSI_HEADERS
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

//****************************************************************************80

double chi_measure ( int dim_num, int n, double z[], int ns, int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_MEASURE determines the pointset quality measure CHI.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Discussion:
//
//    The CHI measure of point distribution quality for a set Z of
//    N points in the DIM_NUM-dimensional unit hypercube is defined as follows:
//
//    Assign every point X in the unit hypercube to the nearest element 
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
//    number of points in the unit hypercube, rather than all of them.
//
//    For an ideally regular mesh, all the CHI(I)'s will be equal.
//    Any deviation from regularity increases the value of some entries
//    of CHI; thus, given two meshes, the one with a lower value of
//    CHI is to be recommended.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    Max Gunzburger
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
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double CHI_MEASURE, the CHI quality measure.
//
{
  double chi;
  double *chi_vec;
  int *closest;
  double dist;
  double *gamma;
  double *h;
  int i;
  int j;
  int k;
  int seed;
  double *x;

  if ( !r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "\n";
    cout << "CHI_MEASURE - Fatal error!\n";
    cout << "  Some of the data is not inside the unit hypercube.\n";
    return r8_huge ( );
  }

  seed = seed_init;

  chi_vec = new double[n];
  closest = new int[1];
  h = new double[n];
  x = new double[dim_num];

  for ( j = 0; j < n; j++ )
  {
    h[j] = 0.0;
  }

  for ( k = 1; k <= ns; k++ )
  {
    r8vec_uniform_01 ( dim_num, &seed, x );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*(dim_num)], 2 );
    }
    h[closest[0]] = r8_max ( h[closest[0]], dist );
  }

  gamma = pointset_spacing ( dim_num, n, z );

  chi = 0.0;

  for ( j = 0; j < n; j++ )
  {
    chi_vec[j] = 2.0 * sqrt ( h[j] ) / gamma[j];
    chi = r8_max ( chi, chi_vec[j] );
  }

  delete [] chi_vec;
  delete [] closest;
  delete [] gamma;
  delete [] h;
  delete [] x;

  return chi;
}
//****************************************************************************80

double d_measure ( int dim_num, int n, double z[], int ns, int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    D_MEASURE determines the pointset quality measure D.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Discussion:
//
//    The D measure of point distribution quality for a set Z of
//    N points in the DIM_NUM-dimensional unit hypercube is defined as follows:
//
//    For each point Z(I) in the pointset, let V(I) be the region
//    defined by the intersection of the unit hypercube with the Voronoi
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
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    Max Gunzburger
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
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double D_MEASURE, the D quality measure.
//
{
  double *a;
  double *centroid;
  int *closest;
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

  if ( !r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "\n";
    cout << "D_MEASURE - Fatal error!\n";
    cout << "  Some of the data is not inside the unit hypercube.\n";
    return r8_huge ( );
  }

  a = new double[dim_num*dim_num];
  centroid = new double[dim_num*n];
  closest = new int[1];
  hit = new int[n];
  moment = new double[dim_num*dim_num*n];
  tri = new double[n];
  x = new double[dim_num];

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
    r8vec_uniform_01 ( dim_num, &seed, x );

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
          moment[i1+i2*dim_num+j*dim_num*dim_num] = moment[i1+i2*dim_num+j*dim_num*dim_num] 
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
      moment[i+i*dim_num+j*dim_num*dim_num] = moment[i+i*dim_num+j*dim_num*dim_num] 
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
  delete [] closest;
  delete [] hit;
  delete [] moment;
  delete [] tri;
  delete [] x;

  return d;
}
//****************************************************************************80

double h_measure ( int dim_num, int n, double z[], int ns, int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    H_MEASURE determines the pointset quality measure H.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Discussion:
//
//    The H measure of dispersion for a set of N points in the DIM_NUM-dimensional
//    unit hypercube is the maximum distance between a point in the unit 
//    hypercube and some point in the set.
//
//    To compute this quantity exactly, for every point X in the unit hypercube,
//    find the nearest element Z of the point set and compute the distance.  
//    H is then the maximum of all these distances.
//
//    To ESTIMATE this quantity, carry out the same process, but only for
//    NS sample points in the region. 
//
//    Under this measure, a mesh with a smaller value of H is preferable.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    Max Gunzburger
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
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double H_MEASURE, the H quality measure.
//
{
  int *closest;
  double dist;
  double h;
  int i;
  int k;
  int seed;
  double *x;

  if ( !r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "\n";
    cout << "H_MEASURE - Fatal error!\n";
    cout << "  Some of the data is not inside the unit hypercube.\n";
    return r8_huge ( );
  }

  seed = seed_init;

  closest = new int[1];
  x = new double[dim_num];

  h = 0.0;

  for ( k = 1; k <= ns; k++ )
  {
    r8vec_uniform_01 ( dim_num, &seed, x );

    find_closest ( dim_num, n, 1, x, z, closest );

    dist = 0.0;
    for ( i = 0; i < dim_num; i++ )
    {
      dist = dist + pow ( x[i] - z[i+closest[0]*(dim_num)], 2 );
    }

    h = r8_max ( h, dist );
  }

  h = sqrt ( h );

  delete [] closest;
  delete [] x;

  return h;
}
//****************************************************************************80

double *pointset_spacing ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    POINTSET_SPACING determines the minimum spacing between points in the set.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    Max Gunzburger
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

double tau_measure ( int dim_num, int n, double z[], int ns, int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TAU_MEASURE determines the pointset quality measure TAU.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Discussion:
//
//    The TAU measure of point distribution quality for a set Z of
//    N points in the DIM_NUM-dimensional unit hypercube is defined as follows:
//
//    For each point Z(I) in the pointset, let V(I) be the region
//    defined by the intersection of the unit hypercube with the Voronoi
//    region associated with Z(I).
//
//    Let T(I) be the trace of the second moment tensor about the point
//    Z(I), associated with the region V(I).  Let T_BAR be the average 
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
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    Max Gunzburger
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
//    Input, int SEED_INIT, the initial value of the random number seed.
//
//    Output, double TAU_MEASURE, a quality measure.
//
{
  double *centroid;
  int *closest;
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

  if ( !r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "\n";
    cout << "TAU_MEASURE - Fatal error!\n";
    cout << "  Some of the data is not inside the unit hypercube.\n";
    return r8_huge ( );
  }

  centroid = new double[dim_num*n];
  closest = new int[1];
  hit = new int[n];
  moment = new double[dim_num*dim_num*n];
  t = new double[n];
  x = new double[dim_num];

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
    r8vec_uniform_01 ( dim_num, &seed, x );

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
            moment[i1+i2*dim_num+j*dim_num*dim_num] = moment[i1+i2*dim_num+j*dim_num*dim_num] 
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
  delete [] closest;
  delete [] hit;
  delete [] moment;
  delete [] t;
  delete [] x;

  return tau;
}
