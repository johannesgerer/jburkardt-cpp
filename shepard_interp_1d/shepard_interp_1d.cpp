# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "shepard_interp_1d.hpp"
# include "r8lib.hpp"

//****************************************************************************80

double *shepard_basis_1d ( int nd, double xd[], int k, double p, int ni, 
  double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the data points.
//
//    Input, int K, the index of the desired basis function,
//    0 <= K <= ND - 1.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double SHEPARD_BASIS_1D[NI], the basis function at the interpolation 
//    points.
// 
{
  double *bk;
  int i;
  int j;
  double s;
  double *w;
  int z;

  w = new double[nd];
  bk = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    if ( p == 0.0 )
    {
      for ( j = 0; j < nd; j++ )
      {
        w[j] = 1.0 / ( double ) ( nd );
      }
    }
    else
    {
      z = -1;
      for ( j = 0; j < nd; j++ )
      {
        w[j] = r8_abs ( xi[i] - xd[j] );
        if ( w[j] == 0.0 )
        {
          z = j;
          break;
        }
      }

      if ( z != -1 )
      {
        for ( j = 0; j < nd; j++ )
        {
          w[j] = 0.0;
        }
        w[z] = 1.0;
      }
      else
      {
        for ( j = 0; j < nd; j++ )
        {
          w[j] = 1.0 / pow ( w[j], p );
        }
        s = r8vec_sum ( nd, w );
        for ( j = 0; j < nd; j++ )
        {
          w[j] = w[j] / s;
        }
      }
    }
    bk[i] = w[k];
  }

  delete [] w;

  return bk;
}
//****************************************************************************80

double *shepard_interp_1d ( int nd, double xd[], double yd[], double p, int ni, 
  double xi[] )

//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_INTERP_1D evaluates a 1D Shepard interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the data points.
//
//    Input, double YD[ND], the data values.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double SHEPARD_INTERP_1D[NI], the interpolated values.
//
{ 
  int i;
  int j;
  int k;
  double s;
  double *w;
  double *yi;
  int z;

  w = new double[nd];
  yi = new double[ni];

  for ( i = 0; i < ni; i++ )
  {
    if ( p == 0.0 )
    {
      for ( j = 0; j < nd; j++ )
      {
        w[j] = 1.0 / ( double ) ( nd );
      }
    }
    else
    {
      z = -1;
      for ( j = 0; j < nd; j++ )
      {
        w[j] = r8_abs ( xi[i] - xd[j] );
        if ( w[j] == 0.0 )
        {
          z = j;
          break;
        }
      }

      if ( z != -1 )
      {
        for ( j = 0; j < nd; j++ )
        {
          w[j] = 0.0;
        }
        w[z] = 1.0;
      }
      else
      {
        for ( j = 0; j < nd; j++ )
        {
          w[j] = 1.0 / pow ( w[j], p );
        }
        s = r8vec_sum ( nd, w );
        for ( j = 0; j < nd; j++ )
        {
          w[j] = w[j] / s;
        }
      }
    }
    yi[i] = r8vec_dot_product ( nd, w, yd );
  }
  delete [] w;

  return yi;
}
