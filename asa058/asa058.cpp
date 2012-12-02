# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>

using namespace std;

# include "asa058.hpp"

//****************************************************************************80

void clustr ( double x[], double d[], double dev[], int b[], double f[], 
  int e[], int i, int j, int n, int nz, int k )

//****************************************************************************80
//
//  Purpose:
//
//    CLUSTR uses the K-means algorithm to cluster data.
//
//  Discussion:
//
//    Given a matrix of I observations on J variables, the
//    observations are allocated to N clusters in such a way that the
//    within-cluster sum of squares is minimised.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by David Sparks.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    David Sparks,
//    Algorithm AS 58:
//    Euclidean Cluster Analysis,
//    Applied Statistics,
//    Volume 22, Number 1, 1973, pages 126-130.
//
//  Parameters:
//
//    Input, double X[I*J], the observed data.
//
//    Input/output, double D[K*J], the cluster centers.
//    On input, the user has chosen these.  On output, they have been
//    updated.
//
//    Output, double DEV[K], the sums of squared deviations
//    of observations from their cluster centers.
//
//    Output, int B[I], indicates the cluster to which
//    each observation has been assigned.
//
//    Workspace, double F[I].
//
//    Output, int E[K], the number of observations assigned
//    to each cluster.
//
//    Input, int I, the number of observations.
//
//    Input, int J, the number of variables.
//
//    Input, int N, the number of clusters.
//
//    Input, int NZ, the minimum number of observations
//    which any cluster is allowed to have.
//
//    Input, int K, the maximum number of clusters.
//
{
  double big = 1.0E+10;
  double da;
  double db;
  double dc;
  double de;
  double fl;
  double fm;
  double fq;
  int ia;
  int ic;
  int id;
  int ie;
  int ig;
  int ih;
  int ii;
  int ij;
  int ik;
  int il;
  int in;
  int ip;
  int ir;
  int is;
  int it;
  int iu;
  int iw;
  int ix;
  int iy;

  for ( ia = 1; ia <= n; ia++ )
  {
    e[ia-1] = 0;
  }
//
//  For each observation, calculate the distance from each cluster
//  center, and assign to the nearest.
//
  for ( ic = 1; ic <= i; ic++ )
  {
    f[ic-1] = 0.0;
    da = big;

    for ( id = 1; id <= n; id++ )
    {
      db = 0.0;
      for ( ie = 1; ie <= j; ie++ )
      {
        dc = x[ic-1+(ie-1)*i] - d[id-1+(ie-1)*k];
        db = db + dc * dc;
      }

      if ( db < da )
      {
        da = db;
        b[ic-1] = id;
      }
    }
    ig = b[ic-1];
    e[ig-1] = e[ig-1] + 1;
  }
//
//  Calculate the mean and sum of squares for each cluster.
//
  for ( ix = 1; ix <= n; ix++ )
  {
    dev[ix-1] = 0.0;
    for ( iy = 1; iy <= j; iy++ )
    {
      d[ix-1+(iy-1)*k] = 0.0;
    }
  }

  for ( ic = 1; ic <= i; ic++ )
  {
    ig = b[ic-1];
    for ( ih = 1; ih <= j; ih++ )
    {
      d[ig-1+(ih-1)*k] = d[ig-1+(ih-1)*k] + x[ic-1+(ih-1)*i];
    }
  }

  for ( ij = 1; ij <= j; ij++ )
  {
    for ( ii = 1; ii <= n; ii++ )
    {
      d[ii-1+(ij-1)*k] = d[ii-1+(ij-1)*k] / ( double ) e[ii-1];
    }
  }

  for ( ij = 1; ij <= j; ij++ )
  {
    for ( ik = 1; ik <= i; ik++ )
    {
      il = b[ik-1];
      da = x[ik-1+(ij-1)*i] - d[il-1+(ij-1)*k];
      db = da * da;
      f[ik-1] = f[ik-1] + db;
      dev[il-1] = dev[il-1] + db;
    }
  }

  for ( ik = 1; ik <= i; ik++ )
  {
    il = b[ik-1];
    fl = e[il-1];
    if ( 2.0 <= fl )
    {
      f[ik-1] = f[ik-1] * fl / ( fl - 1.0 );
    }
  }
//
//  Examine each observation in turn to see if it should be
//  reassigned to a different cluster.
//
  for ( ; ; )
  {
    iw = 0;

    for ( ik = 1; ik <= i; ik++ )
    {
      il = b[ik-1];
      ir = il;
//
//  If the number of cluster points is less than or equal to the
//  specified minimum, NZ, then bypass this iteration.
//
      if ( nz < e[il-1] )
      {
        fl = e[il-1];
        dc = f[ik-1];

        for ( in = 1; in <= n; in++ )
        {
          if ( in != il )
          {
            fm = e[in-1];
            fm = fm / ( fm + 1.0 );

            de = 0.0;
            for ( ip = 1; ip <= j; ip++ )
            {
              da = x[ik-1+(ip-1)*i] - d[in-1+(ip-1)*k];
              de = de + da * da * fm;
            }

            if ( de < dc )
            {
              dc = de;
              ir = in;
            }
          }
        }
//
//  Reassignment is made here if necessary.
//
        if ( ir != il )
        {
          fq = e[ir-1];
          dev[il-1] = dev[il-1] - f[ik-1];
          dev[ir-1] = dev[ir-1] + dc;
          e[ir-1] = e[ir-1] + 1;
          e[il-1] = e[il-1] - 1;

          for ( is = 1; is <= j; is++ )
          {
            d[il-1+(is-1)*k] = ( d[il-1+(is-1)*k] 
                             * fl - x[ik-1+(is-1)*i] ) / ( fl - 1.0 );
            d[ir-1+(is-1)*k] = ( d[ir-1+(is-1)*k]
                             * fq + x[ik-1+(is-1)*i] ) / ( fq + 1.0 );
          }

          b[ik-1] = ir;

          for ( it = 1; it <= i; it++ )
          {
            ij = b[it-1];

            if ( ij == il || ij == ir )
            {
              f[it-1] = 0.0;
              for ( iu = 1; iu <= j; iu++ )
              {
                da = x[it-1+(iu-1)*i] - d[ij-1+(iu-1)*k];
                f[it-1] = f[it-1] + da * da;
              }
              fl = e[ij-1];
              f[it-1] = f[it-1] * fl / ( fl - 1.0 );
            }
          }
          iw = iw + 1;
        }
      }
    }
//
//  If any reassignments were made on this pass, then do another pass.
//
    if ( iw == 0 )
    {
      break;
    }
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
