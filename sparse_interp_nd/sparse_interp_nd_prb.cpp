# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <cmath>
# include <ctime>

using namespace std;

# include "sparse_interp_nd.hpp"
# include "r8lib.hpp"

int main ( );
void test01 ( int m, int sparse_max );
double *f_sinr ( int m, int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    SPARSE_INTERP_ND_PRB tests SPARSE_INTERP_ND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int sparse_max;

  timestamp ( );
  cout << " \n";
  cout << "SPARSE_INTERP_ND_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the SPARSE_INTERP_ND library.\n";
  cout << "  The R8LIB library is also required.\n";

  m = 1;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 2;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 3;
  sparse_max = 9;
  test01 ( m, sparse_max );

  m = 4;
  sparse_max = 7;
  test01 ( m, sparse_max );
//
//  Terminate.
//
  cout << " \n";
  cout << "SPARSE_INTERP_ND_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << " \n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int m, int sparse_max )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01: sequence of sparse interpolants to an M-dimensional function.
//
//  Discussion:
//
//    We have functions that can generate a Lagrange interpolant to data
//    in M dimensions, with specified order or level in each dimension.
//
//    We use the Lagrange function as the inner evaluator for a sparse
//    grid procedure. 
//
//    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
//    to a given function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Local, int M, the spatial dimension.
//
//    Input, int SPARSE_MAX, the maximum sparse grid level to try.
//
//  Local Parameters:
//
//    Local, double A[M], B[M], the upper and lower variable limits 
//    in each dimension.
//
//    Local, double APP_ERROR, the averaged Euclidean norm of the 
//    difference between the sparse interpolant and the exact function at 
//    the interpolation points.
//
//    Local, int C[L_MAX+1], the sparse grid coefficient vector.
//    Results at level L are weighted by C(L).
//
//    Local, int IND[M], the 1D indices defining a Lagrange grid.
//    Each entry is a 1d "level" that specifies the order of a 
//    Clenshaw Curtis 1D grid.
//
//    Local, int L, the current Lagrange grid level.
//
//    Local, int L_MAX, the current sparse grid level.
//
//    Local, int MORE, used to control the enumeration of all the
//    Lagrange grids at a current grid level.
//
//    Local, int ND, the number of points used in a Lagrange grid.
//
//    Local, int ND_TOTAL, the total number of points used in all the
//    Lagrange interpolants for a given sparse interpolant points that occur
//    in more than one Lagrange grid are counted multiple times.
//
//    Local, int NI, the number of interpolant evaluation points.
//
//    Local, int SPARSE_MIN, the minimum sparse grid level to try.
//
//    Local, double XD[M*ND], the data points for a Lagrange grid.
//
//    Local, double XI[M*NI], the interpolant evaluation points.
//
//    Local, double ZD[ND], the data values for a Lagrange grid.
//
//    Local, double ZE[NI], the exact function values at XI.
//
//    Local, double ZI[NI], the sparse interpolant values at XI.
//
//    Local, double ZPI[NI], one set of Lagrange interpolant values at XI.
//
{
  double *a;
  double app_error;
  double *b;
  int *c;
  int h;
  int i;
  int *ind;
  int l;
  int l_max;
  int l_min;
  int more;
  int nd;
  int nd_total;
  int ni;
  int seed;
  int sparse_min;
  int t;
  int *w;
  double *xd;
  double *xi;
  double *zd;
  double *ze;
  double *zi;
  double *zpi;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  Sparse interpolation for a function f(x) of M-dimensional argument.\n";
  cout << "  Use a sequence of sparse grids of levels 0 through SPARSE_MAX.\n";
  cout << "  Invoke a general Lagrange interpolant function to do this.\n";
  cout << "\n";
  cout << "  Compare the exact function and the interpolants at a grid of points.\n";
  cout << "\n";
  cout << "  The \"order\" is the sum of the orders of all the product grids\n";
  cout << "  used to make a particular sparse grid.\n";
//
//  User input.
//
  cout << "\n";
  cout << "  Spatial dimension M = " << m << "\n";
  cout << "  Maximum sparse grid level = " << sparse_max << "\n";
//
//  Define the region.
//
  a = new double[m];
  b = new double[m];

  for ( i = 0; i < m; i++ )
  {
    a[i] = 0.0;
    b[i] = 1.0;
  }
//
//  Define the interpolation evaluation information.
//
  ni = 100;
  seed = 123456789;
  xi = r8mat_uniform_abvec_new ( m, ni, a, b, seed );

  cout << "  Number of interpolation points is NI = " << ni << "\n";

  ze = f_sinr ( m, ni, xi );
//
//  Compute a sequence of sparse grid interpolants of increasing level.
//
  cout << "\n";
  cout << "   L     Order    ApproxError\n";
  cout << "\n";

  ind = new int[m];
  zi = new double[ni];

  sparse_min = 0;

  for ( l_max = sparse_min; l_max <= sparse_max; l_max++ )
  {
    c = new int[l_max+1];
    w = new int[l_max+1];
    smolyak_coefficients ( l_max, m, c, w );
    
    for ( i = 0; i < ni; i++ )
    {
      zi[i] = 0.0;
    }
    nd_total = 0;

    l_min = i4_max ( l_max + 1 - m, 0 );

    for ( l = l_min; l <= l_max; l++ )
    {
      more = 0;
      while ( 1 )
      {
//
//  Define the next product grid.
//
        comp_next ( l, m, ind, more, h, t );
//
//  Count the grid, find its points, evaluate the data there.
//
        nd = lagrange_interp_nd_size2 ( m, ind );
        xd = lagrange_interp_nd_grid2 ( m, ind, a, b, nd );
        zd = f_sinr ( m, nd, xd );
//
//  Use the grid to evaluate the interpolant.
//
        zpi = lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, xi );
//
//  Weighted the interpolant values and add to the sparse grid interpolant.
//
        nd_total = nd_total + nd;
        for ( i = 0; i < ni; i++ )
        {
          zi[i] = zi[i] + c[l] * zpi[i];
        }

        delete [] xd;
        delete [] zd;
        delete [] zpi;

        if ( !more )
        {
          break;
        }
      }
    }
//
//  Compare sparse interpolant and exact function at interpolation points.
//
    app_error = r8vec_norm_affine ( ni, zi, ze ) / ( double ) ( ni );

    cout << "  " << setw(2) << l
         << "  " << setw(8) << nd_total
         << "  " << setw(8) << app_error << "\n";

    delete [] c;
    delete [] w;

  }

  delete [] a;
  delete [] b;
  delete [] ind;
  delete [] xi;
  delete [] ze;
  delete [] zi;

  return;
}
//****************************************************************************80

double *f_sinr ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, double X(M,N), the points.
//
//    Output, double F_SINR[N], the value of the function at each point.
//
{
  int i;
  int j;
  double r;
  double *z;

  z = ( double * ) malloc ( n * sizeof ( double ) );

  for ( j = 0; j < n; j++ )
  {
    r = 0.0;
    for ( i = 0; i < m; i++ )
    {
      r = r + x[i+j*m] * x[i+j*m];
    }
    r = sqrt ( r );
    z[j] = sin ( r );
  }

  return z;
}
