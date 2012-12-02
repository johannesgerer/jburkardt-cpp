# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>

using namespace std;

# include "quality.hpp"

int main ( int argc, char *argv[] );
void handle ( char *input_filename, 
  double *sample_routine ( int dim_num, int n, int *seed ) );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TABLE_QUALITY.
//
//  Discussion:
//
//    TABLE_QUALITY determines quality measures for a given set of points.
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
//    Max Gunzburger, John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Local parameters:
//
//    Local, int DIM_NUM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[DIM_NUM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
{ 
  int i;
  char input_filename[81];

  cout << "\n";
  timestamp ( );

  cout << "\n";
  cout << "TABLE_QUALITY\n";
  cout << "  C++ version:\n";
  cout << "  Compute measures of uniform dispersion for a pointset.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Note: this routine assumes that the input pointset\n";
  cout << "  is contained in the unit hypercube.  If this is not\n";
  cout << "  the case, then you must rewrite the routine\n";
  cout << "    SAMPLE_ROUTINE\n";
  cout << "  so that it properly returns sample points in your\n";
  cout << "  region, with a uniform density, or a probability\n";
  cout << "  density of your choosing.\n";
//
//  If the input file was not specified, get it now.
//
  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "TABLE_QUALITY:\n";
    cout << "  Please enter the name of a file to be analyzed.\n";

    cin.getline ( input_filename, sizeof ( input_filename ) );

    handle ( input_filename, sample_hypercube_uniform );

  }
  else 
  {
    for ( i = 1; i < argc; i++ ) 
    {
      handle ( argv[i], sample_hypercube_uniform );
    }
  } 
//
//  Terminate.
//
  cout << "\n";
  cout << "TABLE_QUALITY:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void handle ( char *input_filename, 
  double *sample_routine ( int dim_num, int n, int *seed ) )

//****************************************************************************80
//
//  Purpose:
//
//    HANDLE handles a single file.
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
//    Max Gunzburger, John Burkardt,
//    Uniformity Measures for Point Samples in Hypercubes.
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//  Local parameters:
//
//    Local, int DIM_NUM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[DIM_NUM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
//    Input, double *SAMPLE_ROUTINE, the name of a routine which
//    is used to produce an DIM_NUM by N array of sample points in the region, 
//    of the form:
//      double *sample_routine ( int dim_num, int n, int *seed )
//
{
  int flag;
  double *gamma;
  double gamma_ave;
  double gamma_max;
  double gamma_min;
  double gamma_std;
  int i;
  int n;
  int dim_num;
  int ns = 100000;
  int nt;
  int seed_init = 123456789;
  int *triangle;
  int *triangle_neighbor;
  int triangle_order;
  double *z; 

  dtable_header_read ( input_filename, &dim_num, &n );
// 
//  Read the point set.
//
  z = dtable_data_read ( input_filename, dim_num, n );
//
//  For 2D datasets, compute the Delaunay triangulation.
//
  if ( dim_num == 2 )
  {
    triangle = new int[3*2*n];
    triangle_neighbor = new int[3*2*n];

    flag = dtris2 ( n, z, &nt, triangle, triangle_neighbor );
    cout << "\n";
    cout << "  Triangulated data generates " << nt << " triangles.\n";
  }
  else
  {
    nt = 0;
  }

  cout << "\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\".\n";
  cout << "  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =      " << dim_num   << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n" << flush;
  cout << "\n";

  if ( dim_num == 2 )
  {
    triangle_order = 3;
    cout << "  The minimum angle measure    Alpha = "
         << alpha_measure ( n, z, triangle_order, nt, triangle ) << "\n";
  }
  else
  {
    cout << "  The minimum angle measure    Alpha = (omitted)\n";
  }
  cout << "  Relative spacing deviation    Beta = "
    << beta_measure ( dim_num, n, z ) << "\n";
  cout << "  The regularity measure         Chi = "
    << chi_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";
  cout << "  2nd moment determinant measure   D = "
    << d_measure ( dim_num, n, z, ns, sample_hypercube_uniform,
    seed_init ) << "\n";
  cout << "  The Voronoi energy measure       E = "
    << e_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";
  cout << "  The mesh ratio               Gamma = "
    << gamma_measure ( dim_num, n, z ) << "\n";
  cout << "  The point distribution norm      H = "
    << h_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";
  cout << "  The COV measure             Lambda = "
    << lambda_measure ( dim_num, n, z ) << "\n";
  cout << "  The point distribution ratio    Mu = "
    << mu_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";
  cout << "  The cell volume deviation       Nu = "
    << nu_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";
  if ( dim_num == 2 )
  {
    triangle_order = 2;
    cout << "  The triangle uniformity measure  Q = "
         << q_measure ( n, z, triangle_order, nt, triangle ) << "\n";
  }
  else
  {
    cout << "  The triangle uniformity measure  Q = (omitted)\n";
  }
  cout << "  The Riesz S = 0 energy measure  R0 = "
    << r0_measure ( dim_num, n, z ) << "\n";
  if ( r8mat_in_01 ( dim_num, n, z ) )
  {
    cout << "  Nonintersecting sphere volume    S = " 
      << sphere_measure ( dim_num, n, z ) << "\n";
  }
  else
  {
    cout << "  Nonintersecting sphere volume    S = (omitted)\n";
  }
  cout << "  2nd moment trace measure       Tau = "
    << tau_measure ( dim_num, n, z, ns, sample_hypercube_uniform, 
    seed_init ) << "\n";

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

  cout << "\n";
  cout << "  Minimum spacing          Gamma_min = " << gamma_min << "\n";
  cout << "  Average spacing          Gamma_ave = " << gamma_ave << "\n";
  cout << "  Maximum spacing          Gamma_max = " << gamma_max << "\n";
  cout << "  Spacing standard deviate Gamma_std = " << gamma_std << "\n";

  delete [] gamma;
  if ( dim_num == 2 )
  {
    delete [] triangle;
    delete [] triangle_neighbor;
  }
  delete [] z;

  return;
}
