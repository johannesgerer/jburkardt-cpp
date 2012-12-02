# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>
# include <cstring>

using namespace std;

# include "quality.hpp"

int main ( );
void test_cvt ( );
void test_halton ( );
void test_sphere ( );
void test005 ( int n, double z[], int triangle_num,
  int triangle_node[] );
void test006 ( int n, double z[], int triangle_num,
  int triangle_node[] );
void test007 ( int dim_num, int n, double z[] );
void test01 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test02 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test03 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test04 ( int dim_num, int n, double z[] );
void test05 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test06 ( int dim_num, int n, double z[] );
void test07 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test08 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test083 ( int n, double z[], int triangle_num, int triangle_node[] );
void test085 ( int dim_num, int n, double z[] );
void test09 ( int dim_num, int n, double z[] );
void test10 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init );
void test11 ( int dim_num, int n, double z[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for QUALITY_PRB.
//
//  Discussion:
//
//    QUALITY_PRB calls the QUALITY routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  char input_filename[80];
  int n;
  int dim_num;
  int ns;
  int seed_init;
  double *z;

  timestamp ( );

  cout << "\n";
  cout << "QUALITY_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the QUALITY library.\n";

  test_cvt ( );

  test_halton ( );

  test_sphere ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "QUALITY_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_cvt ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_CVT carries out tests of a pointset in the unit hypercube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int flag;
  char input_filename[80];
  int n;
  int dim_num;
  int ns;
  int seed_init;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  double *z;

  cout << "\n";
  cout << "TEST_HALTON\n";
  cout << "  Analyze a pointset in the unit hypercube.\n";
  cout << "\n";
  cout << "  We use a built-in sample routine.\n";

  ns = 100000;
  seed_init = 123456789;
  strcpy ( input_filename, "cvt_02_00100.txt" );

  dtable_header_read ( input_filename, &dim_num, &n );

  cout << "\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\"\n";
  cout << "  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =      " << dim_num      << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n";
  cout << "\n";

  z = dtable_data_read ( input_filename, dim_num, n );

  r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5,
    "  5x5 portion of data read from file:" );
//
//  For 2D datasets, compute the Delaunay triangulation.
//
  if ( dim_num == 2 )
  {
    triangle_node = new int[3*2*n];
    triangle_neighbor = new int[3*2*n];

    flag = dtris2 ( n, z, &triangle_num, triangle_node, triangle_neighbor );
    cout << "\n";
    cout << "  Triangulated data generates " << triangle_num << " triangles.\n";
  }
  else
  {
    triangle_num = 0;
  }

  if ( dim_num == 2 )
  {
    test005 ( n, z, triangle_num, triangle_node );
    test006 ( n, z, triangle_num, triangle_node );
  }

  test007 ( dim_num, n, z );
  test01 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test02 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test03 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test04 ( dim_num, n, z );
  test05 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test06 ( dim_num, n, z );
  test07 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test08 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );

  if ( dim_num == 2 )
  {
    test083 ( n, z, triangle_num, triangle_node );
  }

  test085 ( dim_num, n, z );
  test09 ( dim_num, n, z );
  test10 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test11 ( dim_num, n, z );

  if ( dim_num == 2 )
  {
    delete [] triangle_node;
    delete [] triangle_neighbor;
  }
  delete [] z;

  return;
}
//****************************************************************************80

void test_halton ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_HALTON carries out tests of a pointset in the unit hypercube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int flag;
  char input_filename[80];
  int n;
  int dim_num;
  int ns;
  int seed_init;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  double *z;

  cout << "\n";
  cout << "TEST_HALTON\n";
  cout << "  Analyze a pointset in the unit hypercube.\n";
  cout << "\n";
  cout << "  We use a built-in sample routine.\n";

  ns = 100000;
  seed_init = 123456789;
  strcpy ( input_filename, "halton_02_00100.txt" );

  dtable_header_read ( input_filename, &dim_num, &n );

  cout << "\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\"\n";
  cout << "  The sample routine will be SAMPLE_HYPERCUBE_UNIFORM.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =      " << dim_num      << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n";
  cout << "\n";

  z = dtable_data_read ( input_filename, dim_num, n );

  r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5,
    "  5x5 portion of data read from file:" );
//
//  For 2D datasets, compute the Delaunay triangulation.
//
  if ( dim_num == 2 )
  {
    triangle_node = new int[3*2*n];
    triangle_neighbor = new int[3*2*n];

    flag = dtris2 ( n, z, &triangle_num, triangle_node, triangle_neighbor );
    cout << "\n";
    cout << "  Triangulated data generates " << triangle_num << " triangles.\n";
  }
  else
  {
    triangle_num = 0;
  }

  if ( dim_num == 2 )
  {
    test005 ( n, z, triangle_num, triangle_node );
    test005 ( n, z, triangle_num, triangle_node );
  }

  test007 ( dim_num, n, z );
  test01 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test02 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test03 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test04 ( dim_num, n, z );
  test05 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test06 ( dim_num, n, z );
  test07 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test08 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );

  if ( dim_num == 2 )
  {
    test083 ( n, z, triangle_num, triangle_node );
  }

  test085 ( dim_num, n, z );
  test09 ( dim_num, n, z );
  test10 ( dim_num, n, z, ns, sample_hypercube_uniform, seed_init );
  test11 ( dim_num, n, z );

  if ( dim_num == 2 )
  {
    delete [] triangle_node;
    delete [] triangle_neighbor;
  }
  delete [] z;

  return;
}
//****************************************************************************80

void test_sphere ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_SPHERE carries out tests of a pointset in the unit sphere.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int flag;
  char input_filename[80];
  int n;
  int dim_num;
  int ns;
  int seed_init;
  int *triangle_neighbor;
  int *triangle_node;
  int triangle_num;
  double *z;

  cout << "\n";
  cout << "TEST_SPHERE\n";
  cout << "  Analyze a pointset in the unit sphere.\n";
  cout << "\n";
  cout << "  We use a built-in sample routine.\n";

  ns = 100000;
  seed_init = 123456789;
  strcpy ( input_filename, "sphere_02_00100.txt" );

  dtable_header_read ( input_filename, &dim_num, &n );

  cout << "\n";
  cout << "  Measures of uniform point dispersion.\n";
  cout << "\n";
  cout << "  The pointset was read from \"" << input_filename << "\"\n";
  cout << "  The sample routine will be SAMPLE_SPHERE_UNIFORM.\n";
  cout << "\n";
  cout << "  The spatial dimension DIM_NUM =     " << dim_num      << "\n";
  cout << "  The number of points N =         " << n         << "\n";
  cout << "  The number of sample points NS = " << ns        << "\n";
  cout << "  The random number SEED_INIT =    " << seed_init << "\n";
  cout << "\n";

  z = dtable_data_read ( input_filename, dim_num, n );

  r8mat_transpose_print_some ( dim_num, n, z, 1, 1, 5, 5,
    "  5x5 portion of data read from file:" );
//
//  For 2D datasets, compute the Delaunay triangulation.
//
  if ( dim_num == 2 )
  {
    triangle_node = new int[3*2*n];
    triangle_neighbor = new int[3*2*n];

    flag = dtris2 ( n, z, &triangle_num, triangle_node, triangle_neighbor );
    cout << "\n";
    cout << "  Triangulated data generates " << triangle_num << " triangles.\n";
  }
  else
  {
    triangle_num = 0;
  }

  if ( dim_num == 2 )
  {
    test005 ( n, z, triangle_num, triangle_node );
    test006 ( n, z, triangle_num, triangle_node );
  }

  test007 ( dim_num, n, z );
  test01 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test02 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test03 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test04 ( dim_num, n, z );
  test05 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test06 ( dim_num, n, z );
  test07 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test08 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );

  if ( dim_num == 2 )
  {
    test083 ( n, z, triangle_num, triangle_node );
  }

  test085 ( dim_num, n, z );
  test09 ( dim_num, n, z );
  test10 ( dim_num, n, z, ns, sample_sphere_uniform, seed_init );
  test11 ( dim_num, n, z );

  if ( dim_num == 2 )
  {
    delete [] triangle_node;
    delete [] triangle_neighbor;
  }
  delete [] z;

  return;
}
//****************************************************************************80

void test005 ( int n, double z[], int nt, int triangle[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests ALPHA_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  ALPHA_MEASURE computes the ALPHA measure of quality.\n";
  cout << "  The minimum angle measure    ALPHA = "
       << alpha_measure ( n, z, triangle_order, nt, triangle ) << "\n";

  return;
}
//****************************************************************************80

void test006 ( int n, double z[], int nt, int triangle[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests AREA_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  AREA_MEASURE computes the AREA measure of quality.\n";
  cout << "  The area ratio measure        AREA = "
       << area_measure ( n, z, triangle_order, nt, triangle ) << "\n";

  return;
}
//****************************************************************************80

void test007 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests BETA_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST007\n";
  cout << "  BETA_MEASURE computes the BETA measure of quality.\n";
  cout << "  Relative spacing deviation BETA =    "
       << beta_measure ( dim_num, n, z ) << "\n";

  return;
}
//****************************************************************************80

void test01 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests CHI_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST01\n";
  cout << "  CHI_MEASURE computes the CHI measure of quality.\n";
  cout << "  The regularity measure         Chi = "
       << chi_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test02 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests D_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST02\n";
  cout << "  D_MEASURE computes the D measure of quality.\n";
  cout << "  2nd moment determinant measure   D = "
       << d_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test03 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests E_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST03\n";
  cout << "  E_MEASURE computes the E measure of quality.\n";
  cout << "  The Voronoi energy measure       E = "
       << e_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test04 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests GAMMA_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST04\n";
  cout << "  GAMMA_MEASURE computes the Gamma measure of quality.\n";
  cout << "  The mesh ratio               Gamma = "
       << gamma_measure ( dim_num, n, z ) << "\n";

  return;
}
//****************************************************************************80

void test05 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests H_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST05\n";
  cout << "  H_MEASURE computes the H measure of quality.\n";
  cout << "  The point distribution norm      H = "
       << h_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test06 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests LAMBDA_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST06\n";
  cout << "  LAMBDA_MEASURE computes the Lambda measure of quality.\n";
  cout << "  The COV measure             Lambda = "
       << lambda_measure ( dim_num, n, z ) << "\n";

  return;
}
//****************************************************************************80

void test07 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests MU_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST07\n";
  cout << "  MU_MEASURE computes the Mu measure of quality.\n";
  cout << "  The point distribution ratio    Mu = "
       << mu_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test08 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests NU_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST08\n";
  cout << "  NU_MEASURE computes the Nu measure of quality.\n";
  cout << "  The cell volume deviation       Nu = "
       << nu_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test083 ( int n, double z[], int nt, int triangle[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST083 tests Q_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 November 2005
//
//  Author:
//
//    John Burkardt
//
{
  int triangle_order = 3;

  cout << "\n";
  cout << "TEST083\n";
  cout << "  Q_MEASURE computes the Q measure of quality.\n";
  cout << "  The triangle uniformity measure  Q = "
       << q_measure ( n, z, triangle_order, nt, triangle ) << "\n";

  return;
}
//****************************************************************************80

void test085 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests R0_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST085\n";
  cout << "  R0_MEASURE computes the Riesz S = 0 energy measure of quality.\n";
  cout << "  The R0 measure                  R0 = "
       << r0_measure ( dim_num, n, z ) << "\n";

  return;
}
//****************************************************************************80

void test09 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests SPHERE_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST09\n";
  cout << "  SPHERE_MEASURE computes the sphere measure of quality.\n";
  cout << "  Nonintersecting sphere volume    S = "
       << sphere_measure ( dim_num, n, z ) << "\n";

  return;
}
//****************************************************************************80

void test10 ( int dim_num, int n, double z[], int ns,
  double *sample_routine ( int dim_num, int n, int *seed ), int seed_init )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests TAU_MEASURE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  cout << "\n";
  cout << "TEST10\n";
  cout << "  TAU_MEASURE computes the Tau measure of quality.\n";
  cout << "  2nd moment trace measure       Tau = "
       << tau_measure ( dim_num, n, z, ns, sample_routine, seed_init ) << "\n";

  return;
}
//****************************************************************************80

void test11 ( int dim_num, int n, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests POINTSET_SPACING.
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
{
  double *gamma;
  double gamma_ave;
  double gamma_max;
  double gamma_min;
  double gamma_std;
  int i;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  POINTSET_SPACING computes pointset spacing parameters.\n";

  gamma = pointset_spacing ( dim_num, n, z );

  gamma_min = r8vec_min ( n, gamma );
  gamma_max = r8vec_max ( n, gamma );

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
  cout << "  Minimum spacing          GAMMA_MIN = " << gamma_min << "\n";
  cout << "  Average spacing          GAMMA_AVE = " << gamma_ave << "\n";
  cout << "  Maximum spacing          GAMMA_MAX = " << gamma_max << "\n";
  cout << "  Spacing standard dev     GAMMA_STD = " << gamma_std << "\n";

  delete [] gamma;

  return;
}
