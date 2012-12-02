# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <complex>

using namespace std;

# include "kmeans.hpp"

int main ( );

void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );
void test10 ( );
void test11 ( );
void test12 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for KMEANS_PRB.
//
//  Discussion:
//
//    KMEANS_PRB calls the KMEANS tests.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "KMEANS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the KMEANS library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );

  test10 ( );
  test11 ( );
  test12 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "KMEANS_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tries out the HMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test01_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test01_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the HMEANS_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];
  
  hmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, 
    cluster_num, point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tries out the HMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test02_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test02_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Test the HMEANS_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];
  
  hmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy, &seed );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, 
    cluster_num, point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tries out the KMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test03_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test03_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  Test the KMEANS_01 algorithm.\n";
  cout << "  (Applied Statistics Algorithm #58)\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, point, 
    &seed );

  kmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tries out the KMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test04_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test04_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  Test the KMEANS_02 algorithm.\n";
  cout << "  (Applied Statistics Algorithm #136)\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num, point );

  kmeans_02 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tries out the KMEANS_03 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test05_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test05_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  Test the KMEANS_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num,
    point );

  kmeans_03 ( dim_num, point_num, cluster_num, it_max, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tries out the HMEANS_01 + KMEANS_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test06_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test06_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  Test the HMEANS_01 + KMEANS_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;
  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_01 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_01 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num,
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tries out the HMEANS_01 + KMEANS_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test07_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test07_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  Test the HMEANS_01 + KMEANS_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_02 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_02 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_02 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tries out the HMEANS_01 + KMEANS_03 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test08_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test08_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int it_max;
  int it_max_h;
  int it_max_k;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  Test the HMEANS_01 + KMEANS_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;
  it_max_h = 3;
  it_max_k = it_max;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Initialize by using a few steps of HMEANS_02:\n";
  cout << "  Number of HMEANS_01 iterations allowed is " << it_max_h << "\n";
  cout << "  Number of KMEANS_03 iterations allowed is " << it_max_k << "\n";
//
//  Initialize the centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );
//
//  Initialize the clusters.
//
  hmeans_01 ( dim_num, point_num, cluster_num, it_max_h, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of HMEANS_01 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  kmeans_03 ( dim_num, point_num, cluster_num, it_max_k, it_num, point, 
    cluster, cluster_center, cluster_population, cluster_energy );

  cout << "\n";
  cout << "  Number of KMEANS_03 iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tries out the HMEANS_W_01 routine.
//
//  Discussion:
//
//    The weights are all equal, so the results should
//    be identical to those for HMEANS_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test09_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test09_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_equal_100.txt";

  cout << "\n";
  cout << "TEST09\n";
  cout << "  Test the HMEANS_W_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tries out the HMEANS_W_02 routine.
//
//  Discussion:
//
//    The weights are all equal, so the results should
//    be identical to those for HMEANS_02.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test10_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test10_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_equal_100.txt";

  cout << "\n";
  cout << "TEST10\n";
  cout << "  Test the HMEANS_W_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy, &seed );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, cluster_population, 
    cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tries out the KMEANS_W_01 routine.
//
//  Discussion:
//
//   The weights are all equal, so the results should
//    be identical to those for KMEANS_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test11_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test11_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_equal_100.txt";

  cout << "\n";
  cout << "TEST11\n";
  cout << "  Test the KMEANS_W_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tries out the KMEANS_W_03 routine.
//
//  Discussion:
//
//    The weights are all equal, so the results should
//    be identical to those for KMEANS_03.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test12_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test12_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_equal_100.txt";

  cout << "\n";
  cout << "TEST12\n";
  cout << "  Test the KMEANS_W_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num,
    point );

  kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tries out the HMEANS_W_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test13_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test13_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_unequal_100.txt";

  cout << "\n";
  cout << "TEST13\n";
  cout << "  Test the HMEANS_W_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tries out the HMEANS_W_02 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test14_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test14_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_unequal_100.txt";

  cout << "\n";
  cout << "TEST14\n";
  cout << "  Test the HMEANS_W_02 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  hmeans_w_02 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy, &seed );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tries out the KMEANS_W_01 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test15_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test15_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_unequal_100.txt";

  cout << "\n";
  cout << "TEST15\n";
  cout << "  Test the KMEANS_W_01 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_5 ( dim_num, point_num, cluster_num, 
    point, &seed );

  kmeans_w_01 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tries out the KMEANS_W_03 routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 October 2011
//
//  Author:
//
//    John Burkardt
//
{
  string center_filename = "test16_centers.txt";
  int *cluster;
  double *cluster_center;
  double *cluster_energy;
  string cluster_filename = "test16_clusters.txt";
  int cluster_num;
  int *cluster_population;
  double *cluster_variance;
  int dim_num;
  int i;
  int it_max;
  int it_num;
  double *point;
  string point_filename = "points_100.txt";
  int point_num;
  int seed;
  double *weight;
  int weight_dim;
  int weight_num;
  string weight_filename = "weights_unequal_100.txt";

  cout << "\n";
  cout << "TEST16\n";
  cout << "  Test the KMEANS_W_03 algorithm.\n";
//
//  Read the data points.
//
  cout << "\n";
  cout << "  Data points will be read from \"" << point_filename << "\".\n";

  r8mat_header_read ( point_filename, &dim_num, &point_num );

  cout << "\n";
  cout << "  Point spatial dimension = " << dim_num << "\n";
  cout << "  Number of points = " << point_num << "\n";

  point = r8mat_data_read ( point_filename, dim_num, point_num );
//
//  Read the weights.
//
  cout << "\n";
  cout << "  Weights will be read from \"" << weight_filename << "\".\n";

  r8mat_header_read ( weight_filename, &weight_dim, &weight_num );

  if ( weight_dim != 1 )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Spatial dimension of weight array is not 1.\n";
    exit ( 1 );
  }

  if ( weight_num != point_num )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  Number of weights not equal to number of points.\n";
    exit ( 1 );
  }

  weight = r8mat_data_read ( weight_filename, weight_dim, weight_num );
//
//  Clustering parameters.
//
  cluster_num = 5;
  it_max = 20;
  seed = 123456789;

  cluster = i4vec_negone_new ( point_num );
  cluster_energy = new double[cluster_num];
  cluster_population = new int[cluster_num];

  cout << "\n";
  cout << "  Number of iterations allowed is " << it_max << "\n";
//
//  Initialize the cluster centers.
//
  cluster_center = cluster_initialize_1 ( dim_num, point_num, cluster_num, point );

  kmeans_w_03 ( dim_num, point_num, cluster_num, it_max, it_num, 
    point, weight, cluster, cluster_center, cluster_population, 
    cluster_energy );

  cout << "\n";
  cout << "  Number of iterations taken is " << it_num << "\n";

  cluster_variance = cluster_variance_compute ( dim_num, point_num, cluster_num, 
    point, cluster, cluster_center );

  cluster_print_summary ( point_num, cluster_num, 
    cluster_population, cluster_energy, cluster_variance );

  r8mat_write ( center_filename, dim_num, cluster_num, cluster_center );

  cout << "\n";
  cout << "  Cluster centers written to \"" << center_filename << "\".\n";

  i4mat_write ( cluster_filename, 1, point_num, cluster );

  cout << "  Cluster assignments written to \"" << cluster_filename << "\".\n";

  delete [] cluster;
  delete [] cluster_center;
  delete [] cluster_energy;
  delete [] cluster_population;
  delete [] cluster_variance;
  delete [] point;
  delete [] weight;

  return;
}

