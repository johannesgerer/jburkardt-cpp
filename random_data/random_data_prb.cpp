# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "random_data.hpp"

int main ( );
void test005 ( );
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
void test115 ( );
void test12 ( );
void test125 ( );
void test13 ( );
void test14 ( );
void test15 ( );
void test16 ( );
void test17 ( );
void test18 ( );
void test19 ( );
void test20 ( );
void test205 ( );
void test21 ( );
void test22 ( );
void test23 ( );
void test24 ( );
void test245 ( );
void test25 ( );
void test26 ( );
void test264 ( );
void test265 ( );
void test267 ( );
void test27 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    RANDOM_DATA_PRB calls the RANDOM_DATA tests.
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
{
  timestamp ( );
  cout << "\n";
  cout << "RANDOM_DATA_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the RANDOM_DATA library.\n";

  test005 ( );
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
  test115 ( );
  test12 ( );
  test125 ( );
  test13 ( );
  test14 ( );
  test15 ( );
  test16 ( );
  test17 ( );
  test18 ( );
  test19 ( );

  test20 ( );
  test205 ( );
  test21 ( );
  test22 ( );
  test23 ( );
  test24 ( );
  test245 ( );
  test25 ( );
  test26 ( );
  test264 ( );
  test265 ( );
  test267 ( );
  test27 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "RANDOM_DATA_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests BROWNIAN
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
{
  int dim_num;
  int n;
  string output_filename;
  int seed;
  double *x;

  cout << "\n";
  cout << "TEST005:\n";
  cout << "  BAD_IN_SIMPLEX01 is a \"bad\" sampling technique\n";
  cout << "  for the unit simplex.\n";


  for ( dim_num = 2; dim_num <= 3; dim_num++ )
  {
    seed = 123456789;
    n = 10000;
    if ( dim_num == 2 )
    {
      output_filename = "bad_in_triangle.txt";
    }
    else if ( dim_num == 3 )
    {
      output_filename = "bad_in_tetrahedron.txt";
    }

    cout << "\n";
    cout << "  Spatial dimension DIM_NUM =  " << dim_num      << "\n";
    cout << "  Number of points N =         " << n      << "\n";
    cout << "  Initial random number SEED = " << seed   << "\n";

    x = bad_in_simplex01 ( dim_num, n, &seed );

    r8mat_write ( output_filename, dim_num, n, x );

    cout << "\n";
    cout << "  Data written to file \"" << output_filename << "\".\n";

    delete [] x;
  }

  return;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests BROWNIAN
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 100

  string output_filename = "brownian.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST01:\n";
  cout << "  BROWNIAN generates Brownian motion points.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = brownian ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests D_NORMAL_01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed = 123456789;
  int seed_in;
  double x;

  cout << "\n";
  cout << "TEST02:\n";
  cout << "  D_NORMAL_01 generates a single normal \n";
  cout << "  pseudorandom value.\n";
  cout << "\n";
  cout << "     Seed          Seed       D_NORMAL_01\n";
  cout << "    (Input)       (Output)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_normal_01 ( &seed );

    cout                        << "  "
         << setw(12) << seed_in << "  "
         << setw(12) << seed    << "  "
         << setw(12) << x       << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests R8_UNIFORM_01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed = 123456789;
  int seed_in;
  double x;

  cout << "\n";
  cout << "TEST03:\n";
  cout << "  D_UNIFORM_01 generates a single uniform \n";
  cout << "  pseudorandom value.\n";
  cout << "\n";
  cout << "     Seed          Seed       D_UNIFORM_01\n";
  cout << "    (Input)       (Output)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    seed_in = seed;
    x = r8_uniform_01 ( &seed );

    cout                        << "  "
         << setw(12) << seed_in << "  "
         << setw(12) << seed    << "  "
         << setw(12) << x       << "\n";
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests GRID_IN_CUBE01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 85

  int center = 1;
  string output_filename = "grid_in_cube01.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  GRID_IN_CUBE01 generates grid points in the unit hypercube.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  CENTER option =              " << center << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = grid_in_cube01 ( DIM_NUM, N, center, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests HALTON_IN_CIRCLE01_ACCEPT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 400

  string output_filename = "halton_in_circle01_accept.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  HALTON_IN_CIRCLE01_ACCEPT accepts \n";
  cout << "  Halton points in the unit circle.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = halton_in_circle01_accept ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests HALTON_IN_CIRCLE01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 400

  string output_filename = "halton_in_circle01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  HALTON_IN_CIRCLE01_MAP maps \n";
  cout << "  Halton points into the unit circle.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = halton_in_circle01_map ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests HALTON_IN_CUBE01
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 100

  string output_filename = "halton_in_cube01.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  HALTON_IN_CUBE01 generates Halton points\n";
  cout << "  in the unit hypercube.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = halton_in_cube01 ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests HAMMERSLEY_IN_CUBE01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 100

  string output_filename = "hammersley_in_cube01.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  HAMMERSLEY_IN_CUBE01 generates Hammersley points\n";
  cout << "  in the unit hypercube.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = hammersley_in_cube01 ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests NORMAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "normal.txt";
  int i;
  int info;
  int j;
  double mu[DIM_NUM] = { 6.0, 100.0 };
  double r[DIM_NUM*DIM_NUM];
  int seed = 123456789;
  double v[DIM_NUM*DIM_NUM] = { 5.0, 2.0, 2.0, 1.0 };
  double *x;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  NORMAL generates normal points\n";
  cout << "    in M dimensions, using a nonzero mean, and with\n";
  cout << "    user-specified variance-covariance matrix.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  r8vec_print ( DIM_NUM, mu, "  Mean vector MU:" );

  r8mat_print ( DIM_NUM, DIM_NUM, v, "  Variance-covariance matrix V:" );

  for ( i = 0; i < DIM_NUM; i++ )
  {
    for ( j = 0; j < DIM_NUM; j++ )
    {
      r[i+j*DIM_NUM] = v[i+j*DIM_NUM];
    }
  }

  info = dpofa ( r, DIM_NUM, DIM_NUM );

  if ( info != 0 )
  {
    cout << "\n";
    cout << "TEST04 - Fatal error!\n";
    cout << "  Variance-covariance matrix factorization failed.\n";
    exit ( 1 );
  }

  r8mat_print ( DIM_NUM, DIM_NUM, r, "  Cholesky factor R:" );

  x = normal ( DIM_NUM, N, r, mu, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests NORMAL_CIRCULAR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 2000

  string output_filename = "normal_circular.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  NORMAL_CIRCULAR generates points in 2D\n";
  cout << "    distributed according to a circular normal.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = normal_circular ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests NORMAL_SIMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "normal_simple.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST11\n";
  cout << "  NORMAL_SIMPLE generates normal points\n";
  cout << "    in M dimensions, using a zero mean, and with\n";
  cout << "    the identity as the variance-covariance matrix.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = normal_simple ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests UNIFORM_IN_ANNULUS.
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
{
# define DIM_NUM 2
# define N 400

  string output_filename = "uniform_in_annulus.txt";
  double pc[DIM_NUM] = { 10.0, 5.0 };
  double r1 = 1.0;
  double r2 = 3.0;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST115\n";
  cout << "  UNIFORM_IN_ANNULUS generates uniform\n";
  cout << "  points in an annulus by mapping.\n";
  cout << "\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Center PC(1:2) =             " << pc[0]  << "  " << pc[1] << "\n";
  cout << "  Inner radius is R1 =         " << r1     << "\n";
  cout << "  Outer radius is R2 =         " << r2     << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_annulus ( pc, r1, r2, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests UNIFORM_IN_ANNULUS_ACCEPT.
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
{
# define DIM_NUM 2
# define N 400

  string output_filename = "uniform_in_annulus_accept.txt";
  double pc[DIM_NUM] = { 10.0, 5.0 };
  double r1 = 1.0;
  double r2 = 3.0;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  UNIFORM_IN_ANNULUS_ACCEPT generates uniform\n";
  cout << "  points in an annulus by acceptance/rejection.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Center PC(1:2) =             " << pc[0]  << "  " << pc[1] << "\n";
  cout << "  Inner radius is R1 =         " << r1     << "\n";
  cout << "  Outer radius is R2 =         " << r2     << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_annulus_accept ( pc, r1, r2, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests UNIFORM_IN_ANNULUS_SECTOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 400

  string output_filename = "uniform_in_annulus_sector.txt";
  double pc[DIM_NUM] = { 10.0, 5.0 };
  double r1 = 1.0;
  double r2 = 3.0;
  int seed = 123456789;
  double theta1 = 0.0;
  double theta2 = 1.5707964;
  double *x;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  UNIFORM_IN_ANNULUS_SECTOR generates uniform \n";
  cout << "  points in an annular sector by mapping.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Center PC(1:2) =             " << pc[0]  << "  " << pc[1] << "\n";
  cout << "  Inner radius is R1 =         " << r1     << "\n";
  cout << "  Outer radius is R2 =         " << r2     << "\n";
  cout << "  THETA1 =                     " << theta1 << "\n";
  cout << "  THETA2 =                     " << theta2 << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_annulus_sector ( pc, r1, r2, theta1, theta2, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test13 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST13 tests UNIFORM_IN_CIRCLE01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 400

  string output_filename = "uniform_in_circle01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  UNIFORM_IN_CIRCLE01_MAP maps uniform \n";
  cout << "  points into the unit circle.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_circle01_map ( N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test14 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST14 tests UNIFORM_IN_CUBE01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_cube01.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  UNIFORM_IN_CUBE01 generates uniform\n";
  cout << "  points in the unit hypercube.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_cube01 ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test15 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST15 tests UNIFORM_IN_ELLIPSOID_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  double a[DIM_NUM*DIM_NUM] = { 3.0, 1.0, 1.0, 2.0 };
  int fail_num;
  string output_filename = "uniform_in_ellipsoid_map.txt";
  int i;
  int j;
  int k;
  double r = 1.0;
  double r2;
  int seed = 123456789;
  int success_num;
  double *x;

  cout << "\n";
  cout << "TEST15\n";
  cout << "  UNIFORM_IN_ELLIPSOID_MAP maps uniform\n";
  cout << "  points into an ellipsoid.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_ellipsoid_map ( DIM_NUM, N, a, r, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";
//
//  Test the data.
//
  fail_num = 0;
  success_num = 0;

  for ( j = 0; j < N; j++ )
  {

    r2 = 0.0;
    for ( i = 0; i < DIM_NUM; i++ )
    {
      for ( k = 0; k < DIM_NUM; k++ )
      {
        r2 = r2 + x[i+j*DIM_NUM] * a[i+k*DIM_NUM] * x[k+j*DIM_NUM];
      }
    }
    r2 = sqrt ( r2 );

    if ( r < r2 )
    {
      fail_num = fail_num + 1;
    }
    else
    {
      success_num = success_num + 1;
    }
  }

  cout << "\n";
  cout << "  " << fail_num
       << " points failed the ellipsoid test.\n";

  cout << "  " << success_num
       << " points satisfy the ellipsoid test.\n";

  delete [] x;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test16 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST16 tests UNIFORM_IN_PARALLELOGRAM_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_parallelogram_map.txt";
  int seed = 123456789;
  double v1[DIM_NUM] = { 0.75E+00, 0.90E+00 };
  double v2[DIM_NUM] = { 0.00E+00, 0.20E+00 };
  double v3[DIM_NUM] = { 1.10E+00, 0.65E+00 };
  double v4[DIM_NUM];
  double *x;

  cout << "\n";
  cout << "TEST16\n";
  cout << "  UNIFORM_IN_PARALLELOGRAM_MAP maps uniform\n";
  cout << "  points into a parallelogram.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  v4[0] = v3[0] + v2[0] - v1[0];
  v4[1] = v3[1] + v2[1] - v1[1];

  cout << "\n";
  cout << "  V1 = " << v1[0] << ", " << v1[1] << "\n";
  cout << "  V2 = " << v2[0] << ", " << v2[1] << "\n";
  cout << "  V3 = " << v3[0] << ", " << v3[1] << "\n";
  cout << "  V4 = " << v4[0] << ", " << v4[1] << "\n";

  x = uniform_in_parallelogram_map ( v1, v2, v3, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test17 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST17 tests UNIFORM_IN_POLYGON_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000
# define NV 10

  string output_filename = "uniform_in_polygon_map.txt";
  int seed = 123456789;
  double v[DIM_NUM*NV] = {
    0.0E+00, 0.0E+00,
    0.5E+00, 0.3E+00,
    1.0E+00, 0.0E+00,
    0.7E+00, 0.4E+00,
    1.0E+00, 0.6E+00,
    0.6E+00, 0.6E+00,
    0.5E+00, 1.0E+00,
    0.4E+00, 0.6E+00,
    0.0E+00, 0.6E+00,
    0.3E+00, 0.4E+00 };
  double *x;

  cout << "\n";
  cout << "TEST17\n";
  cout << "  UNIFORM_IN_POLYGON_MAP maps uniform\n";
  cout << "  points into a polygon.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  r8mat_print ( DIM_NUM, NV, v, "  Polygonal vertices:" );

  x = uniform_in_polygon_map ( NV, v, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( "polygon_vertices.txt", DIM_NUM, NV, v );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
# undef NV
}
//****************************************************************************80

void test18 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST18 tests UNIFORM_IN_SECTOR_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 300

  string output_filename = "uniform_in_sector_map.txt";
  double r1 = 1.0;
  double r2 = 2.0;
  int seed = 123456789;
  double t1 = 0.78;
  double t2 = 2.35;
  double *x;

  cout << "\n";
  cout << "TEST18\n";
  cout << "  UNIFORM_IN_SECTOR_MAP maps uniform\n";
  cout << "  points into a circular sector.\n";
  cout << "\n";
  cout << "  R1 = " << r1 << "\n";
  cout << "  R2 = " << r2 << "\n";
  cout << "  T1 = " << t1 << "\n";
  cout << "  T2 = " << t2 << "\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_sector_map ( r1, r2, t1, t2, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test19 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST19 tests UNIFORM_IN_SIMPLEX01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_simplex01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST19\n";
  cout << "  UNIFORM_IN_SIMPLEX01_MAP maps uniform\n";
  cout << "  points into the unit simplex\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_simplex01_map ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test20 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST20 tests UNIFORM_IN_SPHERE01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_sphere01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST20\n";
  cout << "  UNIFORM_IN_SPHERE01_MAP maps uniform\n";
  cout << "  points into the unit sphere.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_sphere01_map ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test205 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST205 tests UNIFORM_IN_TETRAHEDRON.
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
{
# define DIM_NUM 3
# define N 1000

  string output_filename = "uniform_in_tetrahedron.txt";
  int seed = 123456789;
  double v[3*4] = {
    1.0,  2.0,  3.0,
    4.0,  1.0,  2.0,
    2.0,  4.0,  4.0,
    3.0,  2.0,  5.0  };
  double *x;

  cout << "\n";
  cout << "TEST205\n";
  cout << "  UNIFORM_IN_TETRAHEDRON returns uniform\n";
  cout << "  points from a tetrahedron.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  r8mat_print ( 3, 4, v, "  Tetrahedron vertices:" );

  x = uniform_in_tetrahedron ( v, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test21 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST21 tests UNIFORM_IN_TRIANGLE_MAP1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_triangle_map1.txt";
  int seed = 123456789;
  double v1[DIM_NUM] = { 0.75E+00, 0.90E+00 };
  double v2[DIM_NUM] = { 0.00E+00, 0.20E+00 };
  double v3[DIM_NUM] = { 0.95E+00, 0.65E+00 };
  double *x;

  cout << "\n";
  cout << "TEST21\n";
  cout << "  UNIFORM_IN_TRIANGLE_MAP1 maps uniform\n";
  cout << "  points into a triangle, by Turk 1 mapping.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  cout << "\n";
  cout << "  V1 = " << v1[0] << ", " << v1[1] << "\n";
  cout << "  V2 = " << v2[0] << ", " << v2[1] << "\n";
  cout << "  V3 = " << v3[0] << ", " << v3[1] << "\n";

  x = uniform_in_triangle_map1 ( v1, v2, v3, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test22 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST22 tests UNIFORM_IN_TRIANGLE_MAP2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 1000

  string output_filename = "uniform_in_triangle_map2.txt";
  int seed = 123456789;
  double v1[DIM_NUM] = { 0.75E+00, 0.90E+00 };
  double v2[DIM_NUM] = { 0.00E+00, 0.20E+00 };
  double v3[DIM_NUM] = { 0.95E+00, 0.65E+00 };
  double *x;

  cout << "\n";
  cout << "TEST22\n";
  cout << "  UNIFORM_IN_TRIANGLE_MAP2 maps uniform\n";
  cout << "  points into a triangle, by Turk 2 mapping.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  cout << "\n";
  cout << "  V1 = " << v1[0] << ", " << v1[1] << "\n";
  cout << "  V2 = " << v2[0] << ", " << v2[1] << "\n";
  cout << "  V3 = " << v3[0] << ", " << v3[1] << "\n";

  x = uniform_in_triangle_map2 ( v1, v2, v3, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test23 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST23 tests UNIFORM_IN_TRIANGLE01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 2000

  string output_filename = "uniform_in_triangle01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST23\n";
  cout << "  UNIFORM_IN_TRIANGLE01_MAP maps uniform\n";
  cout << "  points into the unit triangle.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_in_triangle01_map ( N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test24 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST24 tests UNIFORM_ON_ELLIPSOID_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 200

  double a[DIM_NUM*DIM_NUM] = { 3.0, 1.0, 1.0, 2.0 };
  string output_filename = "uniform_on_ellipsoid_map.txt";
  double r = 1.0;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST24\n";
  cout << "  UNIFORM_ON_ELLIPSOID_MAP maps uniform\n";
  cout << "  points onto an ellipsoid.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_on_ellipsoid_map ( DIM_NUM, N, a, r, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test245 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST245 tests UNIFORM_ON_HEMISPHERE01_PHONG.
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
{
# define DIM_NUM 3
# define N 50

  string output_filename = "uniform_on_hemisphere01_phong.txt";
  int m = 2;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST245\n";
  cout << "  UNIFORM_ON_HEMISPHERE01_PHONG maps uniform\n";
  cout << "  points onto the unit hemisphere with Phong density.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Phong exponent M =           " << m      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_on_hemisphere01_phong ( N, m, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test25 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST25 tests UNIFORM_ON_SIMPLEX01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 50

  string output_filename = "uniform_on_simplex01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST25\n";
  cout << "  UNIFORM_ON_SIMPLEX01_MAP maps uniform \n";
  cout << "  points onto the unit simplex.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_on_simplex01_map ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test26 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST26 tests UNIFORM_ON_SPHERE01_MAP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 50

  string output_filename = "uniform_on_sphere01_map.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST26\n";
  cout << "  UNIFORM_ON_SPHERE01_MAP maps uniform\n";
  cout << "  points onto the unit sphere.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_on_sphere01_map ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test264 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST264 tests UNIFORM_ON_SPHERE01_PATCH_TP.
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
{
# define N 50
# define PI 3.141592653589793

  string output_filename = "uniform_on_sphere01_patch_tp.txt";
  double phi1;
  double phi2;
  int seed = 123456789;
  double theta1;
  double theta2;
  double *tp;

  phi1 = 75.0 * ( PI / 180.0 );
  phi2 = 90.0 * ( PI / 180.0 );
  theta1 =  0.0 * ( PI / 360.0 );
  theta2 = 30.0 * ( PI / 360.0 );

  cout << "\n";
  cout << "TEST264\n";
  cout << "  UNIFORM_ON_SPHERE01_PATCH_TP maps uniform\n";
  cout << "  points onto a TP (THETA,PHI) patch of the unit sphere.\n";
  cout << "\n";
  cout << "  Spatial dimension =          " << 3 << "\n";
  cout << "  Data dimension =             " << 2 << "\n";
  cout << "  Number of points N =         " << N       << "\n";
  cout << "  Latitudinal angle PHI1 =     " << phi1    << "\n";
  cout << "  Latitudinal angle PHI2 =     " << phi2    << "\n";
  cout << "  Longitudinal angle THETA1 =  " << theta1  << "\n";
  cout << "  Longitudinal angle THETA2 =  " << theta2  << "\n";
  cout << "  Initial random number SEED = " << seed    << "\n";

  tp = uniform_on_sphere01_patch_tp ( N, phi1, phi2, theta1, theta2, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, 2, N, tp );

  delete [] tp;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef N
# undef PI
}
//****************************************************************************80

void test265 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST265 tests UNIFORM_ON_SPHERE01_PATCH_XYZ.
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
{
# define DIM_NUM 3
# define N 50
# define PI 3.141592653589793

  string output_filename = "uniform_on_sphere01_patch_xyz.txt";
  double phi1;
  double phi2;
  int seed = 123456789;
  double theta1;
  double theta2;
  double *x;

  phi1 = 75.0 * ( PI / 180.0 );
  phi2 = 90.0 * ( PI / 180.0 );
  theta1 =  0.0 * ( PI / 360.0 );
  theta2 = 30.0 * ( PI / 360.0 );

  cout << "\n";
  cout << "TEST265\n";
  cout << "  UNIFORM_ON_SPHERE01_PATCH_XYZ maps uniform\n";
  cout << "  points onto an XYZ patch of the unit sphere.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM << "\n";
  cout << "  Number of points N =         " << N       << "\n";
  cout << "  Latitudinal angle PHI1 =     " << phi1    << "\n";
  cout << "  Latitudinal angle PHI2 =     " << phi2    << "\n";
  cout << "  Longitudinal angle THETA1 =  " << theta1  << "\n";
  cout << "  Longitudinal angle THETA2 =  " << theta2  << "\n";
  cout << "  Initial random number SEED = " << seed    << "\n";

  x = uniform_on_sphere01_patch_xyz ( N, phi1, phi2, theta1, theta2, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
# undef PI
}
//****************************************************************************80

void test267 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST267 tests UNIFORM_ON_SPHERE01_TRIANGLE_XYZ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2005
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 3
# define N 500

  string output_filename = "uniform_on_sphere01_triangle_xyz.txt";
  int seed = 123456789;
  double *v1;
  double *v2;
  double *v3;
  double *x;

  cout << "\n";
  cout << "TEST267\n";
  cout << "  UNIFORM_ON_SPHERE01_TRIANGLE_XYZ maps uniform\n";
  cout << "  points onto a spherical triangle using XYZ coordinates.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << 3       << "\n";
  cout << "  Number of points N =         " << N       << "\n";
  cout << "  Initial random number SEED = " << seed    << "\n";

  if ( true )
  {
    v1 = uniform_on_sphere01_map ( 3, 1, &seed );
    v2 = uniform_on_sphere01_map ( 3, 1, &seed );
    v3 = uniform_on_sphere01_map ( 3, 1, &seed );
  }
  else
  {
    v1 = r8vec_zero_new ( 3 );
    v1[0] = 1.0;
    v2 = r8vec_zero_new ( 3 );
    v2[1] = 1.0;
    v3 = r8vec_zero_new ( 3 );
    v3[2] = 1.0;
  }

  cout << "\n";
  cout << "  Vertices of spherical triangle:\n";
  cout << "\n";
  cout << "  V1: (" << v1[0] << ", " << v1[1] << ", " << v1[2] << ")\n";
  cout << "  V2: (" << v2[0] << ", " << v2[1] << ", " << v2[2] << ")\n";
  cout << "  V3: (" << v3[0] << ", " << v3[1] << ", " << v3[2] << ")\n";

  x = uniform_on_sphere01_triangle_xyz ( N, v1, v2, v3, &seed );

  cout << "\n";
  cout << "  Final random number SEED =   " << seed   << "\n";

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  delete [] v1;
  delete [] v2;
  delete [] v3;

  return;
# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test27 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST27 tests UNIFORM_WALK
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2008
//
//  Author:
//
//    John Burkardt
//
{
# define DIM_NUM 2
# define N 400

  string output_filename = "uniform_walk.txt";
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST27:\n";
  cout << "  UNIFORM_WALK generates uniform random walk points.\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM =  " << DIM_NUM      << "\n";
  cout << "  Number of points N =         " << N      << "\n";
  cout << "  Initial random number SEED = " << seed   << "\n";

  x = uniform_walk ( DIM_NUM, N, &seed );

  cout << "  Final random number SEED =   " << seed   << "\n";

  scale_to_block01 ( DIM_NUM, N, x );

  r8mat_write ( output_filename, DIM_NUM, N, x );

  delete [] x;

  cout << "\n";
  cout << "  Data written to file \"" << output_filename << "\".\n";

  return;
# undef DIM_NUM
# undef N
}
