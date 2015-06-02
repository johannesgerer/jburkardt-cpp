# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <cstring>

using namespace std;

# include "cvt.hpp"

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

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CVT_PRB.
//
//  Discussion:
//
//    CVT_PRB tests the CVT library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CVT_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CVT library.\n";

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
//
//  Terminate.
//
  cout << "\n";
  cout << "CVT_PRB\n";
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
//    TEST01 tests CVT with uniform initialization and uniform sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST02 repeats test 1, but uses twice as many iterations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with twice the number of iterations.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 80;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
   
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

  return;

# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 repeats test 1 but uses 100 times as many sample points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with 100 times the sample points.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 1000000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
   
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

  return;

# undef DIM_NUM
# undef N
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 repeats test 1 with uniform initialization and Halton sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with Halton sampling.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 1;
  sample_num = 10000;
  strcpy ( sample_string, "halton" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST05 repeats test 1 with uniform initialization and grid sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST05\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with grid sampling.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 2;
  sample_num = 10000;
  strcpy ( sample_string, "grid" );
  seed = 123456789;

  seed_init = seed;
  
  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
   
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST06 repeats test 1 with uniform initialization and C++ RANDOM sampling.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with C++ RANDOM sampling.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = -1;
  sample_num = 10000;
  strcpy ( sample_string, "random" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
   
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST07 tests CVT with a different seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1 with a different seed.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 987654321;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST08 repeats test 1 with a different batch size.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST08\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1 with a different batch size.\n";

  batch = 5;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST09 repeats test 1, but with IT_FIXED = IT_MAX.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Repeat test 1, but with IT_FIXED = IT_MAX.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = it_max;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST10 generates 100 points in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define DIM_NUM 3

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  Compute 100 points in 3D.\n";

  batch = 1000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,  
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print_some ( DIM_NUM, N, r, 1, 1, DIM_NUM, 10, 
    "  First 10 Generators (rows):" );

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test11 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST11 tests CVT.
//
//  Discussion:
//
//    In this test, we initialize the generators to grid points; this is 
//    an unstable CVT solution.  The data would "prefer" to be in a
//    different form.  However, even if we take 2000 steps of CVT iteration,
//    the data is still only slowly progressing towards that other 
//    configuration.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 16
# define DIM_NUM 2

  int batch;
  double energy;
  int i;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  int j;
  int ngrid;
  double r[DIM_NUM*N];
  int rank;
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;
  int tuple[DIM_NUM];

  cout << "\n";
  cout << "TEST11\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "\n";
  cout << "  In this test, we initialize the generators to\n";
  cout << "  grid points; this is an unstable CVT solution.\n";

  batch = 1000;
  init = 4;
  strcpy ( init_string, "user initialization" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 1000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;
//
//  Initialize the tuple generator.
//
  rank = -1;
  ngrid = 4;
  tuple_next_fast ( ngrid, DIM_NUM, rank, tuple );
//
//  Pick points on a grid.
//
  for ( j = 0; j < N; j++ )
  {
    rank = j;
    tuple_next_fast ( ngrid, DIM_NUM, rank, tuple );
    for ( i = 0; i < DIM_NUM; i++ )
    {
      r[i+j*DIM_NUM] = ( double ) ( 2 * tuple[i] - 1 ) 
                     / ( double ) ( 2 * ngrid );
    }
  }
  r8mat_transpose_print ( DIM_NUM, N, r, "  Initial generators (rows):" );

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed,
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

  return;
# undef N
# undef DIM_NUM
}
//****************************************************************************80

void test12 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST12 tests CVT with 'RANDOM' initialization.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 2

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST12\n";
  cout << "  The \"random\" initialization option calls the\n";
  cout << "  system random number generator.  There is some\n";
  cout << "  question about whether this works correctly.\n";
  cout << "\n";
  cout << "  The test is as follows:\n";
  cout << "\n";
  cout << "  CVT call #1:\n";
  cout << "\n";
  cout << "    DIM_NUM      =      2\n";
  cout << "    N         =     10\n";
  cout << "    INIT      =     -1\n";
  cout << "    IT_MAX    =      0\n";
  cout << "    SEED      = 100000\n";
  cout << "\n";
  cout << "    Print output values of SEED and R #1.\n";
  cout << "\n";
  cout << "  CVT call #2: (jump SEED)\n";
  cout << "\n";
  cout << "    DIM_NUM      =      2\n";
  cout << "    N         =     10\n";
  cout << "    INIT      =     -1\n";
  cout << "    IT_MAX    =      0\n";
  cout << "    SEED      = 200000.\n";
  cout << "\n";
  cout << "    Print output values of SEED and R #2.\n";
  cout << "\n";
  cout << "  CVT call #3: (restore SEED)\n";
  cout << "\n";
  cout << "    DIM_NUM      =      2\n";
  cout << "    N         =     10\n";
  cout << "    INIT      =     -1\n";
  cout << "    IT_MAX    =      0\n";
  cout << "    SEED_INIT = 100000\n";
  cout << "\n";
  cout << "    Print output values of SEED and R #3.\n";
  cout << "\n";
  cout << "  We expect that:\n";
  cout << "  * the values of R #1 and R #2 differ;\n";
  cout << "  AND\n";
  cout << "  * the values of R #1 and R #3 agree.\n";
//
//  Run #1.
//
  batch = 1000;
  init = -1;
  strcpy ( init_string, "random" );
  it_max = 0;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 100000;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );
//
//  Run #2.
//
  batch = 1000;
  init = -1;
  strcpy ( init_string, "random" );
  it_max = 0;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 200000;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );
//
//  Run #3.
//
  batch = 1000;
  init = -1;
  strcpy ( init_string, "random" );
  it_max = 0;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 100000;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

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
//    TEST13 tests CVT with the "user" routine.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 100
# define DIM_NUM 2

  int batch;
  double energy;
  string file_out_name = "cvt_circle.txt";
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST13\n";
  cout << "  CVT computes a Centroidal Voronoi Tessellation.\n";
  cout << "  In this example, we call the \"USER\" routine,\n";
  cout << "  which allows the user to define the geometry and\n";
  cout << "  density implicitly, by returning sample points.\n";

  batch = 1000;
  init = 3;
  strcpy ( init_string, "user" );
  it_max = 40;
  it_fixed = 1;
  sample = 3;
  sample_num = 10000;
  strcpy ( sample_string, "user" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_write ( file_out_name, DIM_NUM, N, r );

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
//    TEST14 generates a CVT in the interval [0,1] using 10 points..
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 September 2010
//
//  Author:
//
//    John Burkardt
//
{
# define N 10
# define DIM_NUM 1

  int batch;
  double energy;
  int init;
  char init_string[80];
  double it_diff;
  int it_fixed;
  int it_max;
  int it_num;
  double r[DIM_NUM*N];
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;

  cout << "\n";
  cout << "TEST14\n";
  cout << "  Generate a CVT in the interval [0,1] using 10 points.\n";

  batch = 10000;
  init = 0;
  strcpy ( init_string, "uniform" );
  it_max = 40;
  it_fixed = 1;
  sample = 0;
  sample_num = 10000;
  strcpy ( sample_string, "uniform" );
  seed = 123456789;

  seed_init = seed;

  cvt ( DIM_NUM, N, batch, init, sample, sample_num, it_max, it_fixed, 
    &seed, r, &it_num, &it_diff, &energy );

  cout << "\n";
  cout << "  Dimension DIM_NUM =        "  << DIM_NUM       << "\n";
  cout << "  Number of points N =       "  << N             << "\n";
  cout << "  Initial SEED =             "  << seed_init     << "\n";
  cout << "  Current SEED =             "  << seed          << "\n";
  cout << "  INIT =                    \"" << init_string   << "\".\n";
  cout << "  Max iterations IT_MAX =    "  << it_max        << "\n";
  cout << "  IT_FIXED (fixed samples) = "  << it_fixed      << "\n";
  cout << "  Iterations IT_NUM =        "  << it_num        << "\n";
  cout << "  Difference IT_DIFF =       "  << it_diff       << "\n";
  cout << "  CVT ENERGY =               "  << energy        << "\n";
  cout << "  SAMPLE =                  \"" << sample_string << "\".\n";
  cout << "  Samples SAMPLE_NUM    =    "  << sample_num    << "\n";
  cout << "  Sampling BATCH size =      "  << batch         << "\n";
  cout << "  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";
  
  r8mat_transpose_print ( DIM_NUM, N, r, "  Generators (rows):" );

  return;

# undef DIM_NUM
# undef N
}
