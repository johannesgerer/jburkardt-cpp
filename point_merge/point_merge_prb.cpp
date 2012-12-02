# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "point_merge.hpp"

int main ( );
void test01 ( int m, int n, int n_unique, int seed );
void test02 ( int m, int n, int n_unique, double tol, int seed );
void test03 ( int m, int n, int n_unique, double tol, int seed );
void test04 ( int m, int n, int n_unique, double tol, int seed );
void test05 ( int m, int n, int n_unique, double tol, int seed );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POINT_MERGE_PRB.
//
//  Discussion:
//
//    Compare correctness of the codes.
//
//    Compare speed of the codes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  int n_unique;
  int seed;
  double tol;

  timestamp ( );

  cout << " \n";
  cout << "POINT_MERGE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the POINT_MERGE library.\n";
//
//  TEST01 gives me some confidence that, at least for zero-tolerance,
//  the radial approach is accurate, as compared to the "Unique count"
//  (which cannot be extended to a tolerance version in multidimensions)
//  and the "Tol Unique Count", which is an O(N^2) algorithm.
//
  m = 3;
  n = 10;
  n_unique = 7;
  seed = 123456789;
  test01 ( m, n, n_unique, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  seed = 987654321;
  test01 ( m, n, n_unique, seed );
//
//  In TEST02, I want to compute the same data, but with "blurred"
//  duplicates, and a tolerance version of the radial approach,
//  compared to "Tol Unique Count".
//
  m = 3;
  n = 10;
  n_unique = 7;
  tol = 0.00001;
  seed = 123456789;
  test02 ( m, n, n_unique, tol, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  tol = 0.00001;
  seed = 987654321;
  test02 ( m, n, n_unique, tol, seed );
//
//  In TEST03, I want to measure the time required for a sequence
//  of increasingly hard problems.
//
  m = 3;
  n = 100;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 1000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 10000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test03 ( m, n, n_unique, tol, seed );

  if ( false )
  {
    m = 3;
    n = 100000;
    n_unique = n / 2;
    tol = 0.00001;
    seed = 123456789;
    test03 ( m, n, n_unique, tol, seed );
  }
//
//  In TEST04, repeat TEST02, but now compute the index vector.
//
  m = 3;
  n = 10;
  n_unique = 7;
  tol = 0.00001;
  seed = 123456789;
  test04 ( m, n, n_unique, tol, seed );

  m = 4;
  n = 20;
  n_unique = 11;
  tol = 0.00001;
  seed = 987654321;
  test04 ( m, n, n_unique, tol, seed );
//
//  In TEST05, I want to measure the time required for a sequence
//  of increasingly hard problems.
//
  m = 3;
  n = 100;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 1000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  m = 3;
  n = 10000;
  n_unique = n / 2;
  tol = 0.00001;
  seed = 123456789;
  test05 ( m, n, n_unique, tol, seed );

  if ( false )
  {
    m = 3;
    n = 100000;
    n_unique = n / 2;
    tol = 0.00001;
    seed = 123456789;
    test05 ( m, n, n_unique, tol, seed );
  }
//
//  Terminate.
//
  cout << " \n";
  cout << "POINT_MERGE_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << " \n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int m, int n, int n_unique, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests uniqueness counting with no tolerance.
//
//  Discussion:
//
//    POINT_UNIQUE_COUNT uses an O(N) algorithm.
//    POINT_RADIAL_UNIQUE_COUNT uses an algorithm that should be,
//      in general, O(N);
//    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
//
//    For this test, we just want to make sure the algorithms agree
//    in the counting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double tol;
  int unique_num;

  cout << " \n";
  cout << "TEST01\n";
  cout << "  To count the unique columns in an R8COL, we call\n";
  cout << "  POINT_UNIQUE_COUNT,\n";
  cout << "  POINT_RADIAL_UNIQUE_COUNT, (with random center)\n";
  cout << "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n";
  cout << " \n";
  cout << "  M =     " << m << "\n";
  cout << "  N =     " << n << "\n";
  cout << "  SEED =  " << seed << "\n";

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );

  cout << " \n";
  cout << "  N_UNIQUE =                  " << n_unique << "\n";

  unique_num = point_unique_count ( m, n, a );
  cout << "  POINT_UNIQUE_COUNT =        " << unique_num << "\n";

  unique_num = point_radial_unique_count ( m, n, a, &seed );
  cout << "  POINT_RADIAL_UNIQUE_COUNT = " << unique_num << "\n";

  tol = 0.0;
  unique_num = point_tol_unique_count ( m, n, a, tol );
  cout << "  POINT_TOL_UNIQUE_COUNT =    " << unique_num << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test02 ( int m, int n, int n_unique, double tol, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests uniqueness counting with a tolerance.
//
//  Discussion:
//
//    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
//      in general, O(N);
//    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
//
//    For this test, we just want to make sure the algorithms agree
//    in the counting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  int i;
  int j;
  double *r;
  double r_norm;
  int unique_num;

  cout << " \n";
  cout << "TEST02\n";
  cout << "  To count the unique columns in an R8COL, we call\n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n";
  cout << "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n";
  cout << " \n";
  cout << "  M =     " << m << "\n";
  cout << "  N =     " << n << "\n";
  cout << "  TOL =  " << tol << "\n";
  cout << "  SEED =  " << seed << "\n";

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );
//
//  The form of the tolerance test means that if two vectors are initially
//  equal, they remain "tolerably equal" after the addition of random
//  perturbation vectors whose 2-norm is no greater than TOL/2.
//
  r = new double[m];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  delete [] r;

  r8mat_transpose_print ( m, n, a, "  Blurred matrix:" );

  cout << " \n";
  cout << "  N_UNIQUE =                      " << n_unique << "\n";

  unique_num = point_radial_tol_unique_count ( m, n, a, tol, &seed );
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT = " << unique_num << "\n";

  unique_num = point_tol_unique_count ( m, n, a, tol );
  cout << "  POINT_TOL_UNIQUE_COUNT =        " << unique_num << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test03 ( int m, int n, int n_unique, double tol, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 compares timings for two uniqueness counters.
//
//  Discussion:
//
//    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
//      in general, O(N);
//    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double ctime;
  int i;
  int j;
  double *r;
  double r_norm;
  int unique_num;

  cout << " \n";
  cout << "TEST03\n";
  cout << "  To count the unique columns in an R8COL, we call\n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n";
  cout << "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n";
  cout << " \n";
  cout << "  M =     " << m << "\n";
  cout << "  N =     " << n << "\n";
  cout << "  TOL =  " << tol << "\n";
  cout << "  SEED =  " << seed << "\n";

  a = r8col_duplicates ( m, n, n_unique, &seed );
//
//  The form of the tolerance test means that if two vectors are initially
//  equal, they remain "tolerably equal" after the addition of random
//  perturbation vectors whose 2-norm is no greater than TOL/2.
//
  r = new double[m];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  delete [] r;

  cout << " \n";
  cout << "  N_UNIQUE =                      " << n_unique << "\n";

  ctime = cpu_time ( );
  unique_num = point_radial_tol_unique_count ( m, n, a, tol, &seed );
  ctime = cpu_time ( ) - ctime;
  cout << " \n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT = " << unique_num << "\n";
  cout << "  CPU_TIME = " << ctime << "\n";

  ctime = cpu_time ( );
  unique_num = point_tol_unique_count ( m, n, a, tol );
  ctime = cpu_time ( ) - ctime;
  cout << " \n";
  cout << "  POINT_TOL_UNIQUE_COUNT =        " << unique_num << "\n";
  cout << "  CPU_TIME = " << ctime << "\n";

  delete [] a;

  return;
}
//****************************************************************************80

void test04 ( int m, int n, int n_unique, double tol, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests uniqueness indexing with a tolerance.
//
//  Discussion:
//
//    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
//      in general, O(N);
//    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
//
//    For this test, we just want to make sure the algorithms agree
//    in the counting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double dist;
  int i;
  int j;
  int k;
  double *r;
  double r_norm;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << " \n";
  cout << "TEST04\n";
  cout << "  To index the unique columns in an R8COL, we call\n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n";
  cout << "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n";
  cout << " \n";
  cout << "  M =     " << m << "\n";
  cout << "  N =     " << n << "\n";
  cout << "  TOL =  " << tol << "\n";
  cout << "  SEED =  " << seed << "\n";

  a = r8col_duplicates ( m, n, n_unique, &seed );

  r8mat_transpose_print ( m, n, a, "  Matrix with N_UNIQUE unique columns:" );
//
//  The form of the tolerance test means that if two vectors are initially
//  equal, they remain "tolerably equal" after the addition of random
//  perturbation vectors whose 2-norm is no greater than TOL/2.
//
  r = new double[m];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  delete [] r;

  r8mat_transpose_print ( m, n, a, "  Blurred matrix:" );

  cout << " \n";
  cout << "  N_UNIQUE =                      " << n_unique << "\n";

  undx = new int[n];
  xdnu = new int[n];

  unique_num = point_radial_tol_unique_index ( m, n, a, tol, &seed, undx,
    xdnu );

  cout << " \n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_INDEX\n";
  cout << "  Unique_num = " << unique_num << "\n";

  i4vec_print ( unique_num, undx, "  UNDX:" );

  i4vec_print ( n, xdnu, "  XDNU:" );

  cout << " \n";
  cout << "  List of nonunique points P(J), represented by\n";
  cout << "  point with index I(J).\n";
  cout << " \n";
  cout << "  J, P(J)\n";
  cout << "  I(J), P(I(J))\n";
  cout << "  || P(J) - P(I(J)) || (should be <= TOL)\n";
  cout << " \n";
  for ( j = 0; j < n; j++ )
  {
    k = undx[xdnu[j]];
    if ( j != k )
    {
      cout << " \n";
      cout << "  " << setw(4) << j;
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << a[i+j*m];
      }
      cout << "\n";
      cout << "  " << setw(4) << k;
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << a[i+k*m];
      }
      cout << "\n";
      dist = 0.0;
      for ( i = 0; i < m; i++ )
      {
        dist = dist + pow ( a[i+j*m] - a[i+k*m], 2 );
      }
      dist = sqrt ( dist );
      cout << "          " << setw(10) <<  dist << "\n";
    }
  }
//
//  The interpretation of XDNU is simpler for POINT_TOL_UNIQUE_INDEX.
//
  unique_num = point_tol_unique_index ( m, n, a, tol, xdnu );

  cout << " \n";
  cout << "  POINT_TOL_UNIQUE_INDEX\n";
  cout << "  Unique_num = " << unique_num << "\n";
  cout << " \n";
  cout << "  List of nonunique points P(J), represented by\n";
  cout << "  point with index I(J).\n";
  cout << " \n";
  cout << "  J, P(J)\n";
  cout << "  I(J), P(I(J))\n";
  cout << "  || P(J) - P(I(J)) || (should be <= TOL)\n";
  cout << " \n";
  for ( j = 0; j < n; j++ )
  {
    k = xdnu[j];
    if ( j != k )
    {
      cout << " \n";
      cout << "  " << setw(4) << j;
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << a[i+j*m];
      }
      cout << "\n";
      cout << "  " << setw(4) << k;
      for ( i = 0; i < m; i++ )
      {
        cout << "  " << a[i+k*m];
      }
      cout << "\n";
      dist = 0.0;
      for ( i = 0; i < m; i++ )
      {
        dist = dist + pow ( a[i+j*m] - a[i+k*m], 2 );
      }
      dist = sqrt ( dist );
      cout << "          " << setw(10) <<  dist << "\n";
    }
  }

  delete [] a;
  delete [] undx;
  delete [] xdnu;

  return;
}
//****************************************************************************80

void test05 ( int m, int n, int n_unique, double tol, int seed )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 times uniqueness indexing with a tolerance.
//
//  Discussion:
//
//    POINT_RADIAL_TOL_UNIQUE_COUNT uses an algorithm that should be,
//      in general, O(N);
//    POINT_TOL_UNIQUE_COUNT uses an O(N^2) algorithm.
//
//    For this test, we just want to make sure the algorithms agree
//    in the counting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 July 2010
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double ctime;
  double dist;
  int i;
  int j;
  double *r;
  double r_norm;
  int *undx;
  int unique_num;
  int *xdnu;

  cout << " \n";
  cout << "TEST05\n";
  cout << "  We time the computations in TEST04, calling\n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_COUNT, (with random center)\n";
  cout << "  POINT_TOL_UNIQUE_COUNT, (with zero tolerance)\n";
  cout << " \n";
  cout << "  M =     " << m << "\n";
  cout << "  N =     " << n << "\n";
  cout << "  TOL =  " << tol << "\n";
  cout << "  SEED =  " << seed << "\n";

  a = r8col_duplicates ( m, n, n_unique, &seed );
//
//  The form of the tolerance test means that if two vectors are initially
//  equal, they remain "tolerably equal" after the addition of random
//  perturbation vectors whose 2-norm is no greater than TOL/2.
//
  r = new double[m];

  for ( j = 0; j < n; j++ )
  {
    r8vec_uniform_01 ( m, &seed, r );
    r_norm = r8vec_norm_l2 ( m, r );
    for ( i = 0; i < m; i++ )
    {
      a[i+j*m] = a[i+j*m] + 0.5 * tol * r[i] / r_norm;
    }
  }

  delete [] r;

  cout << " \n";
  cout << "  N_UNIQUE =                      " << n_unique << "\n";

  undx = new int[n];
  xdnu = new int[n];

  ctime = cpu_time ( );
  unique_num = point_radial_tol_unique_index ( m, n, a, tol, &seed, undx,
    xdnu );
  ctime = cpu_time ( ) - ctime;

  cout << " \n";
  cout << "  POINT_RADIAL_TOL_UNIQUE_INDEX\n";
  cout << "  Unique_num = " << unique_num << "\n";
  cout << "  Time = " << ctime << "\n";

  ctime = cpu_time ( );
  unique_num = point_tol_unique_index ( m, n, a, tol, xdnu );
  ctime = cpu_time ( ) - ctime;

  cout << " \n";
  cout << "  POINT_TOL_UNIQUE_INDEX\n";
  cout << "  Unique_num = " << unique_num << "\n";
  cout << "  Time = " << ctime << "\n";

  delete [] a;
  delete [] undx;
  delete [] xdnu;

  return;
}
