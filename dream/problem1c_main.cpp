# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "dream_user.hpp"
# include "pdflib.hpp"
# include "rnglib.hpp"
# include "problem1c_covariance.hpp"

extern Covariance c;

int main ( );
void test01 ( int par_num, int sample_num );
double r8_huge ( );
double *r8mat_covariance ( int m, int n, double x[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PROBLEM1C_MAIN.
//
//  Discussion:
//
//    The coding of PROBLEM1 is tricky enough that I want to be able to
//    try it out independently of the DREAM code.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int chain_num;
  int cr_num;
  int gen_num;
  int pair_num;
  int par_num;
  int sample_num;

  timestamp ( );
  cout << "\n";
  cout << "PROBLEM1C_MAIN:\n";
  cout << "  C++ version\n";
  cout << "  This version uses a CLASS to store the covariance.\n";
//
//  Initialize the random number generator library.
//
  initialize ( );
//
//  By calling PROBLEM_SIZE, we implicitly set up the covariance as well.
//
  problem_size ( chain_num, cr_num, gen_num, pair_num, par_num );

  sample_num = 10000;
  test01 ( par_num, sample_num );
//
//  Terminate.
//
  cout << "\n";
  cout << "PROBLEM1C_MAIN:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( int par_num, int sample_num )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 calls the sampling function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *array;
  double *cov_sample;
  int i;
  int j;
  double z;
  double *zp;
  double *zp_array;
  double *zp_ave;
  double *zp_max;
  double *zp_mean;
  double *zp_min;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Call PRIOR_SAMPLE many times.\n";
  cout << "  Compare statistics to PDF parameters.\n";
  cout << "  Note that the covariance estimate can be very bad\n";
  cout << "  unless the matrix is strongly diagonal.\n";
  cout << "\n";
  cout << "  Parameter dimension is " << par_num << "\n";
  cout << "  Number of samples is " << sample_num << "\n";
//
//  Compute N multinormal samples.
//
  zp_array = ( double * ) malloc ( par_num * sample_num * sizeof ( double ) );
  for ( j = 0; j < sample_num; j++ )
  {
    zp = prior_sample ( par_num );
    for ( i = 0; i < par_num; i++ )
    {
      zp_array[i+j*par_num] = zp[i];
    }
    delete [] zp;
  }

  zp_ave = ( double * ) malloc ( par_num * sizeof ( double ) );
  zp_max = ( double * ) malloc ( par_num * sizeof ( double ) );
  zp_min = ( double * ) malloc ( par_num * sizeof ( double ) );

  for ( i = 0; i < par_num; i++ )
  {
    zp_min[i] =   r8_huge ( );
    zp_max[i] = - r8_huge ( );
    zp_ave[i] = 0.0;
    for ( j = 0; j < sample_num; j++ )
    {
      z = zp_array[i+j*par_num];
      if ( z < zp_min[i] )
      {
        zp_min[i] = z;
      }
      if ( zp_max[i] < z )
      {
        zp_max[i] = z;
      }
      zp_ave[i] = zp_ave[i] + z;
    }
    zp_ave[i] = zp_ave[i] / ( double ) ( sample_num );
  }

  zp_mean = c.mean_get ( );

  cout << "\n";
  cout << " Index       Min            Ave              Max             MU\n";
  cout << "\n";
  for ( i = 0; i < par_num; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << zp_min[i]
         << "  " << setw(14) << zp_ave[i]
         << "  " << setw(14) << zp_max[i]
         << "  " << setw(14) << zp_mean[i] << "\n";
  }

  cov_sample = r8mat_covariance ( par_num, sample_num, zp_array );

  r8mat_print ( par_num, par_num, cov_sample, "  Sample covariance:" );

  array = c.array_get ( );
  r8mat_print ( par_num, par_num, array, "  PDF covariance:" );

  delete [] array;
  delete [] cov_sample;
  delete [] zp_array;
  delete [] zp_ave;
  delete [] zp_max;
  delete [] zp_mean;
  delete [] zp_min;

  return;
}
//****************************************************************************80

double r8_huge ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_HUGE returns a "huge" R8.
//
//  Discussion:
//
//    The value returned by this function is NOT required to be the
//    maximum representable R8.  This value varies from machine to machine,
//    from compiler to compiler, and may cause problems when being printed.
//    We simply want a "very large" but non-infinite number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_HUGE, a "huge" R8 value.
//
{
  double value;

  value = 1.0E+30;

  return value;
}
//****************************************************************************80

double *r8mat_covariance ( int m, int n, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_COVARIANCE computes the sample covariance of a set of vector data.
//
//  Discussion:
//
//    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 June 2013
//
//  Author:
//
//    John Burkardt.
//
//  Parameters:
//
//    Input, int M, the size of a single data vectors.
//
//    Input, int N, the number of data vectors.
//    N should be greater than 1.
//
//    Input, double X[M*N], an array of N data vectors, each
//    of length M.
//
//    Output, double C[M*M], the covariance matrix for the data.
//
{
  double *c;
  int i;
  int j;
  int k;
  double *x_mean;

  c = new double[m*m];
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = 0.0;
    }
  }
//
//  Special case of N = 1.
//
  if ( n == 1 )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+i*m] = 1.0;
    }
    return c;
  }
//
//  Determine the sample means.
//
  x_mean = new double[m];
  for ( i = 0; i < m; i++ )
  {
    x_mean[i] = 0.0;
    for ( j = 0; j < n; j++ )
    {
      x_mean[i] = x_mean[i] + x[i+j*m];
    }
    x_mean[i] = x_mean[i] / ( double ) ( n );
  }
//
//  Determine the sample covariance.
//
  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      for ( k = 0; k < n; k++ )
      {
        c[i+j*m] = c[i+j*m] 
          + ( x[i+k*m] - x_mean[i] ) * ( x[j+k*m] - x_mean[j] );
      }
    }
  }

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      c[i+j*m] = c[i+j*m] / ( double ) ( n - 1 );
    }
  }

  delete [] x_mean;

  return c;
}
