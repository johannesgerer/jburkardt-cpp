# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "asa053.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA053_PRB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA053_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA053 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA053_PRB:\n";
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
//    TEST01 generates a random Wishart variate.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define NP 3

  double d[(NP*(NP+1))/2] = {
    3.0, 
    2.0, 4.0, 
    1.0, 2.0, 5.0 };
  int i;
  int n;
  int np = NP;
  double *sa;
  int seed;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Generate a single Wishart deviate.\n";

  n = 1;
  seed = 123456789;

  cout << "\n";
  cout << "  The number of variables is " << np << "\n";
  cout << "  The number of degrees of freedom is " << n << "\n";

  r8utp_print ( np, d, "  The upper Cholesky factor:" );

  sa = wshrt ( d, n, np, seed );

  r8pp_print ( np, sa, "  The sample matrix:" );

  delete [] sa;

  return;
# undef NP
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 averages many Wishart samples.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define NP 3

  double d[(NP*(NP+1))/2] = {
    3.0, 
    2.0, 4.0, 
    1.0, 2.0, 5.0 };
  int i;
  int j;
  int k;
  int ki;
  int kj;
  int n;
  int np = NP;
  int npp;
  double *sa;
  double *s_average;
  double *sigma;
  int seed;
  int test_num = 100000;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Generate many Wishart deviates.\n";
  cout << "  Compare to D' * D * np / n\n";
  n = 2;
  npp = ( np * ( np + 1 ) ) / 2;
  seed = 123456789;

  cout << "\n";
  cout << "  The number of variables is " << np << "\n";
  cout << "  The number of degrees of freedom is " << n << "\n";

  r8utp_print ( np, d, "  The upper Cholesky factor:" );

  s_average = new double[npp];

  for ( j = 0; j < npp; j++ )
  {
    s_average[j] = 0.0;
  }

  for ( i = 1; i <= test_num; i++ )
  {
    sa = wshrt ( d, n, np, seed );
    for ( j = 0; j < npp; j++ )
    {
      s_average[j] = s_average[j] + sa[j];
    }
    free ( sa );
  }

  for ( j = 0; j < npp; j++ )
  {
    s_average[j] = s_average[j] / ( double ) ( test_num );
  }

  r8pp_print ( np, s_average, "  The averaged matrix:" );
//
//  Compare the result to ( D' * D ) * np / n.
//
  sigma = new double[np*np];

  for ( i = 0; i < np; i++ )
  {
    for ( j = 0; j < np; j++ )
    {
      for ( k = 0; k <= i4_min ( i, j ); k++ )
      {
        ki = k + ( i * ( i + 1 ) ) / 2;
        kj = k + ( j * ( j + 1 ) ) / 2;
        sigma[i+j*np] = sigma[i+j*np] + d[ki] * d[kj];
      }
      sigma[i+j*np] = sigma[i+j*np] * ( double ) np / ( double ) n;
    }
  }

  r8mat_print ( np, np, sigma, "  Expected result:" );

  delete [] s_average;
  delete [] sigma;

  return;
# undef NP
}
