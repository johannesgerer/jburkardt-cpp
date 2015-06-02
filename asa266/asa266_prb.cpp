# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "asa266.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test085 ( );
void test09 ( );
void test10 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ASA266_PRB.
//
//  Discussion:
//
//    ASA266_PRB tests the ASA266 library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "ASA266_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the ASA266 library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test085 ( );
  test09 ( );
  test10 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "ASA266_PRB:\n";
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
//    TEST01 tests ALNORM, NORMP, NPROB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double ccdf1;
  double ccdf2;
  double ccdf3;
  double cdf1;
  double cdf2;
  double cdf3;
  int i;
  int ntest = 16;
  double pdf2;
  double pdf3;
  int upper;
  double x;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  ALNORM,\n";
  cout << "  NORMP, and\n";
  cout << "  NPROB are routines that compute the cumulative\n";
  cout << "  density function for the normal distribution.\n";
  cout << "\n";
  cout << "  X  CDF1  1-CDF1\n";
  cout << "     CDF2  1-CDF2  PDF2\n";
  cout << "     CDF3  1-CDF3  PDF3\n";

  for ( i = 0; i < ntest; i++ )
  {
    x = 3.0 * ( double ) ( i ) / ( double ) ( ntest - 1 );

    upper = 0;
    cdf1 = alnorm ( x, upper );

    upper = 1;
    ccdf1 = alnorm ( x, upper );

    normp ( x, &cdf2, &ccdf2, &pdf2 );

    nprob ( x, &cdf3, &ccdf3, &pdf3 );

    cout << "\n";
    cout << setw(14) << x
         << setw(14) << cdf1
         << setw(14) << ccdf1 << "\n";
    cout << "              "
         << setw(14) << cdf2
         << setw(14) << ccdf2
         << setw(14) << pdf2 << "\n";
    cout << "              "
         << setw(14) << cdf3
         << setw(14) << ccdf3
         << setw(14) << pdf3 << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests PPND, PPND16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  int ifault;
  int ntest = 9;
  double x1;
  double x2;

  ifault = 0;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  PPND,\n";
  cout << "  PPND16 compute the percentage \n";
  cout << "  points of the normal distribution.\n";
  cout << "\n";
  cout << "           CDF     PPND(CDF)   PPND16(CDF)\n";
  cout << "\n";

  for ( i = 1; i <= ntest; i++ )
  {
    cdf = ( double ) ( i ) / ( double ) ( ntest + 1 );
    x1 = ppnd ( cdf, &ifault );
    x2 = ppnd16 ( cdf, &ifault );
    cout << setw(14) << cdf
         << setw(14) << x1
         << setw(14) << x2 << "\n";
  }
  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 tests DIGAMMA, R8_PSI.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ntest = 10;
  double x;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  digamma(X) = d ( Log ( Gamma ( X ) ) ) / dX.\n";
  cout << "\n";
  cout << "  DIGAMMA and\n";
  cout << "  R8_PSI compute the digamma function:\n";
  cout << "\n";
  cout << "             X       DIGAMMA        R8_PSI\n";
  cout << "\n";

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    cout << setw(14) << x
         << setw(14) << digamma ( x )
         << setw(14) << r8_psi ( x ) << "\n";
  }
  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 tests TRIGAMMA.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ifault;
  int ntest = 10;
  double t;
  double x;

  cout << "\n";
  cout << "TEST04\n";
  cout << "  TRIGAMMA computes the trigamma function:\n";
  cout << "    trigamma(X) = d^2 ( Log ( Gamma ( X ) ) ) / dX^2.\n";
  cout << "\n";
  cout << "             X       TRIGAMMA\n";
  cout << "\n";

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    t = trigamma ( x, &ifault );
    cout << setw(14) << x
         << setw(14) << t << "\n";
  }
  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST05 tests ALNGAM, ALOGAM, R8_GAMMA_LOG, LNGAMMA;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ifault;
  double log1;
  double log2;
  double log3;
  double log4;
  int ntest = 10;
  double x;

  ifault = 0;
  
  cout << "\n";
  cout << "TEST05\n";
  cout << "  ALNGAM\n";
  cout << "  ALOGAM,\n";
  cout << "  R8_GAMMA_LOG, and\n";
  cout << "  LNGAMMA compute the logarithm of the gamma function.\n";
  cout << "\n";
  cout << "             X        ALNGAM        ALOGAM    R8_GAMMA_LOG     LNGAMMA\n";
  cout << "\n";

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    log1 = alngam ( x, &ifault );
    log2 = alogam ( x, &ifault );
    log3 = r8_gamma_log ( x );
    log4 = lngamma ( x, &ifault );

    cout << setw(14) << x
         << setw(14) << log1
         << setw(14) << log2
         << setw(14) << log3
         << setw(14) << log4 << "\n";
  }

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST06 tests GAMAIN, GAMMDS, GAMMAD.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double g1;
  double g2;
  double g3;
  int i;
  int ifault;
  int j;
  int ntest = 10;
  double p;
  double x;

  ifault = 0;

  cout << "\n";
  cout << "TEST06\n";
  cout << "  GAMAIN, \n";
  cout << "  GAMMDS and \n";
  cout << "  GAMMAD compute the incomplete Gamma integral.\n";
  cout << "\n";
  cout << "             X             P        GAMMDS        GAMMAD        GAMAIN\n";
  cout << "\n";

  for ( i = 1; i <= ntest; i++ )
  {
    x = ( double ) ( i ) / ( double ) ( ntest );

    cout << "\n";
    for ( j = 1; j <= ntest; j++ )
    {
      p = ( double ) ( j ) / ( double ) ( ntest );
      g1 = gammds ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g1 = -99.0;
      }

      g2 = gammad ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g2 = -99.0;
      }

      g3 = gamain ( x, p, &ifault );
      if ( ifault != 0 )
      {
        g3 = - 99.0;
      }
      cout << setw(14) << x
           << setw(14) << p
           << setw(14) << g1
           << setw(14) << g2
           << setw(14) << g3 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST07 tests PPCHI2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  double gg;
  int i;
  int ifault;
  int j;
  int nitest = 9;
  int njtest = 9;
  double v;
  double x1;

  ifault = 0;

  cout << "\n";
  cout << "TEST07\n";
  cout << "  PPCHI2 computes the percentage points\n";
  cout << "  of the chi squared distribution.\n";
  cout << "\n";
  cout << "      CDF      PPCHI2(CDF)\n";
  cout << "\n";

  for ( j = 1; j <= njtest; j++ )
  {
    v = ( double ) ( j );

    cout << "\n";
    cout << "  For Chi^2 parameter value = " << v << "\n";
    cout << "\n";

    for ( i = 1; i <= nitest; i++ )
    {
      cdf = ( double ) ( i ) / ( double ) ( nitest + 1 );
      gg = alngam ( v / 2.0, &ifault );
      x1 = ppchi2 ( cdf, v, gg, &ifault );
      cout << setw(14) << cdf
           << setw(14) << x1 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST08 tests DIRICHLET_ESTIMATE, DIRICHLET_MEAN, DIRICHLET_VARIANCE.
//
//  Discussion:
//
//    Canned data is used.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
# define ELEM_NUM 3
# define SAMPLE_NUM 23

  double *alpha;
  double alpha_sum;
  double aminus;
  double aplus;
  int elem_i;
  int elem_num = ELEM_NUM;
  double eps;
  double *g;
  int ifault;
  int init;
  double *mean;
  int niter;
  double rlogl;
  double s;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  double *v;
  double vari;
  double *variance;
  double x[SAMPLE_NUM*ELEM_NUM] = {
    0.178, 0.162, 0.083, 0.087, 0.078, 0.040, 0.049, 0.100, 0.075, 0.084,
    0.060, 0.089, 0.050, 0.073, 0.064, 0.085, 0.094, 0.014, 0.060, 0.031,
    0.025, 0.045, 0.0195,
    0.346, 0.307, 0.448, 0.474, 0.503, 0.456, 0.363, 0.317, 0.394, 0.445,
    0.435, 0.418, 0.485, 0.378, 0.562, 0.465, 0.388, 0.449, 0.544, 0.569,
    0.491, 0.613, 0.526,
    0.476, 0.531, 0.469, 0.439, 0.419, 0.504, 0.588, 0.583, 0.531, 0.471,
    0.505, 0.493, 0.465, 0.549, 0.374, 0.450, 0.518, 0.537, 0.396, 0.400,
    0.484, 0.342, 0.4545 };

  cout << "\n";
  cout << "TEST08\n";
  cout << "  For samples of a Dirichlet PDF,\n";
  cout << "  DIRICHLET_ESTIMATE estimates the parameters.\n";
  cout << "  DIRICHLET_MEAN finds the means;\n";
  cout << "  DIRICHLET_VARIANCE finds the variances;\n";

  r8mat_print ( sample_num, elem_num, x, "  Sampled data:" );
//
//  Compute the observed averages.
//
  mean = r8col_mean ( sample_num, elem_num, x );

  variance = r8col_variance ( sample_num, elem_num, x );

  cout << "\n";
  cout << "  Observed means, variances are:\n";
  cout << "\n";
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    cout << setw(6) << elem_i
         << setw(14) << mean[elem_i]
         << setw(14) << variance[elem_i] << "\n";
  }

  init = 1;
  alpha = new double[elem_num];
  g = new double[elem_num];
  v = new double[elem_num*elem_num];

  dirichlet_estimate ( elem_num, sample_num, x, sample_num, 
    init, alpha, rlogl, v, g, niter, s, eps, ifault );

  if ( ifault != 0 )
  {
    cout << "\n";
    cout << "WARNING!\n";
    cout << "  DIRICHLET_ESTIMATE error code:\n";
    cout << "  IFAULT = " << ifault << "\n";
  }

  cout << "\n";
  cout << "  Index, Estimate, Lower Limit, Upper Limit:\n";
  cout << "\n";

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = alpha[elem_i] - 1.96 * sqrt ( vari );
    aplus = alpha[elem_i] + 1.96 * sqrt ( vari );
    cout << setw(6) << elem_i
         << setw(14) << alpha[elem_i]
         << setw(14) << aminus
         << setw(14) << aplus << "\n";
  }

  delete [] mean;
  delete [] variance;

  mean = dirichlet_mean ( elem_num, alpha );

  variance = dirichlet_variance ( elem_num, alpha );

  cout << "\n";
  cout << "  Expected means, variances are:\n";
  cout << "\n";
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    cout << setw(6) << elem_i
         << setw(14) << mean[elem_i]
         << setw(14) << variance[elem_i] << "\n";
  }

  alpha_sum = r8vec_sum ( elem_num, alpha );

  cout << "\n";
  cout << "  Alpha sum is " << alpha_sum << "\n";
  cout << "\n";
  cout << "  NORMALIZED VALUES:\n";
  cout << "  Index, Estimate, Lower Limit, Upper Limit:\n";
  cout << "\n";

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = ( alpha[elem_i] - 1.96 * sqrt ( vari ) ) / alpha_sum;
    aplus = ( alpha[elem_i] + 1.96 * sqrt ( vari ) ) / alpha_sum;
    cout << setw(6) << elem_i
         << setw(14) << alpha[elem_i] / alpha_sum
         << setw(14) << aminus
         << setw(14) << aplus << "\n";
  }

  cout << "\n";
  cout << "  Log likelikhood function = " << rlogl << "\n";

  delete [] alpha;
  delete [] g;
  delete [] mean;
  delete [] v;
  delete [] variance;

  return;
# undef ELEM_NUM
# undef SAMPLE_NUM
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests GAMMA_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int rep;
  int rep_num = 5;
  int seed;
  int test;
  int test_num = 5;
  double x;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  GAMMA_SAMPLE samples a Gamma distribution.\n";

  seed = 123456789;

  for ( test = 1; test <= test_num; test++ )
  {
    a = r8_uniform_ab ( 0.1, 2.0, seed );
    b = r8_uniform_ab ( 0.1, 2.0, seed );
    cout << "\n";
    cout << "  A = " << a << ", B = " << b << "\n";
    for ( rep = 1; rep <= rep_num; rep++ )
    {
      x = gamma_sample ( a, b, seed );
      cout << "  " << setw(2) << rep
           << "  " << setw(14) << x << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST09 tests DIRICHLET_ESTIMATE, _MEAN, _VARIANCE, _SAMPLE.
//
//  Discussion:
//
//    Data is generated by sampling a distribution with known parameters.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
# define ELEM_NUM 3
# define SAMPLE_NUM 1000

  double alpha[ELEM_NUM] = { 3.22, 20.38, 21.68 };
  double alpha_sum;
  double aminus;
  double aplus;
  int elem_i;
  int elem_num = ELEM_NUM;
  double eps;
  double *g;
  int ifault;
  int init;
  double *mean;
  int niter;
  double rlogl;
  double s;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  int seed;
  double *v;
  double vari;
  double *variance;
  double *x_sample;
  double *x;

  seed = 123456789;

  cout << "\n";
  cout << "TEST09\n";
  cout << "  For a Dirichlet distribution,\n";
  cout << "  DIRICHLET_SAMPLE samples;\n";
  cout << "  DIRICHLET_MEAN finds the means;\n";
  cout << "  DIRICHLET_VARIANCE finds the variances;\n";
  cout << "  DIRICHLET_ESTIMATE estimates the parameters.\n";
//
//  Report.
//
  r8vec_print ( elem_num, alpha, "  Distribution parameters:" );

  mean = dirichlet_mean ( elem_num, alpha );

  variance = dirichlet_variance ( elem_num, alpha );

  cout << "\n";
  cout << "  Distribution means, variances are:\n";
  cout << "\n";
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    cout << setw(6) << elem_i
         << setw(14) << mean[elem_i]
         << setw(14) << variance[elem_i] << "\n";
  }
//
//  Sample the distribution.
//
  x_sample = new double[sample_num*elem_num];

  cout << "\n";
  cout << "  Number of samples is " << sample_num << "\n";

  for ( sample_i = 0; sample_i < sample_num; sample_i++ )
  {
    x = dirichlet_sample ( elem_num, alpha, seed );

    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      x_sample[sample_i+elem_i*sample_num] = x[elem_i];
    }
    delete [] x;
  }
//
//  Print some results.
//
  cout << "\n";
  cout << "  First few samples:\n";
  cout << "\n";

  for ( sample_i = 0; sample_i < i4_min ( sample_num, 10 ); sample_i++ )
  {
    cout << setw(6) << sample_i;
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      cout << setw(14) << x_sample[sample_i+elem_i*sample_num];
    }
    cout << "\n";
  }
//
//  Compute means, variances.
//
  delete [] mean;
  delete [] variance;

  mean = r8col_mean ( sample_num, elem_num, x_sample );

  variance = r8col_variance ( sample_num, elem_num, x_sample );

  cout << "\n";
  cout << "  Observed means, variances are:\n";
  cout << "\n";
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    cout << setw(6) << elem_i
         << setw(14) << mean[elem_i]
         << setw(14) << variance[elem_i] << "\n";
  }
//
//  Destroy the values of ALPHA.
//
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    alpha[elem_i] = 0.0;
  }
//
//  Try to recover the values of ALPHA.
//
  init = 1;
  v = new double[elem_num*elem_num];
  g = new double[elem_num];

  dirichlet_estimate ( elem_num, sample_num, x_sample, sample_num, 
    init, alpha, rlogl, v, g, niter, s, eps, ifault );

  if ( ifault != 0 )
  {
    cout << "\n";
    cout << "Warning!\n";
    cout << "  DIRICHLET_ESTIMATE error code:\n";
    cout << "  IFAULT = " << ifault << "\n";
  }

  cout << "\n";
  cout << "  Index, Estimate, Lower Limit, Upper Limit:\n";
  cout << "\n";

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = alpha[elem_i] - 1.96 * sqrt ( vari );
    aplus = alpha[elem_i] + 1.96 * sqrt ( vari );
    cout << setw(6) << elem_i
         << setw(14) << alpha[elem_i]
         << setw(14) << aminus
         << setw(14) << aplus << "\n";
  }

  alpha_sum = r8vec_sum ( elem_num, alpha );

  cout << "\n";
  cout << "  Alpha sum is " << alpha_sum << "\n";
  cout << "\n";
  cout << "  NORMALIZED VALUES:\n";
  cout << "  Index, Estimate, Lower Limit, Upper Limit:\n";
  cout << "\n";

  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    vari = v[elem_i+elem_i*elem_num];
    aminus = ( alpha[elem_i] - 1.96 * sqrt ( vari ) ) / alpha_sum;
    aplus = ( alpha[elem_i] + 1.96 * sqrt ( vari ) ) / alpha_sum;
    cout << setw(6) << elem_i
         << setw(14) << alpha[elem_i] / alpha_sum
         << setw(14) << aminus
         << setw(14) << aplus << "\n";
  }

  cout << "\n";
  cout << "  Log likelikhood function = " << rlogl << "\n";

  delete [] mean;
  delete [] v;
  delete [] variance;
  delete [] x_sample;

  return;
# undef ELEM_NUM
# undef SAMPLE_NUM
}
//****************************************************************************80

void test10 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST10 tests DIRICHLET_MIX_SAMPLE, _MEAN, _VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 June 2013
//
//  Author:
//
//    John Burkardt
//
{
# define COMP_NUM 3
# define COMP_MAX 3
# define ELEM_NUM 3
# define SAMPLE_NUM 200

  double a[ELEM_NUM];
  double alpha[COMP_MAX*ELEM_NUM] = {
    0.05, 0.85, 0.00,
    0.20, 0.10, 0.50,
    0.75, 0.05, 0.50 };
  int comp_max = COMP_MAX;
  int comp_num = COMP_NUM;
  int *comp_sample;
  int comp_i;
  double comp_weight[COMP_NUM] = { 3.0, 2.0, 1.0 };
  int elem_i;
  int elem_num = ELEM_NUM;
  double *mean;
  int sample_i;
  int sample_num = SAMPLE_NUM;
  int seed;
  double *variance;
  double *x;
  double *x_sample;

  seed = 123456789;

  cout << "\n";
  cout << "TEST10\n";
  cout << "  For a Dirichlet mixture distribution,\n";
  cout << "  DIRICHLET_MIX_SAMPLE samples;\n";
  cout << "  DIRICHLET_MIX_MEAN computes means;\n";
  cout << "  DIRICHLET_MIX_VARIANCE computes variances.\n";
//
//  Report.
//
  r8vec_print ( comp_num, comp_weight, "  Component weight:" );

  cout << "\n";
  cout << "  Component  Parameters Means Variances\n";
  for ( comp_i = 0; comp_i < comp_num; comp_i++ )
  {
    cout << "\n";
    cout << setw(6) << comp_i;
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      a[elem_i] = alpha[comp_i+elem_i*comp_max];
    }
    mean = dirichlet_mean ( elem_num, a );
    variance = dirichlet_variance ( elem_num, a );
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      cout << setw(6) << elem_i
           << "  " << setw(10) << alpha[comp_i+elem_i*comp_max]
           << "  " << setw(10) << mean[elem_i]
           << "  " << setw(10) << variance[elem_i] << "\n";
    }
    delete [] mean;
    delete [] variance;
  }

  mean = dirichlet_mix_mean ( comp_max, comp_num, elem_num, alpha, 
    comp_weight );

  r8vec_print ( elem_num, mean, "  Element means:" );
  delete [] mean;
//
//  Sample the distribution.
//
  comp_sample = new int[sample_num];
  x_sample = new double[elem_num*sample_num];
  cout << "\n";
  cout << "  Number of samples is " << sample_num << "\n";

  for ( sample_i = 0; sample_i < sample_num; sample_i++ )
  {
    x = dirichlet_mix_sample ( comp_max, comp_num, elem_num, alpha, 
      comp_weight, seed, comp_i );

    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      x_sample[elem_i+sample_i*elem_num] = x[elem_i];
    }

    comp_sample[sample_i] = comp_i;

    delete [] x;
  }
//
//  Print some results.
//
  cout << "\n";
  cout << "  First few samples:\n";
  cout << "\n";
  cout << "  Sample  Component  X\n";
  cout << "\n";

  for ( sample_i = 0; sample_i < i4_min ( sample_num, 10 ); sample_i++ )
  {
    cout << "  " << setw(2) << sample_i
         << "  " << setw(2) << comp_sample[sample_i];
    for ( elem_i = 0; elem_i < elem_num; elem_i++ )
    {
      cout << "  " << setw(10) << x_sample[elem_i+sample_i*elem_num];
    }
    cout << "\n";
  }
  delete [] comp_sample;
  delete [] x_sample;
//
//  Compute the observed averages.
//
  mean = r8col_mean ( sample_num, elem_num, x_sample );

  variance = r8col_variance ( sample_num, elem_num, x_sample );

  cout << "\n";
  cout << "  Element  Observed mean, variance\n";
  cout << "\n";
  for ( elem_i = 0; elem_i < elem_num; elem_i++ )
  {
    cout << setw(6) << elem_i
         << "  " << setw(10) << mean[elem_i]
         << "  " << setw(10) << variance[elem_i] << "\n";
  }

  delete [] mean;
  delete [] variance;

  return;
# undef COMP_MAX
# undef COMP_NUM
# undef ELEM_NUM
# undef SAMPLE_NUM
}

