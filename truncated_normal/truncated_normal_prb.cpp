# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "truncated_normal.hpp"

int main ( );

void i4_uniform_ab_test ( );

void normal_01_cdf_test ( );
void normal_01_cdf_inv_test ( );
void normal_01_mean_test ( );
void normal_01_moment_test ( );
void normal_01_pdf_test ( );
void normal_01_sample_test ( );
void normal_01_variance_test ( );

void normal_ms_cdf_test ( );
void normal_ms_cdf_inv_test ( );
void normal_ms_mean_test ( );
void normal_ms_moment_test ( );
void normal_ms_moment_central_test ( );
void normal_ms_pdf_test ( );
void normal_ms_sample_test ( );
void normal_ms_variance_test ( );

void r8_choose_test ( );
void r8_factorial2_test ( );
void r8_mop_test ( );
void r8_uniform_01_test ( );

void r8poly_print_test ( );
void r8poly_value_horner_test ( );

void r8vec_linspace_new_test ( );
void r8vec_print_test ( );

void truncated_normal_a_cdf_test ( );
void truncated_normal_a_cdf_inv_test ( );
void truncated_normal_a_mean_test ( );
void truncated_normal_a_moment_test ( );
void truncated_normal_a_pdf_test ( );
void truncated_normal_a_sample_test ( );
void truncated_normal_a_variance_test ( );

void truncated_normal_ab_cdf_test ( );
void truncated_normal_ab_cdf_inv_test ( );
void truncated_normal_ab_mean_test ( );
void truncated_normal_ab_moment_test ( );
void truncated_normal_ab_pdf_test ( );
void truncated_normal_ab_sample_test ( );
void truncated_normal_ab_variance_test ( );

void truncated_normal_b_cdf_test ( );
void truncated_normal_b_cdf_inv_test ( );
void truncated_normal_b_mean_test ( );
void truncated_normal_b_moment_test ( );
void truncated_normal_b_pdf_test ( );
void truncated_normal_b_sample_test ( );
void truncated_normal_b_variance_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRUNCATED_NORMAL_PRB.
//
//  Discussion:
//
//    TRUNCATED_NORMAL_PRB tests the TRUNCATED_NORMAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "TRUNCATED_NORMAL_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the TRUNCATED_NORMAL library.\n";
//
//  Utilities.
//
  i4_uniform_ab_test ( );
  r8_choose_test ( );
  r8_factorial2_test ( );
  r8_mop_test ( );
  r8_uniform_01_test ( );
  r8poly_print_test ( );
  r8poly_value_horner_test ( ); 
  r8vec_linspace_new_test ( );
  r8vec_print_test ( );
//
//  Library.
//
  normal_01_cdf_test ( );
  normal_01_cdf_inv_test ( );
  normal_01_mean_test ( );
  normal_01_moment_test ( );
  normal_01_pdf_test ( );
  normal_01_sample_test ( );
  normal_01_variance_test ( );

  normal_ms_cdf_test ( );
  normal_ms_cdf_inv_test ( );
  normal_ms_mean_test ( );
  normal_ms_moment_test ( );
  normal_ms_moment_central_test ( );
  normal_ms_pdf_test ( );
  normal_ms_sample_test ( );
  normal_ms_variance_test ( );

  truncated_normal_a_cdf_test ( );
  truncated_normal_a_cdf_inv_test ( );
  truncated_normal_a_mean_test ( );
  truncated_normal_a_moment_test ( );
  truncated_normal_a_pdf_test ( );
  truncated_normal_a_sample_test ( );
  truncated_normal_a_variance_test ( );

  truncated_normal_ab_cdf_test ( );
  truncated_normal_ab_cdf_inv_test ( );
  truncated_normal_ab_mean_test ( );
  truncated_normal_ab_moment_test ( );
  truncated_normal_ab_pdf_test ( );
  truncated_normal_ab_sample_test ( );
  truncated_normal_ab_variance_test ( );

  truncated_normal_b_cdf_test ( );
  truncated_normal_b_cdf_inv_test ( );
  truncated_normal_b_mean_test ( );
  truncated_normal_b_moment_test ( );
  truncated_normal_b_pdf_test ( );
  truncated_normal_b_sample_test ( );
  truncated_normal_b_variance_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "TRUNCATED_NORMAL_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void i4_uniform_ab_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    I4_UNIFORM_AB_TEST tests I4_UNIFORM_AB.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 October 2014
//
//  Author:
//
//    John Burkardt
//
{
  int a = -100;
  int b = 200;
  int i;
  int j;
  int seed = 123456789;

  cout << "\n";
  cout << "I4_UNIFORM_AB_TEST\n";
  cout << "  I4_UNIFORM_AB computes pseudorandom values\n";
  cout << "  in an interval [A,B].\n";

  cout << "\n";
  cout << "  The lower endpoint A = " << a << "\n";
  cout << "  The upper endpoint B = " << b << "\n";
  cout << "  The initial seed is " << seed << "\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    j = i4_uniform_ab ( a, b, seed );

    cout << "  " << setw(8) << i
         << "  " << setw(8) << j << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_cdf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_TEST tests NORMAL_01_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double cdf1;
  double cdf2;
  int n_data;
  double x;

  cout << "\n";
  cout << "NORMAL_01_CDF_TEST\n";
  cout << "  NORMAL_01_CDF evaluates the Normal 01 CDF;\n";
  cout << "\n";
  cout << "       X              CDF                       CDF\n";
  cout << "                     (exact)                   (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x, cdf1 );

    if ( n_data == 0 )
    {
      break;
    }

    cdf2 = normal_01_cdf ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(24) << cdf1
         << "  " << setw(24) << cdf2 << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_cdf_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_INV_TEST tests NORMAL_01_CDF_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int n_data;
  double x1;
  double x2;

  cout << "\n";
  cout << "NORMAL_01_CDF_INV_TEST\n";
  cout << "  NORMAL_01_CDF_INV inverts the Normal 01 CDF;\n";
  cout << "\n";
  cout << "      CDF             X                         X\n";
  cout << "                     (exact)                   (computed)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    normal_01_cdf_values ( n_data, x1, cdf );

    if ( n_data == 0 )
    {
      break;
    }

    x2 = normal_01_cdf_inv ( cdf );

    cout << "  " << setw(14) << cdf
         << "  " << setw(24) << x1
         << "  " << setw(24) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_MEAN_TEST tests NORMAL_01_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mean;
  int sample_num;
  int seed = 123456789;
  double *x;
  double xmax;
  double xmin;

  cout << "\n";
  cout << "NORMAL_01_MEAN_TEST\n";
  cout << "  NORMAL_01_MEAN computes the Normal 01 mean;\n";

  mean = normal_01_mean ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";

  sample_num = 1000;
  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_01_sample ( seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void normal_01_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_MOMENT_TEST tests NORMAL_01_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double moment;
  int order;

  cout << "\n";
  cout << "NORMAL_01_MOMENT_TEST\n";
  cout << "  NORMAL_01_MOMENT evaluates Normal 01 moments;\n";
  cout << "\n";
  cout << "      Order              Moment\n";
  cout << "\n";

  for ( order = 0; order <= 10; order++ )
  {
    moment = normal_01_moment ( order );

    cout << "  " << setw(14) << order
         << "  " << setw(14) << moment << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_pdf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_PDF_TEST tests NORMAL_01_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double pdf;
  double x;

  cout << "\n";
  cout << "NORMAL_01_PDF_TEST\n";
  cout << "  NORMAL_01_PDF evaluates the Normal 01 PDF;\n";
  cout << "\n";
  cout << "       X              PDF\n";
  cout << "\n";

  for ( i = - 20; i <= 20; i++ )
  {
    x = ( double ) ( i ) / 10.0;

    pdf = normal_01_pdf ( x );

    cout << "  " << setw(14) << x
         << "  " << setw(14) << pdf << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_sample_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_SAMPLE_TEST tests NORMAL_01_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  double x;

  cout << "\n";
  cout << "NORMAL_01_SAMPLE_TEST\n";
  cout << "  NORMAL_01_SAMPLE returns samples from the normal\n";
  cout << "  distribution with mean 0 and standard deviation 1.\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_01_sample ( seed );
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x << "\n";
  }

  return;
}
//****************************************************************************80

void normal_01_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_VARIANCE_TEST tests NORMAL_01_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int sample_num;
  int seed = 123456789;
  double variance;
  double *x;

  cout << "\n";
  cout << "NORMAL_01_VARIANCE_TEST\n";
  cout << "  NORMAL_01_VARIANCE computes the Normal 01 variance;\n";

  variance = normal_01_variance ( );

  cout << "\n";
  cout << "  PDF variance = " << variance << "\n";

  sample_num = 1000;

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_01_sample ( seed );
  }

  variance = r8vec_variance ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num  << "\n";
  cout << "  Sample variance = " << variance << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void normal_ms_cdf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_CDF_TEST tests NORMAL_MS_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double mu;
  double sigma;
  double x;

  cout << "\n";
  cout << "NORMAL_MS_CDF_TEST\n";
  cout << "  NORMAL_MS_CDF evaluates the Normal MS CDF;\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";
  cout << "\n";
  cout << "       X              CDF\n";
  cout << "\n";

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    cdf = normal_ms_cdf ( x, mu, sigma );
    cout << "  " << setw(14) << x
         << "  " << setw(24) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void normal_ms_cdf_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_CDF_INV_TEST tests NORMAL_MS_CDF_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double mu;
  double sigma;
  double x;
  double x2;

  cout << "\n";
  cout << "NORMAL_MS_CDF_INV_TEST\n";
  cout << "  NORMAL_MS_CDF_INV inverts the Normal MS CDF;\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";

  cout << "\n";
  cout << "       X            CDF           CDF_INV\n";
  cout << "\n";

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    cdf = normal_ms_cdf ( x, mu, sigma );
    x2 = normal_ms_cdf_inv ( cdf, mu, sigma );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << cdf
         << "  " << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void normal_ms_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_MEAN_TEST tests NORMAL_MS_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mean;
  double mu;
  int sample_num;
  int seed = 123456789;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  cout << "\n";
  cout << "NORMAL_MS_MEAN_TEST\n";
  cout << "  NORMAL_MS_MEAN computes the Normal MS mean.\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";

  mean = normal_ms_mean ( mu, sigma );

  cout << "\n";
  cout << "  PDF mean = " << mean << "\n";

  sample_num = 1000;
  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_ms_sample ( mu, sigma, seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void normal_ms_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MOMENT_MS_TEST tests NORMAL_MS_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  double moment1;
  double moment2;
  double mu;
  double mu_test[4] = { 0.0, 2.0, 10.0, 0.0 };
  int order;
  double sigma;
  double sigma_test[4] = { 1.0, 1.0, 2.0, 2.0 };
  int test;
  int test_num = 4;
 
  cout << "\n";
  cout << "NORMAL_MOMENT_MS_TEST\n";
  cout << "  NORMAL_MS_MOMENT evaluates the moments of the Normal MS distribution.\n";

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    cout << "\n";
    cout << "  Mu = " << mu
         << "  Sigma = " << sigma << "\n";
    cout << " Order  Moment\n";
    cout << "\n";
    for ( order = 0; order <= 8; order++ )
    {
      moment1 = normal_ms_moment ( order, mu, sigma );
      moment2 = normal_ms_moment_values ( order, mu, sigma );
      cout << "  " << setw(2) << order
           << "  " << setw(14) << moment1 
           << "  " << setw(14) << moment2 << "\n";
    }
  }

  return;
}
//****************************************************************************80

void normal_ms_moment_central_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_MOMENT_CENTRAL_TEST tests NORMAL_MS_MOMENT_CENTRAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  double moment1;
  double moment2;
  double mu;
  double mu_test[4] = { 0.0, 2.0, 10.0, 0.0 };
  int order;
  double sigma;
  double sigma_test[4] = { 1.0, 1.0, 2.0, 2.0 };
  int test;
  int test_num = 4;
 
  cout << "\n";
  cout << "NORMAL_MS_MOMENT_CENTRAL_TEST\n";
  cout << "  NORMAL_MS_MOMENT_CENTRAL evaluates the central moments of the\n";
  cout << "  Normal MS distribution.\n";

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    cout << "\n";
    cout << "  Mu = " << mu
         << "  Sigma = " << sigma << "\n";
    cout << " Order  Moment\n";
    cout << "\n";
    for ( order = 0; order <= 8; order++ )
    {
      moment1 = normal_ms_moment_central ( order, mu, sigma );
      moment2 = normal_ms_moment_central_values ( order, mu, sigma );
      cout << "  " << setw(2) << order
           << "  " << setw(14) << moment1 
           << "  " << setw(14) << moment2 << "\n";
    }
  }

  return;
}
//****************************************************************************80

void normal_ms_pdf_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_PDF_TEST tests NORMAL_MS_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mu;
  double pdf;
  double sigma;
  double x;

  cout << "\n";
  cout << "NORMAL_MS_PDF_TEST\n";
  cout << "  NORMAL_MS_PDF evaluates the Normal MS PDF;\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";

  cout << "\n";
  cout << "       X              PDF\n";
  cout << "\n";

  for ( i = - 20; i <= +20; i++ )
  {
    x = mu + sigma * ( double ) ( i ) / 10.0;
    pdf = normal_ms_pdf ( mu, sigma, x );
    cout << "  " << setw(14) << x
         << "  " << setw(24) << pdf << "\n";
  }

  return;
}
//****************************************************************************80

void normal_ms_sample_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_SAMPLE_TEST tests NORMAL_MS_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mu;
  int seed;
  double sigma;
  double x;

  cout << "\n";
  cout << "NORMAL_MS_SAMPLE_TEST\n";
  cout << "  NORMAL_MS_SAMPLE returns samples from the Normal MS PDF.\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";

  cout << "\n";

  seed = 123456789;

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_ms_sample ( mu, sigma, seed );
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x << "\n";
  }

  return;
}
//****************************************************************************80

void normal_ms_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_MS_VARIANCE_TEST tests NORMAL_MS_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double mu;
  int sample_num;
  int seed = 123456789;
  double sigma;
  double variance;
  double *x;

  cout << "\n";
  cout << "NORMAL_MS_VARIANCE_TEST\n";
  cout << "  NORMAL_MS_VARIANCE computes the Normal MS variance;\n";

  mu = 100.0;
  sigma = 15.0;

  cout << "\n";
  cout << "  Parameter MU = " << mu << "\n";
  cout << "  Parameteter SIGMA = " << sigma << "\n";

  variance = normal_ms_variance ( mu, sigma );

  cout << "\n";
  cout << "  PDF variance = " << variance << "\n";

  sample_num = 1000;
  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = normal_ms_sample ( mu, sigma, seed );
  }

  variance = r8vec_variance ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample variance = " << variance << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void r8_choose_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_CHOOSE_TEST tests R8_CHOOSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    26 July 2014
//
//  Author:
//
//    John Burkardt
//
{
  double cnk;
  int k;
  int n;

  cout << "\n";
  cout << "R8_CHOOSE_TEST\n";
  cout << "  R8_CHOOSE evaluates C(N,K).\n";
  cout << "\n";
  cout << "         N         K       CNK\n";
 
  for ( n = 0; n <= 5; n++ )
  {
    cout << "\n";
    for ( k = 0; k <= n; k++ )
    {
      cnk = r8_choose ( n, k );
      cout << setw(10) << n << "  "
           << setw(8) << k << "  "
           << setw(14) << cnk << "\n";
    }
  }
 
  return;
}
//****************************************************************************80

void r8_factorial2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL2_TEST tests R8_FACTORIAL2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 February 2015
//
//  Author:
//
//    John Burkardt
//
{
  double f1;
  double f2;
  int n;
  int n_data;
  streamsize ss;
//
//  Save the current precision.
//
  ss = cout.precision ( );

  cout << "\n";
  cout << "R8_FACTORIAL2_TEST\n";
  cout << "  R8_FACTORIAL2 evaluates the double factorial function.\n";
  cout << "\n";
  cout << "    N                Exact                  Computed\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    r8_factorial2_values ( n_data, n, f1 );

    if ( n_data == 0 )
    {
      break;
    }

    f2 = r8_factorial2 ( n );

    cout << "  "
         << setw(4) << n << "  "
         << setprecision(16) << setw(24) << f1 << "  "
         << setprecision(16) << setw(24) << f2 << "\n";
  }
//
//  Restore the default precision.
//
  cout.precision ( ss );

  return;
}
//****************************************************************************80

void r8_mop_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_MOP_TEST tests R8_MOP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2014
//
//  Author:
//
//    John Burkardt
//
{
  int i4;
  int i4_max;
  int i4_min;
  double r8;
  int seed = 123456789;
  int test;

  cout << "\n";
  cout << "R8_MOP_TEST\n";
  cout << "  R8_MOP evaluates (-1.0)^I4 as an R8.\n";
  cout << "\n";
  cout << "    I4  R8_MOP(I4)\n";
  cout << "\n";

  i4_min = -100;
  i4_max = +100;

  for ( test = 1; test <= 10; test++ )
  {
    i4 = i4_uniform_ab ( i4_min, i4_max, seed );
    r8 = r8_mop ( i4 );
    cout << "  "
         << setw(4) << i4 << "  "
         << setw(4) <<r8 << "\n";
  }

  return;
}
//****************************************************************************80

void r8_uniform_01_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01_TEST tests R8_UNIFORM_01.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 September 2005
//
//  Author:
//
//    John Burkardt
//
{
# define N 1000

  int i;
  double max;
  double mean;
  double min;
  int n;
  int seed = 123456789;
  double x[N];
  double variance;

  cout << "\n";
  cout << "R8_UNIFORM_01_TEST\n";
  cout << "  R8_UNIFORM_01 samples a uniform random distribution in [0,1].\n";
  cout << "  distributed random numbers.\n";
  cout << "  Using initial random number seed = " << seed << "\n";

  for ( i = 0; i < N; i++ )
  {
    x[i] = r8_uniform_01 ( seed );
  }

  cout << "\n";
  cout << "  First few values:\n";
  cout << "\n";
  for ( i = 0; i < 10; i++ )
  {
    cout << "  " << setw(6) << i
         << "  " << setw(14) << x[i] << "\n";
  }
  min = r8vec_min ( N, x );
  max = r8vec_max ( N, x );
  mean = r8vec_mean ( N, x );
  variance = r8vec_variance ( N, x );

  cout << "\n";
  cout << "  Number of samples was " << N << "\n";
  cout << "  Minimum value was " << min << "\n";
  cout << "  Maximum value was " << max << "\n";
  cout << "  Average value was " << mean << "\n";
  cout << "  Variance was      " << variance << "\n";

  return;
# undef N
}
//****************************************************************************80

void r8poly_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_PRINT_TEST tests R8POLY_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[6] = { 2.0, -3.4, 56.0, 0.0, 0.78, 9.0 };
  int m = 5;

  cout << "\n";
  cout << "R8POLY_PRINT_TEST\n";
  cout << "  R8POLY_PRINT prints an R8POLY.\n";

  r8poly_print ( m, c, "  The R8POLY:" );

  return;
}
//****************************************************************************80

void r8poly_value_horner_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8POLY_VALUE_HORNER_TEST tests R8POLY_VALUE_HORNER.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 January 2015
//
//  Author:
//
//    John Burkardt
//
{
  double c[5] = { 24.0, -50.0, +35.0, -10.0, 1.0 };
  int i;
  int m = 4;
  int n = 16;
  double p;
  double *x;
  double x_hi;
  double x_lo;

  cout << "\n";
  cout << "R8POLY_VALUE_HORNER_TEST\n";
  cout << "  R8POLY_VALUE_HORNER evaluates a polynomial at\n";
  cout << "  one point, using Horner's method.\n";

  r8poly_print ( m, c, "  The polynomial coefficients:" );

  x_lo = 0.0;
  x_hi = 5.0;
  x = r8vec_linspace_new ( n, x_lo, x_hi );

  cout << "\n";
  cout << "   I    X    P(X)\n";
  cout << "\n";

  for ( i = 0; i < n; i++ )
  {
    p = r8poly_value_horner ( m, c, x[i] );
    cout << "  " << setw(2) << i
         << "  " << setw(8) << x[i]
         << "  " << setw(14) << p << "\n";
  }

  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_linspace_new_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW_TEST tests R8VEC_LINSPACE_NEW.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 June 2012
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int n = 5;
  double *x;

  cout << "\n";
  cout << "R8VEC_LINSPACE_NEW_TEST\n";
  cout << "  For a R8VEC:\n";
  cout << "  R8VEC_LINSPACE_NEW: evenly spaced points between A and B;\n";

  a = 10.0;
  b = 20.0;

  x = r8vec_linspace_new ( n, a, b );
  r8vec_print ( n, x, "  r8vec_linspace ( 5, 10, 20 )" );
  delete [] x;

  return;
}
//****************************************************************************80

void r8vec_print_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_PRINT_TEST tests R8VEC_PRINT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[4] = { 123.456, 0.000005, -1.0E+06, 3.14159265 };
  int n = 4;

  cout << "\n";
  cout << "TEST1335\n";
  cout << "  R8VEC_PRINT prints an R8VEC.\n";

  r8vec_print ( n, a, "  The R8VEC:" );

  return;
}
//****************************************************************************80

void truncated_normal_a_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_A_CDF_TEST tests TRUNCATED_NORMAL_A_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_CDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_A_CDF evaluates\n";
  cout << "  the lower Truncated Normal Cumulative Density Function.\n";
  cout << "\n";
  cout << "        MU       S         A         X        CDF1           CDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_cdf_values ( n_data, mu, sigma, a, x, cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_a_cdf ( x, mu, sigma, a );

    cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << a
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << cdf1
         << "  " << setprecision(16) << setw(24) << cdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_a_cdf_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_CDF_INV_TEST tests TRUNCATED_NORMAL_A_CDF_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;
  double x2;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_CDF_INV_TEST:\n";
  cout << "  TRUNCATED_NORMAL_A_CDF_INV inverts the CDF of\n";
  cout << "  the lower Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << ",+oo)\n";
  cout << "\n";
  cout << "       X            CDF           CDF_INV\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_a_sample ( mu, sigma, a, seed );
    cdf = truncated_normal_a_cdf ( x, mu, sigma, a );
    x2 = truncated_normal_a_cdf_inv ( cdf, mu, sigma, a );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x
         << "  " << setw(14) << cdf
         << "  " << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_a_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_MEAN_TEST tests TRUNCATED_NORMAL_A_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_MEAN_TEST\n";
  cout << "  TRUNCATED_NORMAL_A_MEAN computes the mean\n";
  cout << "  of the lower Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << ",+oo)\n";

  mean = truncated_normal_a_mean ( mu, sigma, a );

  cout << "\n";
  cout << "  PDF mean = " << mean << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_a_sample ( mu, sigma, a, seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void truncated_normal_a_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_MOMENT_TEST tests TRUNCATED_NORMAL_A_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double a_test[6] = {
    0.0, -10.0, 10.0, -10.0, 10.0, -10.0 };
  double moment;
  double mu;
  double mu_test[6] = {
    0.0,  0.0,  0.0,  0.0, 0.0, -5.0 };
  int order;
  double sigma;
  double sigma_test[6] = {
    1.0,  1.0,  1.0, 2.0, 2.0, 1.0 };
  int test;
  int test_num;

  test_num = 6;
 
  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_MOMENT_TEST\n";
  cout << "  TRUNCATED_NORMAL_A_MOMENT evaluates the moments\n";
  cout << "  of the Lower Truncated Normal Distribution.\n";

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    a = a_test[test];
    cout << "\n";
    cout << "  Test = " << test 
         << ", Mu = " << mu 
         << ", Sigma = " << sigma 
         << ", A = " << a << "\n";
    cout << " Order  Moment\n";
    cout << "\n";
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_a_moment ( order, mu, sigma, a );
      cout << "  " << setw(2) << order
           << "  " << setw(14) << moment << "\n";
    }
  }
  return;
}
//****************************************************************************80

void truncated_normal_a_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_A_PDF_TEST tests TRUNCATED_NORMAL_A_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_PDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_A_PDF evaluates the PDF of\n";
  cout << "  the lower Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "        MU       S         A         X        PDF1        PDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_a_pdf_values ( n_data, mu, sigma, a, x, pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_a_pdf ( x, mu, sigma, a );

   cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << a
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << pdf1
         << "  " << setprecision(16)<< setw(24) << pdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_a_sample_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_SAMPLE_TEST tests TRUNCATED_NORMAL_A_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_SAMPLE_TEST:\n";
  cout << "  TRUNCATED_NORMAL_A_SAMPLE samples\n";
  cout << "  the lower Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << ",+oo)\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_a_sample ( mu, sigma, a, seed );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_a_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_VARIANCE_TEST tests TRUNCATED_NORMAL_A_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;

  a = 50.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_A_VARIANCE_TEST\n";
  cout << "  TRUNCATED_NORMAL_A_VARIANCE computes the variance\n";
  cout << "  of the lower Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << ",+oo)\n";

  variance = truncated_normal_a_variance ( mu, sigma, a );

  cout << "\n";
  cout << "  PDF variance = " << variance << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_a_sample ( mu, sigma, a, seed );
  }

  variance = r8vec_variance ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample variance = " << variance << "\n";
 
  delete [] x;

  return;
}
//****************************************************************************80

void truncated_normal_ab_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_AB_CDF_TEST tests TRUNCATED_NORMAL_AB_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_CDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_AB_CDF evaluates\n";
  cout << "  the Truncated Normal Cumulative Density Function.\n";
  cout << "\n";
  cout << "        MU       S         A         B         X        CDF1           CDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_cdf_values ( n_data, mu, sigma, a, b, x, cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_ab_cdf ( x, mu, sigma, a, b );

    cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << cdf1
         << "  " << setprecision(16) << setw(24) << cdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_ab_cdf_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_CDF_INV_TEST tests TRUNCATED_NORMAL_AB_CDF_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double cdf;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;
  double x2;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_CDF_INV_TEST:\n";
  cout << "  TRUNCATED_NORMAL_AB_CDF_INV inverts the CDF of\n";
  cout << "  the Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << "," << b << "]\n";
  cout << "\n";
  cout << "       X            CDF           CDF_INV\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_ab_sample ( mu, sigma, a, b, seed );
    cdf = truncated_normal_ab_cdf ( x, mu, sigma, a, b );
    x2 = truncated_normal_ab_cdf_inv ( cdf, mu, sigma, a, b );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x
         << "  " << setw(14) << cdf
         << "  " << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_ab_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_MEAN_TEST tests TRUNCATED_NORMAL_AB_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_MEAN_TEST\n";
  cout << "  TRUNCATED_NORMAL_AB_MEAN computes the mean\n";
  cout << "  of the Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << "," << b << "]\n";

  mean = truncated_normal_ab_mean ( mu, sigma, a, b );

  cout << "\n";
  cout << "  PDF mean = " << mean << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_ab_sample ( mu, sigma, a, b, seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void truncated_normal_ab_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_MOMENT_TEST tests TRUNCATED_NORMAL_AB_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double a_test[9] = {
    -1.0, 0.0, -1.0, -1.0,  0.0, 0.5, -2.0, -4.0, 4.0  };
  double b;
  double b_test[9] = {
    1.0, 1.0,  0.0,  1.0,  2.0, 2.0,  2.0,  4.0, 7.0 };
  double moment;
  double mu;
  double mu_test[9] = {
    0.0, 0.0,  0.0,  0.0,  1.0, 0.0,  0.0,  0.0, 5.0 };
  int order;
  double sigma;
  double sigma_test[9] = {
    1.0, 1.0,  1.0,  2.0,  1.0, 1.0,  1.0,  1.0, 0.5 };
  int test;
  int test_num;

  test_num = 9;
 
  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_MOMENT_TEST\n";
  cout << "  TRUNCATED_NORMAL_AB_MOMENT evaluates the moments\n";
  cout << "  of the Truncated Normal PDF:\n";

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    a = a_test[test];
    b = b_test[test];
    cout << "\n";
    cout << "  Test = " << test 
         << ", Mu = " << mu 
         << ", Sigma = " << sigma 
         << ", A = " << a
         << ", B = " << b << "\n";
    cout << " Order  Moment\n";
    cout << "\n";
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_ab_moment ( order, mu, sigma, a, b );
      cout << "  " << setw(2) << order
           << "  " << setw(14) << moment << "\n";
    }
  }
  return;
}
//****************************************************************************80

void truncated_normal_ab_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_AB_PDF_TEST tests TRUNCATED_NORMAL_AB_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_PDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_AB_PDF evaluates\n";
  cout << "  the Truncated Normal Probability Density Function.\n";
  cout << "\n";
  cout << "        MU       S         A         B         X        PDF1        PDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_ab_pdf_values ( n_data, mu, sigma, a, b, x, pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_ab_pdf ( x, mu, sigma, a, b );

   cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << a
         << "  " << setw(8) << b
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << pdf1
         << "  " << setprecision(16)<< setw(24) << pdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_ab_sample_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_SAMPLE_TEST tests TRUNCATED_NORMAL_AB_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_SAMPLE_TEST:\n";
  cout << "  TRUNCATED_NORMAL_AB_SAMPLE samples\n";
  cout << "  the Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << "," << b << "]\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_ab_sample ( mu, sigma, a, b, seed );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_ab_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_AB_VARIANCE_TEST tests TRUNCATED_NORMAL_AB_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;

  a = 50.0;
  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_AB_VARIANCE_TEST\n";
  cout << "  TRUNCATED_NORMAL_AB_VARIANCE computes the variance\n";
  cout << "  of the Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval [" << a << "," << b << "]\n";

  variance = truncated_normal_ab_variance ( mu, sigma, a, b );

  cout << "\n";
  cout << "  PDF variance = " << variance << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_ab_sample ( mu, sigma, a, b, seed );
  }

  variance = r8vec_variance ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample variance = " << variance << "\n";
 
  delete [] x;

  return;
}
//****************************************************************************80

void truncated_normal_b_cdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_B_CDF_TEST tests TRUNCATED_NORMAL_B_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double cdf1;
  double cdf2;
  double mu;
  int n_data;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_CDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_B_CDF evaluates\n";
  cout << "  the upper Truncated Normal Cumulative Density Function.\n";
  cout << "\n";
  cout << "        MU       S         B         X        CDF1           CDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_cdf_values ( n_data, mu, sigma, b, x, cdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    cdf2 = truncated_normal_b_cdf ( x, mu, sigma, b );

    cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << b
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << cdf1
         << "  " << setprecision(16) << setw(24) << cdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_b_cdf_inv_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_CDF_INV_TEST tests TRUNCATED_NORMAL_B_CDF_INV.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double cdf;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;
  double x2;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_CDF_INV_TEST:\n";
  cout << "  TRUNCATED_NORMAL_B_CDF_INV inverts the CDF of\n";
  cout << "  the upper Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval (-oo," << b << "]\n";
  cout << "\n";
  cout << "       X            CDF           CDF_INV\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_b_sample ( mu, sigma, b, seed );
    cdf = truncated_normal_b_cdf ( x, mu, sigma, b );
    x2 = truncated_normal_b_cdf_inv ( cdf, mu, sigma, b );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x
         << "  " << setw(14) << cdf
         << "  " << setw(14) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_b_mean_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_MEAN_TEST tests TRUNCATED_NORMAL_B_MEAN.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int i;
  double mean;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double *x;
  double xmax;
  double xmin;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_MEAN_TEST\n";
  cout << "  TRUNCATED_NORMAL_B_MEAN computes the mean\n";
  cout << "  of the upper Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval (-oo," << b << "]\n";

  mean = truncated_normal_b_mean ( mu, sigma, b );

  cout << "\n";
  cout << "  PDF mean = " << mean << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_b_sample ( mu, sigma, b, seed );
  }

  mean = r8vec_mean ( sample_num, x );
  xmax = r8vec_max ( sample_num, x );
  xmin = r8vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void truncated_normal_b_moment_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_MOMENT_TEST tests TRUNCATED_NORMAL_B_MOMENT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double b_test[6] = {
    0.0, 10.0, -10.0, 10.0, -10.0, 10.0 };
  double moment;
  double mu;
  double mu_test[6] = {
    0.0,  0.0,  0.0,  0.0, 0.0, 5.0 };
  int order;
  double sigma;
  double sigma_test[6] = {
    1.0,  1.0,  1.0, 2.0, 2.0, 1.0 };
  int test;
  int test_num;

  test_num = 6;
 
  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_MOMENT_TEST\n";
  cout << "  For the Upper Truncated Normal PDF:\n";
  cout << "  TRUNCATED_NORMAL_B_MOMENT evaluates the moments.\n";

  for ( test = 0; test < test_num; test++ )
  {
    mu = mu_test[test];
    sigma = sigma_test[test];
    b = b_test[test];
    cout << "\n";
    cout << "  Test = " << test 
         << ", Mu = " << mu 
         << ", Sigma = " << sigma 
         << ", B = " << b << "\n";
    cout << " Order  Moment\n";
    cout << "\n";
    for ( order = 0; order <= 8; order++ )
    {
      moment = truncated_normal_b_moment ( order, mu, sigma, b );
      cout << "  " << setw(2) << order
           << "  " << setw(14) << moment << "\n";
    }
  }
  return;
}
//****************************************************************************80

void truncated_normal_b_pdf_test ( )

//****************************************************************************80
//
//  Purpose: 
//
//    TRUNCATED_NORMAL_B_PDF_TEST tests TRUNCATED_NORMAL_B_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    14 September 2013
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double mu;
  int n_data;
  double pdf1;
  double pdf2;
  double sigma;
  double x;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_PDF_TEST:\n";
  cout << "  TRUNCATED_NORMAL_B_PDF evaluates\n";
  cout << "  the upper Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "        MU       S         B         X        PDF1        PDF2\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    truncated_normal_b_pdf_values ( n_data, mu, sigma, b, x, pdf1 );

    if ( n_data == 0 )
    {
      break;
    }
    pdf2 = truncated_normal_b_pdf ( x, mu, sigma, b );

   cout << "  " << setw(8) << mu
         << "  " << setw(8) << sigma
         << "  " << setw(8) << b
         << "  " << setw(8) << x
         << "  " << setprecision(16) << setw(24) << pdf1
         << "  " << setprecision(16)<< setw(24) << pdf2 << "\n";
  }
  return;
}
//****************************************************************************80

void truncated_normal_b_sample_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_SAMPLE_TEST tests TRUNCATED_NORMAL_B_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int i;
  double mu;
  int sample_num = 10;
  int seed;
  double sigma;
  double x;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_SAMPLE_TEST:\n";
  cout << "  TRUNCATED_NORMAL_B_SAMPLE samples\n";
  cout << "  the upper Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval (-oo," << b << "]\n";
  cout << "\n";

  for ( i = 0; i < sample_num; i++ )
  {
    x = truncated_normal_b_sample ( mu, sigma, b, seed );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << x << "\n";
  }

  return;
}
//****************************************************************************80

void truncated_normal_b_variance_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_B_VARIANCE_TEST tests TRUNCATED_NORMAL_B_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  int i;
  double mu;
  int sample_num = 1000;
  int seed;
  double sigma;
  double variance;
  double *x;

  b = 150.0;
  mu = 100.0;
  sigma = 25.0;
  seed = 123456789;

  cout << "\n";
  cout << "TRUNCATED_NORMAL_B_VARIANCE_TEST\n";
  cout << "  TRUNCATED_NORMAL_B_VARIANCE computes the variance\n";
  cout << "  of the upper Truncated Normal Distribution.\n";
  cout << "\n";
  cout << "  The parent normal distribution has\n";
  cout << "    mean =               " << mu << "\n";
  cout << "    standard deviation = " << sigma << "\n";
  cout << "  The parent distribution is truncated to\n";
  cout << "  the interval (-oo," << b << "]\n";

  variance = truncated_normal_b_variance ( mu, sigma, b );

  cout << "\n";
  cout << "  PDF variance = " << variance << "\n";

  x = new double[sample_num];

  for ( i = 0; i < sample_num; i++ )
  {
    x[i] = truncated_normal_b_sample ( mu, sigma, b, seed );
  }

  variance = r8vec_variance ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num << "\n";
  cout << "  Sample variance = " << variance << "\n";
 
  delete [] x;

  return;
}
