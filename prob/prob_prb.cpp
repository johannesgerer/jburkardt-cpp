# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "prob.hpp"

int main ( );

void test001 ( );
void test002 ( );
void test003 ( );
void test004 ( );
void test005 ( );
void test006 ( );
void test007 ( );
void test008 ( );
void test009 ( );

void test010 ( );
void test0105 ( );
void test0106 ( );
void test011 ( );
void test012 ( );
void test013 ( );
void test014 ( );
void test015 ( );
void test016 ( );

void test020 ( );
void test021 ( );
void test022 ( );
void test023 ( );
void test0235 ( );
void test024 ( );
void test025 ( );
void test0251 ( );
void test0252 ( );
void test0253 ( );
void test0254 ( );
void test026 ( );
void test027 ( );
void test0275 ( );
void test0276 ( );
void test028 ( );
void test029 ( );

void test030 ( );
void test031 ( );
void test032 ( );
void test033 ( );
void test034 ( );
void test035 ( );
void test036 ( );
void test037 ( );
void test0375 ( );
void test038 ( );
void test039 ( );
void test0395 ( );

void test040 ( );
void test041 ( );
void test042 ( );
void test043 ( );
void test044 ( );
void test045 ( );
void test046 ( );
void test047 ( );
void test048 ( );
void test049 ( );

void test050 ( );
void test051 ( );
void test052 ( );
void test053 ( );
void test054 ( );
void test055 ( );
void test056 ( );
void test0563 ( );
void test0564 ( );
void test0565 ( );
void test0566 ( );
void test057 ( );
void test058 ( );
void test059 ( );

void test060 ( );
void test061 ( );
void test062 ( );
void test063 ( );
void test064 ( );
void test065 ( );
void test066 ( );
void test067 ( );
void test068 ( );
void test069 ( );

void test070 ( );
void test0705 ( );
void test071 ( );
void test072 ( );
void test073 ( );
void test074 ( );
void test0744 ( );
void test0745 ( );
void test075 ( );
void test076 ( );
void test077 ( );
void test078 ( );
void test079 ( );

void test080 ( );
void test081 ( );
void test082 ( );
void test083 ( );
void test084 ( );
void test085 ( );
void test086 ( );
void test087 ( );
void test088 ( );
void test089 ( );

void test090 ( );
void test091 ( );
void test092 ( );
void test093 ( );
void test094 ( );
void test095 ( );
void test096 ( );
void test0965 ( );
void test097 ( );
void test098 ( );
void test099 ( );

void test100 ( );
void test101 ( );
void test102 ( );
void test103 ( );
void test104 ( );
void test105 ( );
void test106 ( );
void test107 ( );
void test108 ( );
void test109 ( );

void test110 ( );
void test111 ( );
void test112 ( );
void test113 ( );
void test114 ( );
void test1145 ( );
void test1146 ( );
void test115 ( );
void test116 ( );
void test117 ( );
void test118 ( );
void test119 ( );

void test120 ( );
void test123 ( );
void test124 ( );
void test125 ( );
void test126 ( );
void test127 ( );
void test128 ( );
void test129 ( );

void test130 ( );
void test1304 ( );
void test1306 ( );
void test131 ( );
void test132 ( );
void test133 ( );
void test134 ( );
void test1341 ( );
void test1342 ( );
void test1344 ( );
void test135 ( );
void test136 ( );
void test137 ( );
void test138 ( );
void test139 ( );

void test140 ( );
void test141 ( );
void test142 ( );
void test1425 ( );
void test143 ( );
void test144 ( );
void test145 ( );
void test146 ( );
void test147 ( );
void test148 ( );
void test1485 ( );
void test1486 ( );
void test149 ( );

void test150 ( );
void test151 ( );
void test152 ( );
void test153 ( );
void test154 ( );
void test155 ( );
void test1555 ( );
void test156 ( );
void test157 ( );
void test158 ( );
void test159 ( );

void test160 ( );
void test161 ( );
void test162 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for PROB_PRB.
//
//  Discussion:
//
//    PROB_PRB calls sample problems for the PROB routines.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "PROB_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the PROB library.\n";

  test001 ( );
  test002 ( );
  test003 ( );
  test004 ( );
  test005 ( );
  test006 ( );
  test007 ( );
  test008 ( );
  test009 ( );

  test010 ( );
  test0105 ( );
  test0106 ( );
  test011 ( );
  test012 ( );
  test013 ( );
  test014 ( );
  test015 ( );
  test016 ( );

  test020 ( );
  test021 ( );
  test022 ( );
  test023 ( );
  test0235 ( );
  test024 ( );
  test025 ( );
  test0251 ( );
  test0252 ( );
  test0253 ( );
  test0254 ( );
  test026 ( );
  test027 ( );
  test0275 ( );
  test0276 ( );
  test028 ( );
  test029 ( );

  test030 ( );
  test031 ( );
  test032 ( );
  test033 ( );
  test034 ( );
  test035 ( );
  test036 ( );
  test037 ( );
  test0375 ( );
  test038 ( );
  test039 ( );
  test0395 ( );

  test040 ( );
  test041 ( );
  test042 ( );
  test043 ( );
  test044 ( );
  test045 ( );
  test046 ( );
  test047 ( );
  test048 ( );
  test049 ( );

  test050 ( );
  test051 ( );
  test052 ( );
  test053 ( );
  test054 ( );
  test055 ( );
  test056 ( );
  test0563 ( );
  test0564 ( );
  test0565 ( );
  test0566 ( );
  test057 ( );
  test058 ( );
  test059 ( );

  test060 ( );
  test061 ( );
  test062 ( );
  test063 ( );
  test064 ( );
  test065 ( );
  test066 ( );
  test067 ( );
  test068 ( );
  test069 ( );

  test070 ( );
  test0705 ( );
  test071 ( );
  test072 ( );
  test073 ( );
  test074 ( );
  test0744 ( );
  test0745 ( );
  test075 ( );
  test076 ( );
  test077 ( );
  test078 ( );
  test079 ( );

  test080 ( );
  test081 ( );
  test082 ( );
  test083 ( );
  test084 ( );
  test085 ( );
  test086 ( );
  test087 ( );
  test088 ( );
  test089 ( );

  test090 ( );
  test091 ( );
  test092 ( );
  test093 ( );
  test094 ( );
  test095 ( );
  test096 ( );
  test0965 ( );
  test097 ( );
  test098 ( );
  test099 ( );

  test100 ( );
  test101 ( );
  test102 ( );
  test103 ( );
  test104 ( );
  test105 ( );
  test106 ( );
  test107 ( );
  test108 ( );
  test109 ( );

  test110 ( );
  test111 ( );
  test112 ( );
  test113 ( );
  test114 ( );
  test1145 ( );
  test1146 ( );
  test115 ( );
  test116 ( );
  test117 ( );
  test118 ( );
  test119 ( );

  test120 ( );
  test123 ( );
  test124 ( );
  test125 ( );
  test126 ( );
  test127 ( );
  test128 ( );
  test129 ( );

  test130 ( );
  test1304 ( );
  test1306 ( );
  test131 ( );
  test132 ( );
  test133 ( );
  test134 ( );
  test1341 ( );
  test1342 ( );
  test1344 ( );
  test135 ( );
  test136 ( );
  test137 ( );
  test138 ( );
  test139 ( );

  test140 ( );
  test141 ( );
  test142 ( );
  test1425 ( );
  test143 ( );
  test144 ( );
  test145 ( );
  test146 ( );
  test147 ( );
  test148 ( );
  test1485 ( );
  test1486 ( );
  test149 ( );

  test150 ( );
  test151 ( );
  test152 ( );
  test153 ( );
  test154 ( );
  test155 ( );
  test1555 ( );
  test156 ( );
  test157 ( );
  test158 ( );
  test159 ( );

  test160 ( );
  test161 ( );
  test162 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "PROB_PRB\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test001 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST001 tests ANGLE_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int n;
  double x;

  cout << "\n";
  cout << "TEST001\n";
  cout << "  For the ANGLE PDF:\n";
  cout << "  ANGLE_CDF evaluates the CDF;\n";

  n = 5;
  x = 0.50E+00;

  cdf = angle_cdf ( x, n );

  cout << "\n";
  cout << "  Parameter N =     " << n   << "\n";
  cout << "  PDF argument X =   " << x   << "\n";
  cout << "  CDF value =       " << cdf << "\n";

  return;
}
//****************************************************************************80

void test002 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST002 tests ANGLE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double pdf;
  double x;

  cout << "\n";
  cout << "TEST002\n";
  cout << "  For the ANGLE PDF:\n";
  cout << "  ANGLE_PDF evaluates the PDF;\n";

  n = 5;
  x = 0.50E+00;

  pdf = angle_pdf ( x, n );

  cout << "\n";
  cout << "  Parameter N =    " << n   << "\n";
  cout << "  PDF argument X =  " << x   << "\n";
  cout << "  PDF value =      " << pdf << "\n";

  return;
}
//****************************************************************************80

void test003 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST003 tests ANGLE_MEAN;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double mean;
  int n;

  cout << "\n";
  cout << "TEST003\n";
  cout << "  For the ANGLE PDF:\n";
  cout << "  ANGLIT_MEAN computes the mean;\n";

  n = 5;
  mean = angle_mean ( n );

  cout << "\n";
  cout << "  Parameter N = " << n    << "\n";
  cout << "  PDF mean =    " << mean << "\n";

  return;
}
//****************************************************************************80

void test004 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST004 tests ANGLIT_CDF, ANGLIT_CDF_INV, ANGLIT_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST004\n";
  cout << "  For the Anglit PDF:\n";
  cout << "  ANGLIT_CDF evaluates the CDF;\n";
  cout << "  ANGLIT_CDF_INV inverts the CDF.\n";
  cout << "  ANGLIT_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = anglit_sample ( &seed );
    pdf = anglit_pdf ( x );
    cdf = anglit_cdf ( x );
    x2 = anglit_cdf_inv ( cdf );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test005 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST005 tests ANGLIT_MEAN, ANGLIT_SAMPLE, ANGLIT_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST005\n";
  cout << "  For the Anglit PDF:\n";
  cout << "  ANGLIT_MEAN computes the mean;\n";
  cout << "  ANGLIT_SAMPLE samples;\n";
  cout << "  ANGLIT_VARIANCE computes the variance.\n";

  mean     = anglit_mean ( );
  variance = anglit_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = anglit_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test006 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST006 tests ARCSIN_CDF, ARCSIN_CDF_INV, ARCSIN_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST006\n";
  cout << "  For the Arcsin PDF:\n";
  cout << "  ARCSIN_CDF evaluates the CDF;\n";
  cout << "  ARCSIN_CDF_INV inverts the CDF.\n";
  cout << "  ARCSIN_PDF evaluates the PDF;\n";

  a = 1.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !arcsin_check ( a ) )
  {
    cout << "\n";
    cout << "TEST006 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = arcsin_sample ( a, &seed );
    pdf = arcsin_pdf ( x, a );
    cdf = arcsin_cdf ( x, a );
    x2 = arcsin_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test007 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST007 tests ARCSIN_MEAN, ARCSIN_SAMPLE, ARCSIN_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST007\n";
  cout << "  For the Arcsin PDF:\n";
  cout << "  ARCSIN_MEAN computes the mean;\n";
  cout << "  ARCSIN_SAMPLE samples;\n";
  cout << "  ARCSIN_VARIANCE computes the variance.\n";

  for ( i = 1; i <= 2; i++ )
  {
    if ( i == 1 )
    {
      a = 1.0E+00;
    }
    else if ( i == 2 )
    {
      a = 16.0E+00;
    }

    cout << "\n";
    cout << "  PDF parameter A =             " << a        << "\n";

    if ( !arcsin_check ( a ) )
    {
      cout << "\n";
      cout << "TEST007 - Fatal error!\n";
      cout << "  The parameters are not legal.\n";
      return;
    }

    mean = arcsin_mean ( a );
    variance = arcsin_variance ( a );

    cout << "  PDF mean =                    " << mean     << "\n";
    cout << "  PDF variance =                " << variance << "\n";

    for ( j = 0; j < SAMPLE_NUM; j++ )
    {
      x[j] = arcsin_sample ( a, &seed );
    }

    mean = r8vec_mean ( SAMPLE_NUM, x );
    variance = r8vec_variance ( SAMPLE_NUM, x );
    xmax = r8vec_max ( SAMPLE_NUM, x );
    xmin = r8vec_min ( SAMPLE_NUM, x );

    cout << "\n";
    cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
    cout << "  Sample mean =     " << mean     << "\n";
    cout << "  Sample variance = " << variance << "\n";
    cout << "  Sample maximum =  " << xmax     << "\n";
    cout << "  Sample minimum =  " << xmin     << "\n";
  }

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test008 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST008 tests BENFORD_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  double pdf;

  cout << "\n";
  cout << "TEST008\n";
  cout << "  For the Benford PDF:\n";
  cout << "  BENFORD_PDF evaluates the PDF.\n";

  cout << "\n";
  cout << "  N    PDF(N)\n";
  cout << "\n";

  for ( n = 1; n <= 19; n++ )
  {
    pdf = benford_pdf ( n );
    cout                    << "  "
         << setw(6)  << n   << "  "
         << setw(14) << pdf << "\n";
  }

  return;
}
//****************************************************************************80

void test009 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST009 tests BERNOULLI_CDF, BERNOULLI_CDF_INV, BERNOULLI_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST009\n";
  cout << "  For the Bernoulli PDF,\n";
  cout << "  BERNOULLI_CDF evaluates the CDF;\n";
  cout << "  BERNOULLI_CDF_INV inverts the CDF;\n";
  cout << "  BERNOULLI_PDF evaluates the PDF.\n";

  a = 0.75E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !bernoulli_check ( a ) )
  {
    cout << "\n";
    cout << "TEST009 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = bernoulli_sample ( a, &seed );
    pdf = bernoulli_pdf ( x, a );
    cdf = bernoulli_cdf ( x, a );
    x2 = bernoulli_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test010 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST010 tests BERNOULLI_MEAN, BERNOULLI_SAMPLE, BERNOULLI_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST010\n";
  cout << "  For the Bernoulli PDF:\n";
  cout << "  BERNOULLI_MEAN computes the mean;\n";
  cout << "  BERNOULLI_SAMPLE samples;\n";
  cout << "  BERNOULLI_VARIANCE computes the variance.\n";

  a = 0.75E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !bernoulli_check ( a ) )
  {
    cout << "\n";
    cout << "TEST010 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = bernoulli_mean ( a );
  variance = bernoulli_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =       " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = bernoulli_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0105 tests BESSEL_I0_INC, BESSEL_I0_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST0105:\n";
  cout << "  BESSEL_I0 evaluates the Bessel function of the\n";
  cout << "    first kind and order 0;\n";
  cout << "  BESSEL_I0_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      X       Exact F       BESSEL_I0(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i0_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = bessel_i0 ( x );

    cout << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0106 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0106 tests BESSEL_I1_INC, BESSEL_I1_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST0106:\n";
  cout << "  BESSEL_I1 evaluates the Bessel function of the\n";
  cout << "    first kind and order 1;\n";
  cout << "  BESSEL_I1_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      X       Exact F       BESSEL_I1(X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_i1_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = bessel_i1 ( x );

    cout << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test011 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST011 tests BETA, GAMMA;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double beta1;
  double beta2;

  cout << "\n";
  cout << "TEST011\n";
  cout << "  BETA evaluates the Beta function;\n";
  cout << "  GAMMA evaluates the Gamma function.\n";

  a = 2.2;
  b = 3.7;

  beta1 = beta ( a, b );
  beta2 = gamma ( a ) * gamma ( b ) / gamma ( a + b );

  cout << "\n";
  cout << "  Argument A =                   " << a << "\n";
  cout << "  Argument B =                   " << b << "\n";
  cout << "  Beta(A,B) =                    " << beta1 << "\n";
  cout << "  (Expected value = 0.0454 )\n";
  cout << "\n";
  cout << "  Gamma(A)*Gamma(B)/Gamma(A+B) = " << beta2 << "\n";

  return;
}
//****************************************************************************80

void test012 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST012 tests BETA_CDF, BETA_CDF_INV, BETA_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST012\n";
  cout << "  For the Beta PDF:\n";
  cout << "  BETA_CDF evaluates the CDF;\n";
  cout << "  BETA_CDF_INV inverts the CDF.\n";
  cout << "  BETA_PDF evaluates the PDF;\n";

  a = 12.0;
  b = 12.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !beta_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST012 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = beta_sample ( a, b, &seed );
    pdf = beta_pdf ( x, a, b );
    cdf = beta_cdf ( x, a, b );
    x2 = beta_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test013 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST013 tests BETA_INC, BETA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST013:\n";
  cout << "  BETA_INC evaluates the normalized incomplete Beta\n";
  cout << "    function BETA_INC(A,B,X).\n";
  cout << "  BETA_INC_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "         A         B         X       Exact F       BETA_INC(A,B,X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    beta_inc_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = beta_inc ( a, b, x );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test014 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST014 tests BETA_MEAN, BETA_SAMPLE, BETA_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST014\n";
  cout << "  For the Beta PDF:\n";
  cout << "  BETA_MEAN computes the mean;\n";
  cout << "  BETA_SAMPLE samples;\n";
  cout << "  BETA_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !beta_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST014 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = beta_mean ( a, b );
  variance = beta_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = beta_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test015 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST015 tests BETA_BINOMIAL_CDF, BETA_BINOMIAL_CDF_INV, BETA_BINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST015\n";
  cout << "  For the Beta Binomial PDF:\n";
  cout << "  BETA_BINOMIAL_CDF evaluates the CDF;\n";
  cout << "  BETA_BINOMIAL_CDF_INV inverts the CDF.\n";
  cout << "  BETA_BINOMIAL_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;
  c = 4;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !beta_binomial_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST015 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = beta_binomial_sample ( a, b, c, &seed );
    pdf = beta_binomial_pdf ( x, a, b, c );
    cdf = beta_binomial_cdf ( x, a, b, c );
    x2 = beta_binomial_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test016 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST016 tests BETA_BINOMIAL_MEAN, BETA_BINOMIAL_SAMPLE, BETA_BINOMIAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST016\n";
  cout << "  For the Beta Binomial PDF:\n";
  cout << "  BETA_BINOMIAL_MEAN computes the mean;\n";
  cout << "  BETA_BINOMIAL_SAMPLE samples;\n";
  cout << "  BETA_BINOMIAL_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;
  c = 4;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !beta_binomial_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST016 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = beta_binomial_mean ( a, b, c );
  variance = beta_binomial_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = beta_binomial_sample ( a, b, c, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test020 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST020 tests BINOMIAL_CDF, BINOMIAL_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double b;
  double fx;
  double fx2;
  int n_data;
  int x;

  cout << "\n";
  cout << "TEST020:\n";
  cout << "  BINOMIAL_CDF evaluates the cumulative distribution\n";
  cout << "    function for the discrete binomial probability\n";
  cout << "    density function.\n";
  cout << "  BINOMIAL_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  A is the number of trials;\n";
  cout << "  B is the probability of success on one trial;\n";
  cout << "  X is the number of successes;\n";
  cout << "  BINOMIAL_CDF is the probability of having up to X\n";
  cout << "  successes.\n";
  cout << "\n";
  cout << "      A     B         X   Exact F     BINOMIAL_CDF(A,B,X)\n";

  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    binomial_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = binomial_cdf ( x, a, b );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test021 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST021 tests BINOMIAL_CDF, BINOMIAL_CDF_INV, BINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST021\n";
  cout << "  For the Binomial PDF:\n";
  cout << "  BINOMIAL_CDF evaluates the CDF;\n";
  cout << "  BINOMIAL_CDF_INV inverts the CDF.\n";
  cout << "  BINOMIAL_PDF evaluates the PDF;\n";

  a = 5;
  b = 0.65;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !binomial_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST021 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = binomial_sample ( a, b, &seed );
    pdf = binomial_pdf ( x, a, b );
    cdf = binomial_cdf ( x, a, b );
    x2 = binomial_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test022 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST022 tests BINOMIAL_COEF, BINOMIAL_COEF_LOG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int cnk1;
  double cnk2_log;
  double cnk2;
  int k;
  int n;

  cout << "\n";
  cout << "TEST022\n";
  cout << "  BINOMIAL_COEF evaluates binomial coefficients.\n";
  cout << "  BINOMIAL_COEF_LOG evaluates the logarithm.\n";

  cout << "\n";
  cout << "    N     K       C(N,K)\n";
  cout << "\n";

  for ( n = 0; n <= 4; n++ )
  {
    for ( k = 0; k <= n; k++ )
    {
      cnk1 = binomial_coef ( n, k );
      cnk2_log = binomial_coef_log ( n, k );
      cnk2 = exp ( cnk2_log );
      cout << "  "
           << setw(6)  << n    << "  "
           << setw(6)  << k    << "  "
           << setw(6)  << cnk1 << "  "
           << setw(14) << cnk2 << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test023 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST023 tests BINOMIAL_MEAN, BINOMIAL_SAMPLE, BINOMIAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST023\n";
  cout << "  For the Binomial PDF:\n";
  cout << "  BINOMIAL_MEAN computes the mean;\n";
  cout << "  BINOMIAL_SAMPLE samples;\n";
  cout << "  BINOMIAL_VARIANCE computes the variance;\n";

  a = 5;
  b = 0.30;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !binomial_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST023 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = binomial_mean ( a, b );
  variance = binomial_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = binomial_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0235 ( )

//****************************************************************************80*
//
//  Purpose:
//
//    TEST0235 tests BIRTHDAY_CDF, BIRTHDAY_CDF_INV, BIRTHDAY_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int n;
  int n2;
  double pdf;

  cout << "\n";
  cout << "TEST0235\n";
  cout << "  For the Birthday PDF,\n";
  cout << "  BIRTHDAY_CDF evaluates the CDF;\n";
  cout << "  BIRTHDAY_CDF_INV inverts the CDF.\n";
  cout << "  BIRTHDAY_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       N            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( n = 1; n <= 30; n++ )
  {
    pdf = birthday_pdf ( n );

    cdf = birthday_cdf ( n );

    n2 = birthday_cdf_inv ( cdf );

    cout << "  " << setw(8)  << n
         << "  " << setw(14) << pdf
         << "  " << setw(14) << cdf
         << "  " << setw(8)  << n2 << "\n";
  }
  return;
}
//****************************************************************************80

void test024 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST024 tests BRADFORD_CDF, BRADFORD_CDF_INV, BRADFORD_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST024\n";
  cout << "  For the Bradford PDF:\n";
  cout << "  BRADFORD_CDF evaluates the CDF;\n";
  cout << "  BRADFORD_CDF_INV inverts the CDF.\n";
  cout << "  BRADFORD_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !bradford_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST024 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = bradford_sample ( a, b, c, &seed );
    pdf = bradford_pdf ( x, a, b, c );
    cdf = bradford_cdf ( x, a, b, c );
    x2 = bradford_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test025 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST025 tests BRADFORD_MEAN, BRADFORD_SAMPLE, BRADFORD_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST025\n";
  cout << "  For the Bradford PDF:\n";
  cout << "  BRADFORD_MEAN computes the mean;\n";
  cout << "  BRADFORD_SAMPLE samples;\n";
  cout << "  BRADFORD_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !bradford_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST025 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = bradford_mean ( a, b, c );
  variance = bradford_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = bradford_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0251 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0251 tests BUFFON_LAPLACE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  int j;
  int k;
  double l;
  double pdf;

  cout << "\n";
  cout << "TEST0251\n";
  cout << "  BUFFON_LAPLACE_PDF evaluates the Buffon-Laplace PDF, the probability\n";
  cout << "  that, on a grid of cells of width A and height B,\n";
  cout << "  a needle of length L, dropped at random, will cross\n";
  cout << "  at least one grid line.\n";
  cout << "\n";
  cout << "      A         B         L        PDF\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    a = ( double ) ( i );
    for ( j = 1; j <= 5; j++ )
    {
      b = ( double ) ( j );
      for ( k = 0; k <= 5; k++ )
      {
        l = ( double ) ( k ) * r8_min ( a, b ) / 5.0;
        pdf = buffon_laplace_pdf ( a, b, l );
        cout << "  " << setw(8) << a
             << "  " << setw(8) << b
             << "  " << setw(8) << l
             << "  " << setw(14) << pdf << "\n";
      }
      cout << "\n";
    }
  }
  return;
}
//****************************************************************************80

void test0252 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0252 tests BUFFON_LAPLACE_SIMULATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double a;
  double b;
  double err;
  int hits;
  double l;
  double pi = 3.141592653589793238462643;
  double pi_est;
  int test;
  int trial_num;
  int trial_num_test[TEST_NUM] = { 10, 100, 10000, 1000000 };

  a = 1.0;
  b = 1.0;
  l = 1.0;

  cout << "\n";
  cout << "TEST0252\n";
  cout << "  BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle dropping\n";
  cout << "  experiment.  On a grid of cells of width A and height B,\n";
  cout << "  a needle of length L is dropped at random.  We count\n";
  cout << "  the number of times it crosses at least one grid line,\n";
  cout << "  and use this to estimate the value of PI.\n";

  cout << "\n";
  cout << "  Cell width A =    " << a << "\n";
  cout << "  Cell height B =   " << b << "\n";
  cout << "  Needle length L = " << l << "\n";
  cout << "\n";
  cout << "    Trials      Hits          Est(Pi)     Err\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    trial_num = trial_num_test[test];

    hits = buffon_laplace_simulate ( a, b, l, trial_num );

    if ( 0 < hits )
    {
      pi_est = ( 2.0 * l * ( a + b ) - l * l ) * ( double ) trial_num
        / ( a * b * ( double ) hits );
    }
    else
    {
      pi_est = r8_huge ( );
    }

    err = r8_abs ( pi_est - pi );

    cout << "  " << setw(8) << trial_num
         << "  " << setw(8) << hits
         << "  " << setw(14) << pi_est
         << "  " << setw(14) << err << "\n";
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test0253 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0253 tests BUFFON_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  int i;
  int j;
  int k;
  double l;
  double pdf;

  cout << "\n";
  cout << "TEST0253\n";
  cout << "  BUFFON_PDF evaluates the Buffon PDF, the probability\n";
  cout << "  that, on a grid of cells of width A,\n";
  cout << "  a needle of length L, dropped at random, will cross\n";
  cout << "  at least one grid line.\n";
  cout << "\n";
  cout << "      A         L        PDF\n";
  cout << "\n";

  for ( i = 1; i <= 5; i++ )
  {
    a = ( double ) ( i );
    for ( k = 0; k <= 5; k++ )
    {
      l = ( double ) ( k ) * a / 5.0;
      pdf = buffon_pdf ( a, l );
      cout << "  " << setw(8) << a
           << "  " << setw(8) << l
           << "  " << setw(14) << pdf << "\n";
    }
    cout << "\n";
  }
  return;
}
//****************************************************************************80

void test0254 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0254 tests BUFFON_SIMULATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define TEST_NUM 4

  double a;
  double err;
  int hits;
  double l;
  double pi = 3.141592653589793238462643;
  double pi_est;
  int test;
  int trial_num;
  int trial_num_test[TEST_NUM] = { 10, 100, 10000, 1000000 };

  a = 1.0;
  l = 1.0;

  cout << "\n";
  cout << "TEST0254\n";
  cout << "  BUFFON_SIMULATE simulates a Buffon needle dropping\n";
  cout << "  experiment.  On a grid of cells of width A,\n";
  cout << "  a needle of length L is dropped at random.  We count\n";
  cout << "  the number of times it crosses at least one grid line,\n";
  cout << "  and use this to estimate the value of PI.\n";

  cout << "\n";
  cout << "  Cell width A =    " << a << "\n";
  cout << "  Needle length L = " << l << "\n";
  cout << "\n";
  cout << "    Trials      Hits          Est(Pi)     Err\n";
  cout << "\n";

  for ( test = 0; test < TEST_NUM; test++ )
  {
    trial_num = trial_num_test[test];

    hits = buffon_simulate ( a, l, trial_num );

    if ( 0 < hits )
    {
      pi_est = ( 2.0 * l * ( double ) trial_num ) / ( a * ( double ) hits );
    }
    else
    {
      pi_est = r8_huge ( );
    }

    err = r8_abs ( pi_est - pi );

    cout << "  " << setw(8) << trial_num
         << "  " << setw(8) << hits
         << "  " << setw(14) << pi_est
         << "  " << setw(14) << err << "\n";
  }
  return;
# undef TEST_NUM
}
//****************************************************************************80

void test026 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST026 tests BURR_CDF, BURR_CDF_INV, BURR_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  double d;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST026\n";
  cout << "  For the Burr PDF:\n";
  cout << "  BURR_CDF evaluates the CDF;\n";
  cout << "  BURR_CDF_INV inverts the CDF.\n";
  cout << "  BURR_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;
  d = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";
  cout << "  PDF parameter D =      " << d << "\n";

  if ( !burr_check ( a, b, c, d ) )
  {
    cout << "\n";
    cout << "TEST026 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = burr_sample ( a, b, c, d, &seed );
    pdf = burr_pdf ( x, a, b, c, d );
    cdf = burr_cdf ( x, a, b, c, d );
    x2 = burr_cdf_inv ( cdf, a, b, c, d );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test027 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST027 tests BURR_MEAN, BURR_SAMPLE, BURR_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  double d;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST027\n";
  cout << "  For the Burr PDF:\n";
  cout << "  BURR_MEAN computes the mean;\n";
  cout << "  BURR_SAMPLE samples;\n";
  cout << "  BURR_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;
  d = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";
  cout << "  PDF parameter D =      " << d << "\n";

  if ( !burr_check ( a, b, c, d ) )
  {
    cout << "\n";
    cout << "TEST027 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = burr_mean ( a, b, c, d );
  variance = burr_variance ( a, b, c, d );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = burr_sample ( a, b, c, d, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0275 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0275 tests CARDIOID_CDF, CARDIOID_CDF_INV, CARDIOID_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2005
//
//  Author:
//
//    John Burkardt
//
{
  double a = 0.0;
  double b = 0.25;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST0275\n";
  cout << "  For the Cardioid PDF:\n";
  cout << "  CARDIOID_CDF evaluates the CDF;\n";
  cout << "  CARDIOID_CDF_INV inverts the CDF.\n";
  cout << "  CARDIOID_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "  PDF parameter A = " << a << "\n";
  cout << "  PDF parameter B = " << b << "\n";

  if ( !cardioid_check ( a, b ) )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 0; i < 10; i++ )
  {
    x = cardioid_sample ( a, b, &seed );
    pdf = cardioid_pdf ( x, a, b );
    cdf = cardioid_cdf ( x, a, b );
    x2 = cardioid_cdf_inv ( cdf, a, b );
    cout << "  " << setw(12) << x
         << "  " << setw(12) << pdf
         << "  " << setw(12) << cdf
         << "  " << setw(12) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0276 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0276 tests CARDIOID_MEAN, CARDIOID_SAMPLE, CARDIOID_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a = 0.0;
  double b = 0.25;
  int i;
  int imax;
  int imin;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST0276\n";
  cout << "  For the Cardioid PDF:\n";
  cout << "  CARDIOID_MEAN computes the mean;\n";
  cout << "  CARDIOID_SAMPLE samples;\n";
  cout << "  CARDIOID_VARIANCE computes the variance.\n";

  cout << "\n";
  cout << "  PDF parameter A = " << a << "\n";
  cout << "  PDF parameter B = " << b << "\n";

  if ( !cardioid_check ( a, b ) )
  {
    cout << "\n";
    cout << "Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = cardioid_mean ( a, b );
  variance = cardioid_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =                    " << mean << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cardioid_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test028 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST028 tests CAUCHY_CDF, CAUCHY_CDF_INV, CAUCHY_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST028\n";
  cout << "  For the Cauchy PDF:\n";
  cout << "  CAUCHY_CDF evaluates the CDF;\n";
  cout << "  CAUCHY_CDF_INV inverts the CDF.\n";
  cout << "  CAUCHY_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !cauchy_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST028 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = cauchy_sample ( a, b, &seed );
    pdf = cauchy_pdf ( x, a, b );
    cdf = cauchy_cdf ( x, a, b );
    x2 = cauchy_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test029 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST029 tests CAUCHY_MEAN, CAUCHY_SAMPLE, CAUCHY_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST029\n";
  cout << "  For the Cauchy PDF:\n";
  cout << "  CAUCHY_MEAN computes the mean;\n";
  cout << "  CAUCHY_SAMPLE samples;\n";
  cout << "  CAUCHY_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !cauchy_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST029 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = cauchy_mean ( a, b );
  variance = cauchy_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cauchy_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test030 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST030 tests CHI_CDF, CHI_CDF_INV, CHI_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST030\n";
  cout << "  For the Chi PDF:\n";
  cout << "  CHI_CDF evaluates the CDF;\n";
  cout << "  CHI_CDF_INV inverts the CDF.\n";
  cout << "  CHI_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !chi_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST030 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = chi_sample ( a, b, c, &seed );
    pdf = chi_pdf ( x, a, b, c );
    cdf = chi_cdf ( x, a, b, c );
    x2 = chi_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test031 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST031 tests CHI_MEAN, CHI_SAMPLE, CHI_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST031\n";
  cout << "  For the Chi PDF:\n";
  cout << "  CHI_MEAN computes the mean;\n";
  cout << "  CHI_SAMPLE samples;\n";
  cout << "  CHI_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !chi_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST031 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = chi_mean ( a, b, c );
  variance = chi_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = chi_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test032 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST032 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double a2;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST032:\n";
  cout << "  CHI_SQUARE_CDF evaluates the cumulative\n";
  cout << "    distribution function for the chi-square central\n";
  cout << "    probability density function.\n";
  cout << "  CHI_SQUARE_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      A     X   Exact F     CHI_SQUARE_CDF(A,X)\n";

  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    chi_square_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    a2 = ( double ) a;

    fx2 = chi_square_cdf ( x, a2 );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test033 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST033 tests CHI_SQUARE_CDF, CHI_SQUARE_CDF_INV, CHI_SQUARE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST033\n";
  cout << "  For the Chi Square PDF:\n";
  cout << "  CHI_SQUARE_CDF evaluates the CDF;\n";
  cout << "  CHI_SQUARE_CDF_INV inverts the CDF.\n";
  cout << "  CHI_SQUARE_PDF evaluates the PDF;\n";

  a = 4.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !chi_square_check ( a ) )
  {
    cout << "\n";
    cout << "TEST033 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = chi_square_sample ( a, &seed );
    pdf = chi_square_pdf ( x, a );
    cdf = chi_square_cdf ( x, a );
    x2 = chi_square_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test034 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST034 tests CHI_SQUARE_MEAN, CHI_SQUARE_SAMPLE, CHI_SQUARE_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST034\n";
  cout << "  For the Chi Square PDF:\n";
  cout << "  CHI_SQUARE_MEAN computes the mean;\n";
  cout << "  CHI_SQUARE_SAMPLE samples;\n";
  cout << "  CHI_SQUARE_VARIANCE computes the variance.\n";

  a = 10.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !chi_square_check ( a ) )
  {
    cout << "\n";
    cout << "TEST034 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = chi_square_mean ( a );
  variance = chi_square_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = chi_square_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test035 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST035 tests CHI_SQUARE_NONCENTRAL_MEAN, CHI_SQUARE_NONCENTRAL_SAMPLE, CHI_SQUARE_NONCENTRAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  int j;
  double mean;
  int seed;
  int seed_init = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST035\n";
  cout << "  For the Chi Square Noncentral PDF:\n";
  cout << "  CHI_SQUARE_NONCENTRAL_MEAN computes the mean;\n";
  cout << "  CHI_SQUARE_NONCENTRAL_SAMPLE samples;\n";
  cout << "  CHI_SQUARE_NONCENTRAL_VARIANCE computes the variance;\n";

  a = 3.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !chi_square_noncentral_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST035 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = chi_square_noncentral_mean ( a, b );
  variance = chi_square_noncentral_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  seed = seed_init;

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = chi_square_noncentral_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Initial seed =     " << seed_init  << "\n";
  cout << "  Final seed =       " << seed       << "\n";

  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test036 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST036 tests CIRCLE_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  cout << "\n";
  cout << "TEST036\n";
  cout << "  For the Circle PDF:\n";
  cout << "  CIRCLE_SAMPLE samples;\n";

  a = 10.0;
  b = 4.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circle_sample ( a, b, c, &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    delete [] y;
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     "
    << setw(12) << mean[0] << "  "
    << setw(12) << mean[1] << "\n";
  cout << "  Sample variance = "
    << setw(12) << variance[0] << "  "
    << setw(12) << variance[1] << "\n";
  cout << "  Sample maximum =  "
    << setw(12) << xmax[0] << "  "
    << setw(12) << xmax[1] << "\n";
  cout << "  Sample minimum =  "
    << setw(12) << xmin[0] << "  "
    << setw(12) << xmin[1] << "\n";

  delete [] mean;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test037 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST037 tests CIRCULAR_NORMAL_01_*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  cout << "\n";
  cout << "TEST037\n";
  cout << "  For the Circular Normal 01 PDF:\n";
  cout << "  CIRCULAR_NORMAL_01_MEAN computes the mean;\n";
  cout << "  CIRCULAR_NORMAL_01_SAMPLE samples;\n";
  cout << "  CIRCULAR_NORMAL_01_VARIANCE computes the variance.\n";

  mean = circular_normal_01_mean ( );
  variance = circular_normal_01_variance ( );

  cout << "\n";
  cout << "  PDF mean =     "
    << setw(12) << mean[0] << "  "
    << setw(12) << mean[1] << "\n";
  cout << "  PDF variance = "
    << setw(12) << variance[0] << "  "
    << setw(12) << variance[1] << "\n";

  delete [] mean;
  delete [] variance;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circular_normal_01_sample ( &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    delete [] y;
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     "
    << setw(12) << mean[0] << "  "
    << setw(12) << mean[1] << "\n";
  cout << "  Sample variance = "
    << setw(12) << variance[0] << "  "
    << setw(12) << variance[1] << "\n";
  cout << "  Sample maximum =  "
    << setw(12) << xmax[0] << "  "
    << setw(12) << xmax[1] << "\n";
  cout << "  Sample minimum =  "
    << setw(12) << xmin[0] << "  "
    << setw(12) << xmin[1] << "\n";

  delete [] mean;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0375 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0375 tests CIRCULAR_NORMAL_*.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 January 2012
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a[2];
  double b;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  double x[2*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  a[0] = 1.0;
  a[1] = 5.0;
  b = 0.75;

  cout << "\n";
  cout << "TEST0375\n";
  cout << "  For the Circular Normal PDF:\n";
  cout << "  CIRCULAR_NORMAL_MEAN computes the mean;\n";
  cout << "  CIRCULAR_NORMAL_SAMPLE samples;\n";
  cout << "  CIRCULAR_NORMAL_VARIANCE computes the variance.\n";

  mean = circular_normal_mean ( a, b );
  variance = circular_normal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     "
    << setw(12) << mean[0] << "  "
    << setw(12) << mean[1] << "\n";
  cout << "  PDF variance = "
    << setw(12) << variance[0] << "  "
    << setw(12) << variance[1] << "\n";

  delete [] mean;
  delete [] variance;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = circular_normal_sample ( a, b, &seed );
    x[0+j*2] = y[0];
    x[1+j*2] = y[1];
    delete [] y;
  }

  mean = r8row_mean ( 2, SAMPLE_NUM, x );
  variance = r8row_variance ( 2, SAMPLE_NUM, x );
  xmax = r8row_max ( 2, SAMPLE_NUM, x );
  xmin = r8row_min ( 2, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     "
    << setw(12) << mean[0] << "  "
    << setw(12) << mean[1] << "\n";
  cout << "  Sample variance = "
    << setw(12) << variance[0] << "  "
    << setw(12) << variance[1] << "\n";
  cout << "  Sample maximum =  "
    << setw(12) << xmax[0] << "  "
    << setw(12) << xmax[1] << "\n";
  cout << "  Sample minimum =  "
    << setw(12) << xmin[0] << "  "
    << setw(12) << xmin[1] << "\n";

  delete [] mean;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test038 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST038 tests COSINE_CDF, COSINE_CDF_INV, COSINE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST038\n";
  cout << "  For the Cosine PDF:\n";
  cout << "  COSINE_CDF evaluates the CDF;\n";
  cout << "  COSINE_CDF_INV inverts the CDF.\n";
  cout << "  COSINE_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 1.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !cosine_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST038 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = cosine_sample ( a, b, &seed );
    pdf = cosine_pdf ( x, a, b );
    cdf = cosine_cdf ( x, a, b );
    x2 = cosine_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test039 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST039 tests COSINE_MEAN, COSINE_SAMPLE, COSINE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST039\n";
  cout << "  For the Cosine PDF:\n";
  cout << "  COSINE_MEAN computes the mean;\n";
  cout << "  COSINE_SAMPLE samples;\n";
  cout << "  COSINE_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 1.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !cosine_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST039 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = cosine_mean ( a, b );
  variance = cosine_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = cosine_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0395 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0395 tests COUPON_COMPLETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  int box_num;
  double cdf;
  double pdf;
  int type_num;

  cout << "\n";
  cout << "TEST0395\n";
  cout << "  COUPON_COMPLETE_PDF evaluates the coupon collector's\n";
  cout << "  complete collection pdf.\n";
  cout << "\n";

  for ( type_num = 2; type_num <= 4; type_num++ )
  {
    cout << "\n";
    cout << "  Number of coupon types is " << type_num << "\n";
    cout << "\n";
    cout << "   BOX_NUM      PDF             CDF\n";
    cout << "\n";
    cdf = 0.0;
    for ( box_num = 1; box_num <= 20; box_num++ )
    {
      pdf = coupon_complete_pdf ( type_num, box_num );
      cdf = cdf + pdf;
      cout << "  " << setw(8) << box_num
           << "  " << setw(14) << pdf
           << "  " << setw(14) << cdf << "\n";
    }
  }

  return;
}
//****************************************************************************80

void test040 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST040 tests COUPON_SIMULATE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N_TRIAL 10
# define MAX_TYPE 25

  double average;
  int coupon[MAX_TYPE];
  double expect;
  int i;
  int n_coupon;
  int n_type;
  int seed = 123456789;

  cout << "\n";
  cout << "TEST040:\n";
  cout << "  COUPON_SIMULATE simulates the couponn";
  cout << "  collector's problem.\n";
  cout << "\n";

  for ( n_type = 5; n_type <= MAX_TYPE; n_type = n_type + 5 )
  {
    cout << "\n";
    cout << "  Number of coupon types is " << n_type << "\n";
    expect = ( double ) ( n_type ) * log ( ( double ) ( n_type ) );
    cout << "  Expected wait is about " << expect << "\n";
    cout << "\n";

    average = 0.0;
    for ( i = 1; i <= N_TRIAL; i++ )
    {
      coupon_simulate ( n_type, &seed, coupon, &n_coupon );

      cout << "  "
           << setw(6) << i        << "  "
           << setw(6) << n_coupon << "\n";

      average = average + ( double ) ( n_coupon );
    }

    average = average / ( double ) ( N_TRIAL );
    cout << "\n";
    cout << "  Average wait was " << average << "\n";
  }

  return;
# undef N_TRIAL
# undef MAX_TRIAL
}
//****************************************************************************80

void test041 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST041 tests DERANGED_CDF, DERANGED_CDF_INV, DERANGED_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST041\n";
  cout << "  For the Deranged PDF:\n";
  cout << "  DERANGED_CDF evaluates the CDF;\n";
  cout << "  DERANGED_CDF_INV inverts the CDF.\n";
  cout << "  DERANGED_PDF evaluates the PDF;\n";

  a = 7;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !deranged_check ( a ) )
  {
    cout << "\n";
    cout << "TEST041 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = deranged_sample ( a, &seed );
    pdf = deranged_pdf ( x, a );
    cdf = deranged_cdf ( x, a );
    x2 = deranged_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test042 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST042 tests DERANGED_CDF, DERANGED_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST042\n";
  cout << "  For the Deranged PDF:\n";
  cout << "  DERANGED_CDF evaluates the CDF;\n";
  cout << "  DERANGED_PDF evaluates the PDF;\n";

  a = 7;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !deranged_check ( a ) )
  {
    cout << "\n";
    cout << "TEST042 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 0; x <= a; x++ )
  {
    pdf = deranged_pdf ( x, a );
    cdf = deranged_cdf ( x, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test043 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST043 tests DERANGED_MEAN, DERANGED_SAMPLE, DERANGED_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST043\n";
  cout << "  For the Deranged PDF:\n";
  cout << "  DERANGED_MEAN computes the mean;\n";
  cout << "  DERANGED_SAMPLE samples;\n";
  cout << "  DERANGED_VARIANCE computes the variance.\n";

  a = 7;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !deranged_check ( a ) )
  {
    cout << "\n";
    cout << "TEST042 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = deranged_mean ( a );
  variance = deranged_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = deranged_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test044 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST044 tests DIGAMMA, PSI_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST044:\n";
  cout << "  DIGAMMA evaluates the DIGAMMA or PSI function.\n";
  cout << "  PSI_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      X   Exact F     DIGAMMA(X)\n";

  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    psi_values ( &n_data, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    if ( x <= 0.0 )
    {
      continue;
    }
    fx2 = digamma ( x );

    cout << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test045 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST045 tests DIPOLE_CDF, DIPOLE_CDF_INV, DIPOLE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define PI 3.14159265358979323
# define TEST_NUM 3

  double a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int test_i;
  double test_a[TEST_NUM] = { 0.0, PI/4.0, PI/2.0 };
  double test_b[TEST_NUM] = { 1.0, 0.5,    0.0 };
  double x;
  double x2;

  cout << "\n";
  cout << "TEST045\n";
  cout << "  For the Cosine PDF:\n";
  cout << "  DIPOLE_CDF evaluates the CDF;\n";
  cout << "  DIPOLE_CDF_INV inverts the CDF.\n";
  cout << "  DIPOLE_PDF evaluates the PDF;\n";

  for ( test_i = 0; test_i < TEST_NUM; test_i++ )
  {
    a = test_a[test_i];
    b = test_b[test_i];

    cout << "\n";
    cout << "  PDF parameter A =      " << a << "\n";
    cout << "  PDF parameter B =      " << b << "\n";

    if ( !dipole_check ( a, b ) )
    {
      cout << "\n";
      cout << "TEST045 - Fatal error!\n";
      cout << "  The parameters are not legal.\n";
      return;
    }

    cout << "\n";
    cout << "       X            PDF           CDF            CDF_INV\n";
    cout << "\n";

    for ( i = 1; i <= 10; i++ )
    {
      x = dipole_sample ( a, b, &seed );
      pdf = dipole_pdf ( x, a, b );
      cdf = dipole_cdf ( x, a, b );
      x2 = dipole_cdf_inv ( cdf, a, b );

      cout << "  "
           << setw(12) << x   << "  "
           << setw(12) << pdf << "  "
           << setw(12) << cdf << "  "
           << setw(12) << x2  << "\n";
    }

  }

  return;
# undef PI
# undef TEST_NUM
}
//****************************************************************************80

void test046 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST046 tests DIPOLE_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000
# define PI 3.14159265358979323
# define TEST_NUM 3

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double test_a[TEST_NUM] = { 0.0, PI/4.0, PI/2.0 };
  double test_b[TEST_NUM] = { 1.0, 0.5,    0.0 };
  int test_i;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST046\n";
  cout << "  For the Cosine PDF:\n";
  cout << "  DIPOLE_SAMPLE samples;\n";

  for ( test_i = 0; test_i < TEST_NUM; test_i++ )
  {
    a = test_a[test_i];
    b = test_b[test_i];

    cout << "\n";
    cout << "  PDF parameter A =      " << a << "\n";
    cout << "  PDF parameter B =      " << b << "\n";

    if ( !dipole_check ( a, b ) )
    {
      cout << "\n";
      cout << "TEST046 - Fatal error!\n";
      cout << "  The parameters are not legal.\n";
      return;
    }

    for ( i = 0; i < SAMPLE_NUM; i++ )
    {
      x[i] = dipole_sample ( a, b, &seed );
    }

    mean = r8vec_mean ( SAMPLE_NUM, x );
    variance = r8vec_variance ( SAMPLE_NUM, x );
    xmax = r8vec_max ( SAMPLE_NUM, x );
    xmin = r8vec_min ( SAMPLE_NUM, x );

    cout << "\n";
    cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
    cout << "  Sample mean =     " << mean     << "\n";
    cout << "  Sample variance = " << variance << "\n";
    cout << "  Sample maximum =  " << xmax     << "\n";
    cout << "  Sample minimum =  " << xmin     << "\n";
  }

  return;
# undef SAMPLE_NUM
# undef PI
# undef TEST_NUM
}
//****************************************************************************80

void test047 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST047 tests DIRICHLET_MEAN, DIRICHLET_SAMPLE, DIRICHLET_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3
# define SAMPLE_NUM 1000

  double a[N] = { 0.250, 0.500, 1.250 };
  int i;
  int j;
  double *mean;
  double *m2;
  int seed = 123456789;
  double *variance;
  double x[N*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  cout << "\n";
  cout << "TEST047\n";
  cout << "  For the Dirichlet PDF:\n";
  cout << "  DIRICHLET_MEAN computes the mean;\n";
  cout << "  DIRICHLET_SAMPLE samples;\n";
  cout << "  DIRICHLET_VARIANCE computes the variance;\n";

  cout << "\n";
  cout << "  Number of components N = " << N << "\n";

  r8vec_print ( N, a, "  PDF parameter A:" );

  if ( !dirichlet_check ( N, a ) )
  {
    cout << "\n";
    cout << "TEST047 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = dirichlet_mean ( N, a );
  variance = dirichlet_variance ( N, a );
  r8vec_print ( N, mean, "  PDF mean:" );
  r8vec_print ( N, variance, "  PDF variance:" );

  delete [] mean;
  delete [] variance;

  m2 = dirichlet_moment2 ( N, a );

  r8mat_print ( N, N, m2, "  Second moment matrix:" );

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = dirichlet_sample ( N, a, &seed );
    for ( i = 0; i < N; i++ )
    {
      x[i+j*N] = y[i];
    }
    delete [] y;
  }

  mean = r8row_mean ( N, SAMPLE_NUM, x );
  variance = r8row_variance ( N, SAMPLE_NUM, x );
  xmax = r8row_max ( N, SAMPLE_NUM, x );
  xmin = r8row_min ( N, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "\n";
  cout << "  Component Mean, Variance, Min, Max:\n";
  cout << "\n";

  for ( i = 0; i < N; i++ )
  {
    cout << "  "
         << setw(6)  << i           << "  "
         << setw(12) << mean[i]     << "  "
         << setw(12) << variance[i] << "  "
         << setw(12) << xmax[i]     << "  "
         << setw(12) << xmin[i]     << "\n";
  }

  delete [] mean;
  delete [] m2;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
# undef N
# undef SAMPLE_NUM
}
//****************************************************************************80

void test048 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST048 tests DIRICHLET_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double a[N] = { 0.250, 0.500, 1.250 };
  int i;
  double pdf;
  double x[N] = { 0.500, 0.125, 0.375 };

  cout << "\n";
  cout << "TEST048\n";
  cout << "  For the Dirichlet PDF:\n";
  cout << "  DIRICHLET_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "  Number of components N = " << N << "\n";

  r8vec_print ( N, a, "  PDF parameter A:" );

  if ( !dirichlet_check ( N, a ) )
  {
    cout << "\n";
    cout << "TEST048 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  r8vec_print ( N, x, "  PDF argument X:" );

  pdf = dirichlet_pdf ( x, N, a );

  cout << "\n";
  cout << "  PDF value = " << pdf  << "\n";

  return;
# undef N
}
//****************************************************************************80

void test049 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST049 tests DIRICHLET_MIX_MEAN, DIRICHLET_MIX_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define COMP_NUM 2
# define ELEM_NUM 3
# define SAMPLE_NUM 1000

  double a[ELEM_NUM*COMP_NUM] = {
    0.250, 0.500, 1.250,
    1.500, 0.500, 2.000 };
  int comp;
  int comp_i;
  double comp_weight[COMP_NUM] = { 1.0, 2.0 };
  int elem_i;
  int i;
  int j;
  double *mean;
  double pdf;
  int seed = 123456789;
  double *variance;
  double x[ELEM_NUM*SAMPLE_NUM];
  double *xmax;
  double *xmin;
  double *y;

  cout << "\n";
  cout << "TEST049\n";
  cout << "  For the Dirichlet Mixture PDF:\n";
  cout << "  DIRICHLET_MIX_SAMPLE samples;\n";
  cout << "  DIRICHLET_MIX_MEAN computes the mean;\n";

  cout << "\n";
  cout << "  Number of elements ELEM_NUM =   " << ELEM_NUM << "\n";
  cout << "  Number of components COMP_NUM = " << COMP_NUM << "\n";
  r8mat_print ( ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):" );
  r8vec_print ( COMP_NUM, comp_weight, "  Component weights" );

  if ( !dirichlet_mix_check ( COMP_NUM, ELEM_NUM, a, comp_weight ) )
  {
    cout << "\n";
    cout << "TEST049 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = dirichlet_mix_mean ( COMP_NUM, ELEM_NUM, a, comp_weight );

  r8vec_print ( ELEM_NUM, mean, "  PDF mean" );

  delete [] mean;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = dirichlet_mix_sample ( COMP_NUM, ELEM_NUM, a, comp_weight, &seed,
      &comp );

    for ( i = 0; i < ELEM_NUM; i++ )
    {
      x[i+j*ELEM_NUM] = y[i];
    }
    delete [] y;
  }

  mean = r8row_mean ( ELEM_NUM, SAMPLE_NUM, x );
  variance = r8row_variance ( ELEM_NUM, SAMPLE_NUM, x );
  xmax = r8row_max ( ELEM_NUM, SAMPLE_NUM, x );
  xmin = r8row_min ( ELEM_NUM, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "\n";
  cout << "  Component Mean, Variance, Max, Min:\n";
  cout << "\n";

  for ( i = 0; i < ELEM_NUM; i++ )
  {
    cout << "  "
         << setw(6)  << i           << "  "
         << setw(12) << mean[i]     << "  "
         << setw(12) << variance[i] << "  "
         << setw(12) << xmax[i]     << "  "
         << setw(12) << xmin[i]     << "\n";
  }

  delete [] mean;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
}
//****************************************************************************80

void test050 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST050 tests DIRICHLET_MIX_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define COMP_NUM 2
# define ELEM_NUM 3

  double a[ELEM_NUM*COMP_NUM] = {
    0.250, 0.500, 1.250,
    1.500, 0.500, 2.000 };
  int comp_i;
  double comp_weight[COMP_NUM] = { 1.0, 2.0 };
  int elem_i;
  double pdf;
  double x[ELEM_NUM] = { 0.500, 0.125, 0.375 };

  cout << "\n";
  cout << "TEST050\n";
  cout << "  For the Dirichlet mixture PDF:\n";
  cout << "  DIRICHLET_MIX_PDF evaluates the PDF.\n";

  cout << flush;

  cout << "\n";
  cout << "  Number of elements ELEM_NUM =   " << ELEM_NUM << "\n";
  cout << "  Number of components COMP_NUM = " << COMP_NUM << "\n";
  r8mat_print ( ELEM_NUM, COMP_NUM, a, "  PDF parameters A(ELEM,COMP):" );
  r8vec_print ( COMP_NUM, comp_weight, "  Component weights" );

  cout << flush;

  if ( !dirichlet_mix_check ( COMP_NUM, ELEM_NUM, a, comp_weight ) )
  {
    cout << "\n";
    cout << "TEST050 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  r8vec_print ( ELEM_NUM, x, "  PDF argument X:" );

  cout << flush;

  pdf = dirichlet_mix_pdf ( x, COMP_NUM, ELEM_NUM, a, comp_weight );

  cout << "\n";
  cout << "  PDF value =           " << pdf << "\n";

  return;
# undef COMP_NUM
# undef ELEM_NUM
}
//****************************************************************************80

void test051 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST051 tests BETA_PDF, DIRICHLET_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2

  double a;
  double aval;
  double avec[N];
  double b;
  double bval;
  int i;
  double pdf;
  double x;
  double xval;
  double xvec[N];

  xval = 0.25;
  aval = 2.50;
  bval = 3.50;

  cout << "\n";
  cout << "TEST051\n";
  cout << "  BETA_PDF evaluates the Beta PDF.\n";
  cout << "  DIRICHLET_PDF evaluates the Dirichlet PDF.\n";
  cout << "\n";
  cout << "  For N = 2, Dirichlet = Beta.\n";

  avec[0] = aval;
  avec[1] = bval;

  cout << "\n";
  cout << "  Number of components N = " << N << "\n";
  r8vec_print ( N, avec, "  PDF parameters A(1:N):" );

  if ( !dirichlet_check ( N, avec ) )
  {
    cout << "\n";
    cout << "TEST051 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  xvec[0] = xval;
  xvec[1] = 1.0 - xval;

  r8vec_print ( N, xvec, "  PDF arguments X(1:N):" );

  pdf = dirichlet_pdf ( xvec, N, avec );

  cout << "\n";
  cout << "  Dirichlet PDF value =  " << pdf << "\n";

  x = xval;
  a = aval;
  b = bval;

  pdf = beta_pdf ( x, a, b );

  cout << "  Beta PDF value =       " << pdf << "\n";

  return;
# undef N
}
//****************************************************************************80

void test052 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST052 tests DISCRETE_CDF, DISCRETE_CDF_INV, DISCRETE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define A 6

  double b[A] = { 1.0, 2.0, 6.0, 2.0, 4.0, 1.0 };
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST052\n";
  cout << "  For the Discrete PDF:\n";
  cout << "  DISCRETE_CDF evaluates the CDF;\n";
  cout << "  DISCRETE_CDF_INV inverts the CDF.\n";
  cout << "  DISCRETE_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "  PDF parameter A =      " << A << "\n";
  r8vec_print ( A, b, "  PDF parameter B:" );

  if ( !discrete_check ( A, b ) )
  {
    cout << "\n";
    cout << "TEST052 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }


  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = discrete_sample ( A, b, &seed );
    pdf = discrete_pdf ( x, A, b );
    cdf = discrete_cdf ( x, A, b );
    x2 = discrete_cdf_inv ( cdf, A, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
# undef A
}
//****************************************************************************80

void test053 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST053 tests DISCRETE_MEAN, DISCRETE_SAMPLE, DISCRETE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define A 6
# define SAMPLE_NUM 1000

  double b[A] = { 1.0, 2.0, 6.0, 2.0, 4.0, 1.0 };
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST053\n";
  cout << "  For the Discrete PDF:\n";
  cout << "  DISCRETE_MEAN computes the mean;\n";
  cout << "  DISCRETE_SAMPLE samples;\n";
  cout << "  DISCRETE_VARIANCE computes the variance;\n";

  cout << "\n";
  cout << "  PDF parameter A =      " << A << "\n";
  r8vec_print ( A, b, "  PDF parameter B:" );

  if ( !discrete_check ( A, b ) )
  {
    cout << "\n";
    cout << "TEST053 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = discrete_mean ( A, b );
  variance = discrete_variance ( A, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = discrete_sample ( A, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef A
# undef SAMPLE_NUM
}
//****************************************************************************80

void test054 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST054 tests EMPIRICAL_DISCRETE_CDF, EMPIRICAL_DISCRETE_CDF_INV, EMPIRICAL_DISCRETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define A 6

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST054\n";
  cout << "  For the Empirical Discrete PDF:\n";
  cout << "  EMPIRICAL_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  EMPIRICAL_DISCRETE_CDF_INV inverts the CDF.\n";
  cout << "  EMPIRICAL_DISCRETE_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "  PDF parameter A = " << A   << "\n";
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    cout << "\n";
    cout << "TEST054 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = empirical_discrete_sample ( A, b, c, &seed );
    pdf = empirical_discrete_pdf ( x, A, b, c );
    cdf = empirical_discrete_cdf ( x, A, b, c );
    x2 = empirical_discrete_cdf_inv ( cdf, A, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
# undef A
}
//****************************************************************************80

void test055 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST055 tests EMPIRICAL_DISCRETE_MEAN, EMPIRICAL_DISCRETE_SAMPLE, EMPIRICAL_DISCRETE_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2008
//
//  Author:
//
//    John Burkardt
//
{
# define A 6
# define SAMPLE_NUM 1000

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST055\n";
  cout << "  For the Empirical Discrete PDF:\n";
  cout << "  EMPIRICAL_DISCRETE_MEAN computes the mean;\n";
  cout << "  EMPIRICAL_DISCRETE_SAMPLE samples;\n";
  cout << "  EMPIRICAL_DISCRETE_VARIANCE computes the variance.\n";

  cout << "\n";
  cout << "  PDF parameter A = " << A   << "\n";
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    cout << "\n";
    cout << "TEST055 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = empirical_discrete_mean ( A, b, c );
  variance = empirical_discrete_variance ( A, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = empirical_discrete_sample ( A, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef A
# undef SAMPLE_NUM
}
//****************************************************************************80

void test056 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST056 tests EMPIRICAL_DISCRETE_CDF, EMPIRICAL_DISCRETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define A 6

  double b[A] = { 1.0, 1.0, 3.0, 2.0, 1.0, 2.0 };
  double c[A] = { 0.0, 1.0, 2.0, 4.5, 6.0, 10.0 };
  double cdf;
  int i;
  double pdf;
  double x;

  cout << "\n";
  cout << "TEST056\n";
  cout << "  For the Empirical Discrete PDF:\n";
  cout << "  EMPIRICAL_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  EMPIRICAL_DISCRETE_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "  PDF parameter A = " << A   << "\n";
  r8vec_print ( A, b, "  PDF parameter B = " );
  r8vec_print ( A, c, "  PDF parameter C = " );

  if ( !empirical_discrete_check ( A, b, c ) )
  {
    cout << "\n";
    cout << "TEST056 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( i = -2; i <= 12; i++ )
  {
    x = ( double ) i;
    pdf = empirical_discrete_pdf ( x, A, b, c );
    cdf = empirical_discrete_cdf ( x, A, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
# undef A
}
//****************************************************************************80

void test0563 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0563 tests ENGLISH_SENTENCE_LENGTH_CDF, ENGLISH_SENTENCE_LENGTH_CDF_INV and ENGLISH_SENTENCE_LENGTH_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST0563\n";
  cout << "  For the English Sentence Length PDF:\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_CDF evaluates the CDF;\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_CDF_INV inverts the CDF.\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = english_sentence_length_sample ( &seed );

    pdf = english_sentence_length_pdf ( x );

    cdf = english_sentence_length_cdf ( x );

    x2 = english_sentence_length_cdf_inv ( cdf );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << pdf
         << "  " << setw(12) << cdf
         << "  " << setw(12) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0564 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0564 tests ENGLISH_SENTENCE_LENGTH_MEAN, ENGLISH_SENTENCE_LENGTH_SAMPLE and ENGLISH_SENTENCE_LENGTH_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 August 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST0564\n";
  cout << "  For the English Sentence Length PDF:\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_MEAN computes the mean;\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_SAMPLE samples;\n";
  cout << "  ENGLISH_SENTENCE_LENGTH_VARIANCE computes the variance.\n";

  mean = english_sentence_length_mean ( );
  variance = english_sentence_length_variance ( );

  cout << "\n";
  cout << "  PDF mean =                    " << mean << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = english_sentence_length_sample ( &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0565 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0565 tests ENGLISH_WORD_LENGTH_CDF, ENGLISH_WORD_LENGTH_CDF_INV and ENGLISH_WORD_LENGTH_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST0565\n";
  cout << "  For the English Word Length PDF:\n";
  cout << "  ENGLISH_WORD_LENGTH_CDF evaluates the CDF;\n";
  cout << "  ENGLISH_WORD_LENGTH_CDF_INV inverts the CDF.\n";
  cout << "  ENGLISH_WORD_LENGTH_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = english_word_length_sample ( &seed );

    pdf = english_word_length_pdf ( x );

    cdf = english_word_length_cdf ( x );

    x2 = english_word_length_cdf_inv ( cdf );

    cout << "  " << setw(12) << x
         << "  " << setw(12) << pdf
         << "  " << setw(12) << cdf
         << "  " << setw(12) << x2 << "\n";
  }

  return;
}
//****************************************************************************80

void test0566 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0566 tests ENGLISH_WORD_LENGTH_MEAN, ENGLISH_WORD_LENGTH_SAMPLE and ENGLISH_WORD_LENGTH_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 July 2006
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST0566\n";
  cout << "  For the English Word Length PDF:\n";
  cout << "  ENGLISH_WORD_LENGTH_MEAN computes the mean;\n";
  cout << "  ENGLISH_WORD_LENGTH_SAMPLE samples;\n";
  cout << "  ENGLISH_WORD_LENGTH_VARIANCE computes the variance.\n";

  mean = english_word_length_mean ( );
  variance = english_word_length_variance ( );

  cout << "\n";
  cout << "  PDF mean =                    " << mean << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = english_word_length_sample ( &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM << "\n";
  cout << "  Sample mean =     " << mean << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax << "\n";
  cout << "  Sample minimum =  " << xmin << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test057 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST057 tests ERLANG_CDF, ERLANG_CDF_INV, ERLANG_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST057\n";
  cout << "  For the Erlang PDF:\n";
  cout << "  ERLANG_CDF evaluates the CDF;\n";
  cout << "  ERLANG_CDF_INV inverts the CDF.\n";
  cout << "  ERLANG_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !erlang_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST057 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = erlang_sample ( a, b, c, &seed );
    pdf = erlang_pdf ( x, a, b, c );
    cdf = erlang_cdf ( x, a, b, c );
    x2 = erlang_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test058 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST058 tests ERLANG_MEAN, ERLANG_SAMPLE, ERLANG_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST058\n";
  cout << "  For the Erlang PDF:\n";
  cout << "  ERLANG_MEAN computes the mean;\n";
  cout << "  ERLANG_SAMPLE samples;\n";
  cout << "  ERLANG_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !erlang_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST058 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = erlang_mean ( a, b, c );
  variance = erlang_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = erlang_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test059 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST059 tests ERROR_F and ERROR_F_INVERSE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2010
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int seed;
  double x;
  double y;
  double z;

  cout << "\n";
  cout << "TEST059\n";
  cout << "  ERROR_F evaluates ERF(X).\n";
  cout << "  ERROR_F_INVERSE inverts ERF(X).\n";
  cout << "\n";
  cout << "X   -> Y = error_F(X) -> Z = error_f_inverse(Y)\n";
  cout << "\n";

  seed = 123456789;

  x = 1.0;

  for ( i = 1; i <= 20; i++ )
  {
    x = normal_01_sample ( &seed );
    y = error_f ( x );
    z = error_f_inverse ( y );
    cout << "  " << setw(14) << x
         << "  " << setw(14) << x
         << "  " << setw(14) << z << "\n";
  }
  return;
}
//****************************************************************************80

void test060 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST060 tests EXPONENTIAL_01_CDF, EXPONENTIAL_01_CDF_INV, EXPONENTIAL_01_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST060\n";
  cout << "  For the Exponential 01 PDF:\n";
  cout << "  EXPONENTIAL_01_CDF evaluates the CDF;\n";
  cout << "  EXPONENTIAL_01_CDF_INV inverts the CDF.\n";
  cout << "  EXPONENTIAL_01_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = exponential_01_sample ( &seed );
    pdf = exponential_01_pdf ( x );
    cdf = exponential_01_cdf ( x );
    x2 = exponential_01_cdf_inv ( cdf );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test061 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST061 tests EXPONENTIAL_01_MEAN, EXPONENTIAL_01_SAMPLE, EXPONENTIAL_01_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST061\n";
  cout << "  For the Exponential 01 PDF:\n";
  cout << "  EXPONENTIAL_01_MEAN computes the mean;\n";
  cout << "  EXPONENTIAL_01_SAMPLE samples;\n";
  cout << "  EXPONENTIAL_01_VARIANCE computes the variance.\n";

  mean     = exponential_01_mean ( );
  variance = exponential_01_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = exponential_01_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test062 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST062 tests EXPONENTIAL_CDF, EXPONENTIAL_CDF_INV, EXPONENTIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST062\n";
  cout << "  For the Exponential PDF:\n";
  cout << "  EXPONENTIAL_CDF evaluates the CDF;\n";
  cout << "  EXPONENTIAL_CDF_INV inverts the CDF.\n";
  cout << "  EXPONENTIAL_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !exponential_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST062 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = exponential_sample ( a, b, &seed );
    pdf = exponential_pdf ( x, a, b );
    cdf = exponential_cdf ( x, a, b );
    x2 = exponential_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test063 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST063 tests EXPONENTIAL_MEAN, EXPONENTIAL_SAMPLE, EXPONENTIAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST063\n";
  cout << "  For the Exponential PDF:\n";
  cout << "  EXPONENTIAL_MEAN computes the mean;\n";
  cout << "  EXPONENTIAL_SAMPLE samples;\n";
  cout << "  EXPONENTIAL_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !exponential_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST063 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = exponential_mean ( a, b );
  variance = exponential_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = exponential_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test064 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST064 tests EXTREME_VALUES_CDF, EXTREME_VALUES_CDF_INV, EXTREME_VALUES_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST064\n";
  cout << "  For the Extreme Values PDF:\n";
  cout << "  EXTREME_VALUES_CDF evaluates the CDF;\n";
  cout << "  EXTREME_VALUES_CDF_INV inverts the CDF.\n";
  cout << "  EXTREME_VALUES_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !extreme_values_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST064 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = extreme_values_sample ( a, b, &seed );
    pdf = extreme_values_pdf ( x, a, b );
    cdf = extreme_values_cdf ( x, a, b );
    x2 = extreme_values_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test065 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST065 tests EXTREME_VALUES_MEAN, EXTREME_VALUES_SAMPLE, EXTREME_VALUES_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST065\n";
  cout << "  For the Extreme Values PDF:\n";
  cout << "  EXTREME_VALUES_MEAN computes the mean;\n";
  cout << "  EXTREME_VALUES_SAMPLE samples;\n";
  cout << "  EXTREME_VALUES_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !extreme_values_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST065 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = extreme_values_mean ( a, b );
  variance = extreme_values_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = extreme_values_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test066 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST066 tests F_CDF, F_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST066:\n";
  cout << "  F_CDF evaluates the cumulative\n";
  cout << "    distribution function for the F\n";
  cout << "    probability density function.\n";
  cout << "  F_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "      A     B     X   Exact F     F_CDF(A,B,X)\n";

  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    f_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }
    fx2 = f_cdf ( x, a, b );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test067 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST067 tests F_CDF, F_PDF, F_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  int m;
  int n;
  double pdf;
  int seed = 123456789;
  double x;

  cout << "\n";
  cout << "TEST067\n";
  cout << "  For the F PDF:\n";
  cout << "  F_CDF evaluates the CDF;\n";
  cout << "  F_PDF evaluates the PDF;\n";
  cout << "  F_SAMPLE samples the PDF;\n";

  m = 1;
  n = 1;

  cout << "\n";
  cout << "  PDF parameter M = " << m   << "\n";
  cout << "  PDF parameter N = " << n   << "\n";

  if ( !f_check ( m, n ) )
  {
    cout << "\n";
    cout << "TEST067 - Fatal error!\n";
    cout << "  The parameter values are illegal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = f_sample ( m, n, &seed );
    pdf = f_pdf ( x, m, n );
    cdf = f_cdf ( x, m, n );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test068 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST068 tests F_MEAN, F_SAMPLE, F_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST068\n";
  cout << "  For the F PDF:\n";
  cout << "  F_MEAN computes the mean;\n";
  cout << "  F_SAMPLE samples;\n";
  cout << "  F_VARIANCE computes the variance;\n";

  m = 8;
  n = 6;

  cout << "\n";
  cout << "  PDF parameter M = " << m << "\n";
  cout << "  PDF parameter N = " << n << "\n";

  if ( !f_check ( m, n ) )
  {
    cout << "\n";
    cout << "TEST068 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = f_mean ( m, n );
  variance = f_variance ( m, n );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = f_sample ( m, n, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test069 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST069 tests FACTORIAL_LOG, GAMMA_LOG_INT;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double f;
  double g;
  int i;
  double x;

  cout << "\n";
  cout << "TEST069\n";
  cout << "  FACTORIAL_LOG evaluates the log of the factorial function;\n";
  cout << "  GAMMA_LOG_INT evaluates the log for integer argument.\n";

  cout << "\n";
  cout << "  I, GAMMA_LOG_INT(I+1) FACTORIAL_LOG(I)\n";
  cout << "\n";

  for ( i = 1; i <= 20; i++ )
  {
    g = gamma_log_int ( i+1 );

    f = factorial_log ( i );

    cout << "  "
         << setw(6)  << i << "  "
         << setw(12) << g << "  "
         << setw(12) << f << "\n";
  }

  return;
}
//****************************************************************************80

void test070 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST070 tests FACTORIAL_STIRLING, I_FACTORIAL;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  double value;

  cout << "\n";
  cout << "TEST070\n";
  cout << "  FACTORIAL_STIRLING computes Stirling's\n";
  cout << "    approximate factorial function;\n";
  cout << "  I_FACTORIAL evaluates the factorial function;\n";
  cout << "\n";
  cout << "  N      Stirling     N!\n";
  cout << "\n";

  for ( i = 0; i <= 20; i++ )
  {
    value = factorial_stirling ( i );

    cout << "  "
         << setw(6) << i << "  "
         << setw(12) << value << "  "
         << setw(20) << i4_factorial ( i ) << "\n";
  }

  return;
}
//****************************************************************************80

void test0705 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0705 tests FISHER_PDF and FISHER_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 March 2007
//
//  Author:
//
//    John Burkardt
//
{
  int j;
  double kappa;
  double mu[3];
  int n = 10;
  double pdf;
  int seed;
  int test;
  int test_num = 3;
  double *x;

  cout << "\n";
  cout << "TEST0705\n";
  cout << "  For the Fisher PDF:\n";
  cout << "  FISHER_PDF evaluates the PDF.\n";
  cout << "  FISHER_SAMPLE samples the PDF.\n";

  for ( test = 1; test <= test_num; test++ )
  {
    if ( test == 1 )
    {
      kappa = 0.0;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }
    else if ( test == 2 )
    {
      kappa = 0.5;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }
    else if ( test == 3 )
    {
      kappa = 10.0;
      mu[0] = 1.0;
      mu[1] = 0.0;
      mu[2] = 0.0;
    }

    cout << "\n";
    cout << "  PDF parameters:\n";
    cout << "    Concentration parameter KAPPA = " << kappa << "\n";
    cout << "    Direction MU(1:3) = "
         << "  " << mu[0]
         << "  " << mu[1]
         << "  " << mu[2] << "\n";

    cout << "\n";
    cout << "      X                         PDF\n";
    cout << "\n";

    seed = 123456789;
    x = fisher_sample ( kappa, mu, n, &seed );

    for ( j = 0; j < n; j++ )
    {
      pdf = fisher_pdf ( x+j*3, kappa, mu );

      cout << "  " << setw(10) << x[0+j*3]
           << "  " << setw(10) << x[1+j*3]
           << "  " << setw(10) << x[2+j*3]
           << "  " << setw(14) << pdf << "\n";
    }
    delete [] x;
  }
  return;
}
//****************************************************************************80

void test071 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST071 tests FISK_CDF, FISK_CDF_INV, FISK_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST071\n";
  cout << "  For the Fisk PDF:\n";
  cout << "  FISK_CDF evaluates the CDF;\n";
  cout << "  FISK_CDF_INV inverts the CDF.\n";
  cout << "  FISK_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !fisk_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST071 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = fisk_sample ( a, b, c, &seed );
    pdf = fisk_pdf ( x, a, b, c );
    cdf = fisk_cdf ( x, a, b, c );
    x2 = fisk_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test072 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST072 tests FISK_MEAN, FISK_SAMPLE, FISK_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST072\n";
  cout << "  For the Fisk PDF:\n";
  cout << "  FISK_MEAN computes the mean;\n";
  cout << "  FISK_SAMPLE samples;\n";
  cout << "  FISK_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !fisk_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST072 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = fisk_mean ( a, b, c );
  variance = fisk_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = fisk_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test073 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST073 tests FOLDED_NORMAL_CDF, FOLDED_NORMAL_CDF_INV, FOLDED_NORMAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST073\n";
  cout << "  For the Folded Normal PDF:\n";
  cout << "  FOLDED_NORMAL_CDF evaluates the CDF;\n";
  cout << "  FOLDED_NORMAL_CDF_INV inverts the CDF.\n";
  cout << "  FOLDED_NORMAL_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !folded_normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST073 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = folded_normal_sample ( a, b, &seed );
    pdf = folded_normal_pdf ( x, a, b );
    cdf = folded_normal_cdf ( x, a, b );
    x2 = folded_normal_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test074 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST074 tests FOLDED_NORMAL_MEAN, FOLDED_NORMAL_SAMPLE, FOLDED_NORMAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST074\n";
  cout << "  For the Folded Normal PDF:\n";
  cout << "  FOLDED_NORMAL_MEAN computes the mean;\n";
  cout << "  FOLDED_NORMAL_SAMPLE samples;\n";
  cout << "  FOLDED_NORMAL_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !folded_normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST074 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = folded_normal_mean ( a, b );
  variance = folded_normal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = folded_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0744 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0744 tests FRECHET_CDF, FRECHET_CDF_INV, FRECHET_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double alpha;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST0744\n";
  cout << "  For the Frechet PDF:\n";
  cout << "  FRECHET_CDF evaluates the CDF;\n";
  cout << "  FRECHET_CDF_INV inverts the CDF.\n";
  cout << "  FRECHET_PDF evaluates the PDF;\n";

  alpha = 3.0;

  cout << "\n";
  cout << "  PDF parameter ALPHA =  " << alpha << "\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = frechet_sample ( alpha, &seed );
    pdf = frechet_pdf ( x, alpha );
    cdf = frechet_cdf ( x, alpha );
    x2 = frechet_cdf_inv ( cdf, alpha );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test0745 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0745 tests FRECHET_MEAN, FRECHET_SAMPLE, FRECHET_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2008
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double alpha;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST0745\n";
  cout << "  For the Frechet PDF:\n";
  cout << "  FRECHET_MEAN computes the mean;\n";
  cout << "  FRECHET_SAMPLE samples;\n";
  cout << "  FRECHET_VARIANCE computes the variance;\n";

  alpha = 3.0;

  cout << "\n";
  cout << "  PDF parameter ALPHA =  " << alpha << "\n";

  mean = frechet_mean ( alpha );
  variance = frechet_variance ( alpha );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = frechet_sample ( alpha, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test075 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST075 tests GAMMA, GAMMA_LOG, GAMMA_LOG_INT, D_FACTORIAL.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double g1;
  double g2;
  double g3;
  double g4;
  int i;
  double x;

  cout << "\n";
  cout << "TEST075\n";
  cout << "  GAMMA evaluates the Gamma function;\n";
  cout << "  GAMMA_LOG evaluates the log of the Gamma function;\n";
  cout << "  GAMMA_LOG_INT evaluates the log for integer argument;\n";
  cout << "  D_FACTORIAL evaluates the factorial function.\n";

  cout << "\n";
  cout << "  X, GAMMA(X), Exp(GAMMA_LOG(X)), Exp(GAMMA_LOG_INT(X)) ";
  cout << "D_FACTORIAL(X+1)\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = ( double ) i;
    g1 = gamma ( x );
    g2 = exp ( gamma_log ( x ) );
    g3 = exp ( gamma_log_int ( i ) );
    g4 = i4_factorial ( i - 1 );
    cout << "  "
         << setw(6)  << x  << "  "
         << setw(14) << g1 << "  "
         << setw(14) << g2 << "  "
         << setw(14) << g3 << "  "
         << setw(14) << g4 << "\n";
  }

  return;
}
//****************************************************************************80

void test076 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST076 tests GAMMA_INC, GAMMA_INC_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST076:\n";
  cout << "  GAMMA_INC evaluates the normalized incomplete Gamma\n";
  cout << "    function GAMMA_INC(A,B,X).\n";
  cout << "  GAMMA_INC_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   A      X       Exact F       GAMMA_INC(A,X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    gamma_inc_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = gamma_inc ( a, x );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test077 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST077 tests GAMMA_CDF, GAMMA_PDF, GAMMA_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  cout << "\n";
  cout << "TEST077\n";
  cout << "  For the Gamma PDF:\n";
  cout << "  GAMMA_CDF evaluates the CDF;\n";
  cout << "  GAMMA_PDF evaluates the PDF;\n";
  cout << "  GAMMA_SAMPLE samples the PDF;\n";

  a = 1.0;
  b = 1.5;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A = " << a   << "\n";
  cout << "  PDF parameter B = " << b   << "\n";
  cout << "  PDF parameter B = " << c   << "\n";

  if ( !gamma_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST077 - Fatal error!\n";
    cout << "  The parameter values are illegal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = gamma_sample ( a, b, c, &seed );
    pdf = gamma_pdf ( x, a, b, c );
    cdf = gamma_cdf ( x, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test078 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST078 tests GAMMA_MEAN, GAMMA_SAMPLE, GAMMA_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST078\n";
  cout << "  For the Gamma PDF:\n";
  cout << "  GAMMA_MEAN computes the mean;\n";
  cout << "  GAMMA_SAMPLE samples;\n";
  cout << "  GAMMA_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 3.0;
  c = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !gamma_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST078 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = gamma_mean ( a, b, c );
  variance = gamma_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gamma_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test079 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST079 tests GENLOGISTIC_CDF, GENLOGISTIC_CDF_INV, GENLOGISTIC_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST079\n";
  cout << "  For the Genlogistic PDF:\n";
  cout << "  GENLOGISTIC_CDF evaluates the CDF;\n";
  cout << "  GENLOGISTIC_CDF_INV inverts the CDF.\n";
  cout << "  GENLOGISTIC_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !genlogistic_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST079 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = genlogistic_sample ( a, b, c, &seed );
    pdf = genlogistic_pdf ( x, a, b, c );
    cdf = genlogistic_cdf ( x, a, b, c );
    x2 = genlogistic_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test080 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST080 tests GENLOGISTIC_MEAN, GENLOGISTIC_SAMPLE, GENLOGISTIC_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST080\n";
  cout << "  For the Genlogistic PDF:\n";
  cout << "  GENLOGISTIC_MEAN computes the mean;\n";
  cout << "  GENLOGISTIC_SAMPLE samples;\n";
  cout << "  GENLOGISTIC_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !genlogistic_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST080 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = genlogistic_mean ( a, b, c );
  variance = genlogistic_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = genlogistic_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test081 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST081 tests GEOMETRIC_CDF, GEOMETRIC_CDF_INV, GEOMETRIC_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST081\n";
  cout << "  For the Geometric PDF:\n";
  cout << "  GEOMETRIC_CDF evaluates the CDF;\n";
  cout << "  GEOMETRIC_CDF_INV inverts the CDF.\n";
  cout << "  GEOMETRIC_PDF evaluates the PDF;\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !geometric_check ( a ) )
  {
    cout << "\n";
    cout << "TEST081 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = geometric_sample ( a, &seed );
    pdf = geometric_pdf ( x, a );
    cdf = geometric_cdf ( x, a );
    x2 = geometric_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test082 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST082 tests GEOMETRIC_MEAN, GEOMETRIC_SAMPLE, GEOMETRIC_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST082\n";
  cout << "  For the Geometric PDF:\n";
  cout << "  GEOMETRIC_MEAN computes the mean;\n";
  cout << "  GEOMETRIC_SAMPLE samples;\n";
  cout << "  GEOMETRIC_VARIANCE computes the variance.\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !geometric_check ( a ) )
  {
    cout << "\n";
    cout << "TEST082 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = geometric_mean ( a );
  variance = geometric_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = geometric_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test083 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST083 tests GEOMETRIC_CDF, GEOMETRIC_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST083\n";
  cout << "  For the Geometric PDF:\n";
  cout << "  GEOMETRIC_CDF evaluates the CDF;\n";
  cout << "  GEOMETRIC_PDF evaluates the PDF;\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !geometric_check ( a ) )
  {
    cout << "\n";
    cout << "TEST083 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 0; x <= 10; x++ )
  {
    pdf = geometric_pdf ( x, a );
    cdf = geometric_cdf ( x, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test084 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST084 tests GOMPERTZ_CDF, GOMPERTZ_CDF_INV, GOMPERTZ_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST084\n";
  cout << "  For the Gompertz PDF:\n";
  cout << "  GOMPERTZ_CDF evaluates the CDF;\n";
  cout << "  GOMPERTZ_CDF_INV inverts the CDF.\n";
  cout << "  GOMPERTZ_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !gompertz_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST084 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = gompertz_sample ( a, b, &seed );
    pdf = gompertz_pdf ( x, a, b );
    cdf = gompertz_cdf ( x, a, b );
    x2 = gompertz_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test085 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST085 tests GOMPERTZ_MEAN, GOMPERTZ_SAMPLE, GOMPERTZ_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST085\n";
  cout << "  For the Gompertz PDF:\n";
  cout << "  GOMPERTZ_MEAN computes the mean;\n";
  cout << "  GOMPERTZ_SAMPLE samples;\n";
  cout << "  GOMPERTZ_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !gompertz_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST085 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gompertz_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test086 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST086 tests GUMBEL_CDF, GUMBEL_CDF_INV, GUMBEL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST086\n";
  cout << "  For the Gumbel PDF:\n";
  cout << "  GUMBEL_CDF evaluates the CDF;\n";
  cout << "  GUMBEL_CDF_INV inverts the CDF.\n";
  cout << "  GUMBEL_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = gumbel_sample ( &seed );
    pdf = gumbel_pdf ( x );
    cdf = gumbel_cdf ( x );
    x2 = gumbel_cdf_inv ( cdf );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test087 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST087 tests GUMBEL_MEAN, GUMBEL_SAMPLE, GUMBEL_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST087\n";
  cout << "  For the Gumbel PDF:\n";
  cout << "  GUMBEL_MEAN computes the mean;\n";
  cout << "  GUMBEL_SAMPLE samples;\n";
  cout << "  GUMBEL_VARIANCE computes the variance.\n";

  mean     = gumbel_mean ( );
  variance = gumbel_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = gumbel_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test088 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST088 tests HALF_NORMAL_CDF, HALF_NORMAL_CDF_INV, HALF_NORMAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST088\n";
  cout << "  For the Half Normal PDF:\n";
  cout << "  HALF_NORMAL_CDF evaluates the CDF;\n";
  cout << "  HALF_NORMAL_CDF_INV inverts the CDF.\n";
  cout << "  HALF_NORMAL_PDF evaluates the PDF;\n";

  a = 0.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !half_normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST088 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = half_normal_sample ( a, b, &seed );
    pdf = half_normal_pdf ( x, a, b );
    cdf = half_normal_cdf ( x, a, b );
    x2 = half_normal_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test089 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST089 tests HALF_NORMAL_MEAN, HALF_NORMAL_SAMPLE, HALF_NORMAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST089\n";
  cout << "  For the Half Normal PDF:\n";
  cout << "  HALF_NORMAL_MEAN computes the mean;\n";
  cout << "  HALF_NORMAL_SAMPLE samples;\n";
  cout << "  HALF_NORMAL_VARIANCE computes the variance;\n";

  a = 0.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !half_normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST089 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = half_normal_mean ( a, b );
  variance = half_normal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = half_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test090 ( )

//****************************************************************************80
//
//  Purpose:
//
//  TEST090 tests HYPERGEOMETRIC_CDF, HYPERGEOMETRIC_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int l;
  int m;
  int n;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST090\n";
  cout << "  For the Hypergeometric PDF:\n";
  cout << "  HYPERGEOMETRIC_CDF evaluates the CDF.\n";
  cout << "  HYPERGEOMETRIC_PDF evaluates the PDF.\n";

  n = 100;
  m = 70;
  l = 1000;

  cout << "\n";
  cout << "  Total number of balls L =         " << l << "\n";
  cout << "  Number of white balls M =         " << m << "\n";
  cout << "  Number of balls taken N =         " << n << "\n";

  if ( !hypergeometric_check ( n, m, l ) )
  {
    cout << "\n";
    cout << "TEST090 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  x = 7;

  pdf = hypergeometric_pdf ( x, n, m, l );

  cdf = hypergeometric_cdf ( x, n, m, l );

  cout << "  PDF argument X =                " << x   << "\n";
  cout << "  PDF value =                   = " << pdf << "\n";
  cout << "  CDF value =                   = " << cdf << "\n";

  return;
}
//****************************************************************************80

void test091 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST091 tests HYPERGEOMETRIC_MEAN, HYPERGEOMETRIC_SAMPLE, HYPERGEOMETRIC_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  int j;
  int l;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST091\n";
  cout << "  For the Hypergeometric PDF:\n";
  cout << "  HYPERGEOMETRIC_MEAN computes the mean;\n";
  cout << "  HYPERGEOMETRIC_SAMPLE samples;\n";
  cout << "  HYPERGEOMETRIC_VARIANCE computes the variance.\n";

  n = 100;
  m = 70;
  l = 1000;

  cout << "\n";
  cout << "  Total number of balls L =         " << l << "\n";
  cout << "  Number of white balls M =         " << m << "\n";
  cout << "  Number of balls taken N =         " << n << "\n";

  if ( !hypergeometric_check ( n, m, l ) )
  {
    cout << "\n";
    cout << "TEST090 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = hypergeometric_mean ( n, m, l );
  variance = hypergeometric_variance ( n, m, l );

  cout << "  PDF mean =                    " << mean << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  cout << "\n";
  cout << "THIS CALL IS TAKING FOREVER!\n";
  return;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = hypergeometric_sample ( n, m, l, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test092 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST092 tests R8_CEILING.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int ival;
  double rval;

  cout << "\n";
  cout << "TEST092\n";
  cout << "  R8_CEILING rounds an R8 up.\n";
  cout << "\n";
  cout << "       X           R8_CEILING(X)\n";
  cout << "\n";

  for ( i = -6; i <= 6; i++ )
  {
    rval = ( double ) ( i ) / 5.0;
    ival = r8_ceiling ( rval );
    cout << "  "
         << setw(14) << rval << "  "
         << setw(6)  << ival << "\n";
  }

  return;
}
//****************************************************************************80

void test093 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST093 tests INVERSE_GAUSSIAN_CDF, INVERSE_GAUSSIAN_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;

  cout << "\n";
  cout << "TEST093\n";
  cout << "  For the Inverse Gaussian PDF:\n";
  cout << "  INVERSE_GAUSSIAN_CDF evaluates the CDF;\n";
  cout << "  INVERSE_GAUSSIAN_PDF evaluates the PDF;\n";

  a = 5.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !inverse_gaussian_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST093 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = inverse_gaussian_sample ( a, b, &seed );
    pdf = inverse_gaussian_pdf ( x, a, b );
    cdf = inverse_gaussian_cdf ( x, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test094 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST094 tests INVERSE_GAUSSIAN_MEAN, INVERSE_GAUSSIAN_SAMPLE, INVERSE_GAUSSIAN_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST094\n";
  cout << "  For the Inverse Gaussian PDF:\n";
  cout << "  INVERSE_GAUSSIAN_MEAN computes the mean;\n";
  cout << "  INVERSE_GAUSSIAN_SAMPLE samples;\n";
  cout << "  INVERSE_GAUSSIAN_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !inverse_gaussian_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST094 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = inverse_gaussian_mean ( a, b );
  variance = inverse_gaussian_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = inverse_gaussian_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test095 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST095 tests LAPLACE_CDF, LAPLACE_CDF_INV, LAPLACE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST095\n";
  cout << "  For the Laplace PDF:\n";
  cout << "  LAPLACE_CDF evaluates the CDF;\n";
  cout << "  LAPLACE_CDF_INV inverts the CDF.\n";
  cout << "  LAPLACE_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !laplace_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST095 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = laplace_sample ( a, b, &seed );
    pdf = laplace_pdf ( x, a, b );
    cdf = laplace_cdf ( x, a, b );
    x2 = laplace_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test096 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST096 tests LAPLACE_MEAN, LAPLACE_SAMPLE, LAPLACE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST096\n";
  cout << "  For the Laplace PDF:\n";
  cout << "  LAPLACE_MEAN computes the mean;\n";
  cout << "  LAPLACE_SAMPLE samples;\n";
  cout << "  LAPLACE_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !laplace_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST096 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = laplace_mean ( a, b );
  variance = laplace_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = laplace_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test0965 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST0965 tests LEVY_CDF, LEVY_CDF_INV, LEVY_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST0965\n";
  cout << "  For the Levy PDF:\n";
  cout << "  LEVY_CDF evaluates the CDF;\n";
  cout << "  LEVY_CDF_INV inverts the CDF.\n";
  cout << "  LEVY_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = levy_sample ( a, b, &seed );
    pdf = levy_pdf ( x, a, b );
    cdf = levy_cdf ( x, a, b );
    x2 = levy_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test097 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST097 tests LOGISTIC_CDF, LOGISTIC_CDF_INV, LOGISTIC_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST097\n";
  cout << "  For the Logistic PDF:\n";
  cout << "  LOGISTIC_CDF evaluates the CDF;\n";
  cout << "  LOGISTIC_CDF_INV inverts the CDF.\n";
  cout << "  LOGISTIC_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !logistic_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST097 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = logistic_sample ( a, b, &seed );
    pdf = logistic_pdf ( x, a, b );
    cdf = logistic_cdf ( x, a, b );
    x2 = logistic_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test098 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST098 tests LOGISTIC_MEAN, LOGISTIC_SAMPLE, LOGISTIC_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST098\n";
  cout << "  For the Logistic PDF:\n";
  cout << "  LOGISTIC_MEAN computes the mean;\n";
  cout << "  LOGISTIC_SAMPLE samples;\n";
  cout << "  LOGISTIC_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !logistic_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST098 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = logistic_mean ( a, b );
  variance = logistic_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = logistic_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test099 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST099 tests LOG_NORMAL_CDF, LOG_NORMAL_CDF_INV, LOG_NORMAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST099\n";
  cout << "  For the Log Normal PDF:\n";
  cout << "  LOG_NORMAL_CDF evaluates the CDF;\n";
  cout << "  LOG_NORMAL_CDF_INV inverts the CDF.\n";
  cout << "  LOG_NORMAL_PDF evaluates the PDF;\n";

  a = 10.0;
  b = 2.25;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !log_normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST099 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = log_normal_sample ( a, b, &seed );
    pdf = log_normal_pdf ( x, a, b );
    cdf = log_normal_cdf ( x, a, b );
    x2 = log_normal_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test100 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST100 tests LOG_NORMAL_MEAN, LOG_NORMAL_SAMPLE, LOG_NORMAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST100\n";
  cout << "  For the LogNormal PDF:\n";
  cout << "  LOG_NORMAL_MEAN computes the mean;\n";
  cout << "  LOG_NORMAL_SAMPLE samples;\n";
  cout << "  LOG_NORMAL_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST100 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = log_normal_mean ( a, b);
  variance = log_normal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = log_normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test101 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST101 tests LOG_SERIES_CDF, LOG_SERIES_CDF_INV, LOG_SERIES_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST101\n";
  cout << "  For the Log Series PDF:\n";
  cout << "  LOG_SERIES_CDF evaluates the CDF;\n";
  cout << "  LOG_SERIES_CDF_INV inverts the CDF.\n";
  cout << "  LOG_SERIES_PDF evaluates the PDF;\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !log_series_check ( a ) )
  {
    cout << "\n";
    cout << "TEST101 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = log_series_sample ( a, &seed );
    pdf = log_series_pdf ( x, a );
    cdf = log_series_cdf ( x, a );
    x2 = log_series_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test102 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST102 tests LOG_SERIES_CDF, LOG_SERIES_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST102\n";
  cout << "  For the Log Series PDF:\n";
  cout << "  LOG_SERIES_CDF evaluates the CDF;\n";
  cout << "  LOG_SERIES_PDF evaluates the PDF;\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !log_series_check ( a ) )
  {
    cout << "\n";
    cout << "TEST101 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 1; x <= 10; x++ )
  {
    pdf = log_series_pdf ( x, a );
    cdf = log_series_cdf ( x, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test103 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST103 tests LOG_SERIES_MEAN, LOG_SERIES_SAMPLE, LOG_SERIES_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST103\n";
  cout << "  For the Log Series PDF:\n";
  cout << "  LOG_SERIES_MEAN computes the mean;\n";
  cout << "  LOG_SERIES_SAMPLE samples;\n";
  cout << "  LOG_SERIES_VARIANCE computes the variance.\n";

  a = 0.25E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !log_series_check ( a ) )
  {
    cout << "\n";
    cout << "TEST103 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = log_series_mean ( a );
  variance = log_series_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = log_series_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test104 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST104 tests LOG_UNIFORM_CDF, LOG_UNIFORM_CDF_INV, LOG_UNIFORM_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST104\n";
  cout << "  For the Log Uniform PDF:\n";
  cout << "  LOG_UNIFORM_CDF evaluates the CDF;\n";
  cout << "  LOG_UNIFORM_CDF_INV inverts the CDF.\n";
  cout << "  LOG_UNIFORM_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 20.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !log_uniform_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST104 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = log_uniform_sample ( a, b, &seed );
    pdf = log_uniform_pdf ( x, a, b );
    cdf = log_uniform_cdf ( x, a, b );
    x2 = log_uniform_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test105 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST105 tests LOG_UNIFORM_MEAN, LOG_UNIFORM_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST105\n";
  cout << "  For the Log Uniform PDF:\n";
  cout << "  LOG_UNIFORM_MEAN computes the mean;\n";
  cout << "  LOG_UNIFORM_SAMPLE samples;\n";

  a = 2.0;
  b = 20.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !log_uniform_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST105 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = log_uniform_mean ( a, b);

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = log_uniform_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test106 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST106 tests LORENTZ_CDF, LORENTZ_CDF_INV, LORENTZ_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//

{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST106\n";
  cout << "  For the Lorentz PDF:\n";
  cout << "  LORENTZ_CDF evaluates the CDF;\n";
  cout << "  LORENTZ_CDF_INV inverts the CDF.\n";
  cout << "  LORENTZ_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = lorentz_sample ( &seed );
    pdf = lorentz_pdf ( x );
    cdf = lorentz_cdf ( x );
    x2 = lorentz_cdf_inv ( cdf );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test107 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST107 tests LORENTZ_MEAN, LORENTZ_SAMPLE, LORENTZ_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST107\n";
  cout << "  For the Lorentz PDF:\n";
  cout << "  LORENTZ_MEAN computes the mean;\n";
  cout << "  LORENTZ_SAMPLE samples;\n";
  cout << "  LORENTZ_VARIANCE computes the variance.\n";

  mean     = lorentz_mean ( );
  variance = lorentz_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = lorentz_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test108 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST108 tests MAXWELL_CDF, MAXWELL_CDF_INV, MAXWELL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST108\n";
  cout << "  For the Maxwell PDF:\n";
  cout << "  MAXWELL_CDF evaluates the CDF;\n";
  cout << "  MAXWELL_CDF_INV inverts the CDF.\n";
  cout << "  MAXWELL_PDF evaluates the PDF;\n";

  a = 2.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !maxwell_check ( a ) )
  {
    cout << "\n";
    cout << "TEST108 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }
  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = maxwell_sample ( a, &seed );
    pdf = maxwell_pdf ( x, a );
    cdf = maxwell_cdf ( x, a );
    x2 = maxwell_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test109 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST109 tests MAXWELL_MEAN, MAXWELL_SAMPLE, MAXWELL_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST109\n";
  cout << "  For the Maxwell PDF:\n";
  cout << "  MAXWELL_MEAN computes the mean;\n";
  cout << "  MAXWELL_SAMPLE samples;\n";
  cout << "  MAXWELL_VARIANCE computes the variance.\n";

  a = 2.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !maxwell_check ( a ) )
  {
    cout << "\n";
    cout << "TEST109 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = maxwell_mean ( a );
  variance = maxwell_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = maxwell_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test110 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST110 tests MULTINOMIAL_COEF1, MULTINOMIAL_COEF2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define MAXFACTOR 5

  int factor[MAXFACTOR];
  int i;
  int j;
  int n;
  int ncomb1;
  int ncomb2;
  int nfactor;

  cout << "\n";
  cout << "TEST110\n";
  cout << "  MULTINOMIAL_COEF1 computes multinomial\n";
  cout << "    coefficients using the Gamma function;\n";
  cout << "  MULTINOMIAL_COEF2 computes multinomial\n";
  cout << "    coefficients directly.\n";

  cout << "\n";
  cout << "  Line 10 of the BINOMIAL table:\n";
  cout << "\n";

  n = 10;
  nfactor = 2;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;
    factor[1] = n - i;

    ncomb1 = multinomial_coef1 ( nfactor, factor );

    ncomb2 = multinomial_coef2 ( nfactor, factor );

    cout << "  "
         << setw(2) << factor[0] << "  "
         << setw(2) << factor[1] << "  "
         << setw(5) << ncomb1    << "  "
         << setw(5) << ncomb2    << "\n";
  }

  cout << "\n";
  cout << "  Level 5 of the TRINOMIAL coefficients:\n";

  n = 5;
  nfactor = 3;

  for ( i = 0; i <= n; i++ )
  {
    factor[0] = i;

    cout << "\n";

    for ( j = 0; j <= n - factor[0]; j++ )
    {
      factor[1] = j;
      factor[2] = n - factor[0] - factor[1];

      ncomb1 = multinomial_coef1 ( nfactor, factor );

      ncomb2 = multinomial_coef2 ( nfactor, factor );

      cout << "  "
           << setw(2) << factor[0] << "  "
           << setw(2) << factor[1] << "  "
           << setw(2) << factor[2] << "  "
           << setw(5) << ncomb1    << "  "
           << setw(5) << ncomb2    << "\n";
    }
  }

  return;
# undef MAXFACTOR
}
//****************************************************************************80

void test111 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST111 tests MULTINOMIAL_MEAN, MULTINOMIAL_SAMPLE, MULTINOMIAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define B 3
# define SAMPLE_NUM 1000

  int a;
  double c[B] = { 0.125, 0.500, 0.375 };
  int i;
  int j;
  double *mean;
  int seed = 123456789;
  double *variance;
  int x[B*SAMPLE_NUM];
  int *xmax;
  int *xmin;
  int *y;

  cout << "\n";
  cout << "TEST111\n";
  cout << "  For the Multinomial PDF:\n";
  cout << "  MULTINOMIAL_MEAN computes the mean;\n";
  cout << "  MULTINOMIAL_SAMPLE samples;\n";
  cout << "  MULTINOMIAL_VARIANCE computes the variance;\n";

  a = 5;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << B << "\n";
  r8vec_print ( B, c, "  PDF parameter C:" );

  if ( !multinomial_check ( a, B, c ) )
  {
    cout << "\n";
    cout << "TEST111 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = multinomial_mean ( a, B, c );
  variance = multinomial_variance ( a, B, c );
  r8vec_print ( B, mean, "  PDF mean:" );
  r8vec_print ( B, variance, "  PDF variance:" );

  delete [] mean;
  delete [] variance;

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    y = multinomial_sample ( a, B, c, &seed );
    for ( i = 0; i < B; i++ )
    {
      x[i+j*B] = y[i];
    }
    delete [] y;
  }

  mean = i4row_mean ( B, SAMPLE_NUM, x );
  variance = i4row_variance ( B, SAMPLE_NUM, x );
  xmax = i4row_max ( B, SAMPLE_NUM, x );
  xmin = i4row_min ( B, SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "\n";
  cout << "  Component Mean, Variance, Min, Max:\n";
  cout << "\n";

  for ( i = 0; i < B; i++ )
  {
    cout << "  "
         << setw(6)  << i+1         << "  "
         << setw(12) << mean[i]     << "  "
         << setw(12) << variance[i] << "  "
         << setw(12) << xmin[i]     << "  "
         << setw(12) << xmax[i]     << "\n";
  }

  delete [] mean;
  delete [] variance;
  delete [] xmax;
  delete [] xmin;

  return;
# undef B
# undef SAMPLE_NUM
}
//****************************************************************************80

void test112 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST112 tests MULTINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define B 3

  int a;
  double c[B] = { 0.1, 0.5, 0.4 };
  int i;
  double pdf;
  int x[B] = { 0, 2, 3 };

  cout << "\n";
  cout << "TEST112\n";
  cout << "  For the Multinomial PDF:\n";
  cout << "  MULTINOMIAL_PDF evaluates the PDF;\n";

  a = 5;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << B << "\n";
  r8vec_print ( B, c, "  PDF parameter C:" );

  if ( !multinomial_check ( a, B, c ) )
  {
    cout << "\n";
    cout << "TEST112 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  i4vec_print ( B, x, "  PDF argument X:" );

  pdf = multinomial_pdf ( x, a, B, c );

  cout << "\n";
  cout << "  PDF value = " << pdf<< "\n";

  return;
# undef B
# undef SAMPLE_NUM
}
//****************************************************************************80

void test113 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST113 tests NAKAGAMI_CDF, NAKAGAMI_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  double pdf;
  double x;

  cout << "\n";
  cout << "TEST113\n";
  cout << "  For the Nakagami PDF:\n";
  cout << "  NAKAGAMI_CDF evaluates the CDF;\n";
  cout << "  NAKAGAMI_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !nakagami_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST113 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  x = 1.25;
  pdf = nakagami_pdf ( x, a, b, c );
  cdf = nakagami_cdf ( x, a, b, c );

  cout << "  "
       << setw(12) << x   << "  "
       << setw(12) << pdf << "  "
       << setw(12) << cdf << "\n";

  return;
}
//****************************************************************************80

void test114 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST114 tests NAKAGAMI_MEAN, NAKAGAMI_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double mean;
  double variance;
  double x;

  cout << "\n";
  cout << "TEST114\n";
  cout << "  For the Nakagami PDF:\n";
  cout << "  NAKAGAMI_MEAN evaluates the mean;\n";
  cout << "  NAKAGAMI_VARIANCE evaluates the variance;\n";

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !nakagami_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST114 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = nakagami_mean ( a, b, c );
  variance = nakagami_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =      " << mean     << "\n";
  cout << "  PDF variance =  " << variance << "\n";

  return;
}
//****************************************************************************80

void test1145 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1145 tests NEGATIVE_BINOMIAL_CDF, NEGATIVE_BINOMIAL_CDF_INV, NEGATIVE_BINOMIAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  double b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST1145\n";
  cout << "  For the Negative Binomial PDF:\n";
  cout << "  NEGATIVE_BINOMIAL_CDF evaluates the CDF;\n";
  cout << "  NEGATIVE_BINOMIAL_CDF_INV inverts the CDF.\n";
  cout << "  NEGATIVE_BINOMIAL_PDF evaluates the PDF;\n";

  a = 2;
  b = 0.25;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !negative_binomial_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST1145 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = negative_binomial_sample ( a, b, &seed );
    pdf = negative_binomial_pdf ( x, a, b );
    cdf = negative_binomial_cdf ( x, a, b );
    x2 = negative_binomial_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test1146 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1146 tests NEGATIVE_BINOMIAL_MEAN, NEGATIVE_BINOMIAL_SAMPLE, NEGATIVE_BINOMIAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST1146\n";
  cout << "  For the Negative Binomial PDF:\n";
  cout << "  NEGATIVE_BINOMIAL_MEAN computes the mean;\n";
  cout << "  NEGATIVE_BINOMIAL_SAMPLE samples;\n";
  cout << "  NEGATIVE_BINOMIAL_VARIANCE computes the variance;\n";

  a = 2;
  b = 0.75;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !negative_binomial_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST1146 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = negative_binomial_mean ( a, b );
  variance = negative_binomial_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = negative_binomial_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test115 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST115 tests NORMAL_01_CDF, NORMAL_01_CDF_INV, NORMAL_01_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST115\n";
  cout << "  For the Normal 01 PDF:\n";
  cout << "  NORMAL_01_CDF evaluates the CDF;\n";
  cout << "  NORMAL_01_CDF_INV inverts the CDF.\n";
  cout << "  NORMAL_01_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_01_sample ( &seed );
    pdf = normal_01_pdf ( x );
    cdf = normal_01_cdf ( x );
    x2 = normal_01_cdf_inv ( cdf );

    cout << "  "
         << setw(24) << setprecision(16) << x   << "  "
         << setw(12) << setprecision(6)  << pdf << "  "
         << setw(12) << setprecision(6)  << cdf << "  "
         << setw(24) << setprecision(16) << x2  << "\n";
  }

  cout << setprecision(6);

  return;
}
//****************************************************************************80

void test116 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST116 tests NORMAL_01_MEAN, NORMAL_01_SAMPLE, NORMAL_01_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST116\n";
  cout << "  For the Normal 01 PDF:\n";
  cout << "  NORMAL_01_MEAN computes the mean;\n";
  cout << "  NORMAL_01_SAMPLE samples;\n";
  cout << "  NORMAL_01_VARIANCE computes the variance;\n";

  mean = normal_01_mean ( );
  variance = normal_01_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = normal_01_sample ( &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test117 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST117 tests NORMAL_CDF, NORMAL_CDF_INV, NORMAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST117\n";
  cout << "  For the Normal PDF:\n";
  cout << "  NORMAL_CDF evaluates the CDF;\n";
  cout << "  NORMAL_CDF_INV inverts the CDF.\n";
  cout << "  NORMAL_PDF evaluates the PDF;\n";

  a = 100.0;
  b = 15.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST117 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = normal_sample ( a, b, &seed );
    pdf = normal_pdf ( x, a, b );
    cdf = normal_cdf ( x, a, b );
    x2 = normal_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test118 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST118 tests NORMAL_MEAN, NORMAL_SAMPLE, NORMAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST118\n";
  cout << "  For the Normal PDF:\n";
  cout << "  NORMAL_MEAN computes the mean;\n";
  cout << "  NORMAL_SAMPLE samples;\n";
  cout << "  NORMAL_VARIANCE computes the variance;\n";

  a = 100.0;
  b = 15.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !normal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST118 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = normal_mean ( a, b );
  variance = normal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = normal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test119 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST119 tests PARETO_CDF, PARETO_CDF_INV, PARETO_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST119\n";
  cout << "  For the Pareto PDF:\n";
  cout << "  PARETO_CDF evaluates the CDF;\n";
  cout << "  PARETO_CDF_INV inverts the CDF.\n";
  cout << "  PARETO_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !pareto_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST119 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = pareto_sample ( a, b, &seed );
    pdf = pareto_pdf ( x, a, b );
    cdf = pareto_cdf ( x, a, b );
    x2 = pareto_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test120 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST120 tests PARETO_MEAN, PARETO_SAMPLE, PARETO_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST120\n";
  cout << "  For the Pareto PDF:\n";
  cout << "  PARETO_MEAN computes the mean;\n";
  cout << "  PARETO_SAMPLE samples;\n";
  cout << "  PARETO_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !pareto_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST120 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = pareto_mean ( a, b );
  variance = pareto_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = pareto_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test123 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST123 tests PEARSON_05_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double pdf;
  double x;

  cout << "\n";
  cout << "TEST123\n";
  cout << "  For the Pearson 05 PDF:\n";
  cout << "  PEARSON_05_PDF evaluates the PDF.\n";

  x = 5.0;

  a = 1.0;
  b = 2.0;
  c = 3.0;

  cout << "\n";
  cout << "  PDF parameter A = " << a   << "\n";
  cout << "  PDF parameter B = " << b   << "\n";
  cout << "  PDF parameter C = " << c   << "\n";

  if ( !pearson_05_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST123 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  pdf = pearson_05_pdf ( x, a, b, c );

  cout << "\n";
  cout << "  PDF argument X =  " << x   << "\n";
  cout << "  PDF value =       " << pdf << "\n";

  return;
}
//****************************************************************************80

void test124 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST124 tests PLANCK_PDF, PLANCK_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  cout << "\n";
  cout << "TEST124\n";
  cout << "  For the Planck PDF:\n";
  cout << "  PLANCK_PDF evaluates the PDF.\n";
  cout << "  PLANCK_SAMPLE samples the PDF.\n";

  a = 2.0E+00;
  b = 3.0E+00;

  cout << "\n";
  cout << "  PDF parameter A = " << a << "\n";
  cout << "  PDF parameter B = " << b << "\n";

  if ( !planck_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST124 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = planck_sample ( a, b, &seed );

    pdf = planck_pdf ( x, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "\n";
  }

  return;
}
//****************************************************************************80

void test125 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST125 tests PLANCK_MEAN, PLANCK_SAMPLE, PLANCK_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST125\n";
  cout << "  For the Planck PDF:\n";
  cout << "  PLANCK_MEAN computes the mean;\n";
  cout << "  PLANCK_SAMPLE samples;\n";
  cout << "  PLANCK_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !planck_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST125 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = planck_mean ( a, b );
  variance = planck_variance ( a, b );

  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = planck_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test126 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST126 tests POISSON_CDF, POISSON_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double fx;
  double fx2;
  int n_data;
  int x;

  cout << "\n";
  cout << "TEST126:\n";
  cout << "  POISSON_CDF evaluates the cumulative distribution\n";
  cout << "    function for the discrete Poisson probability\n";
  cout << "    density function.\n";
  cout << "  POISSON_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  A is the expected mean number of successes per unit time;\n";
  cout << "  X is the number of successes;\n";
  cout << "  POISSON_CDF is the probability of having up to X\n";
  cout << "  successes in unit time.\n";
  cout << "\n";
  cout << "   A          X   Exact F     POISSON_CDF(A,X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    poisson_cdf_values ( &n_data, &a, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = poisson_cdf ( x, a );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test127 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST127 tests POISSON_CDF, POISSON_CDF_INV, POISSON_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST127\n";
  cout << "  For the Poisson PDF:\n";
  cout << "  POISSON_CDF evaluates the CDF;\n";
  cout << "  POISSON_CDF_INV inverts the CDF.\n";
  cout << "  POISSON_PDF evaluates the PDF;\n";

  a = 10.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !poisson_check ( a ) )
  {
    cout << "\n";
    cout << "TEST127 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = poisson_sample ( a, &seed );
    pdf = poisson_pdf ( x, a );
    cdf = poisson_cdf ( x, a );
    x2 = poisson_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test128 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST128 tests POISSON_MEAN, POISSON_SAMPLE, POISSON_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST128\n";
  cout << "  For the Poisson PDF,\n";
  cout << "  POISSON_SAMPLE samples the Poisson PDF.\n";
  cout << "  POISSON_SAMPLE samples the Poisson PDF.\n";
  cout << "  POISSON_SAMPLE samples the Poisson PDF.\n";

  a = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";

  if ( !poisson_check ( a ) )
  {
    cout << "\n";
    cout << "TEST128 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = poisson_mean ( a );
  variance = poisson_variance ( a );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = poisson_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test129 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST129 tests POWER_CDF, POWER_CDF_INV, POWER_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST129\n";
  cout << "  For the Power PDF:\n";
  cout << "  POWER_CDF evaluates the CDF;\n";
  cout << "  POWER_CDF_INV inverts the CDF.\n";
  cout << "  POWER_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !power_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST129 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = power_sample ( a, b, &seed );
    pdf = power_pdf ( x, a, b );
    cdf = power_cdf ( x, a, b );
    x2 = power_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test130 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST130 tests POWER_MEAN, POWER_SAMPLE, POWER_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST130\n";
  cout << "  For the Power PDF:\n";
  cout << "  POWER_MEAN computes the mean;\n";
  cout << "  POWER_SAMPLE samples;\n";
  cout << "  POWER_VARIANCE computes the variance;\n";

  a = 2.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !power_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST130 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = power_mean ( a, b );
  variance = power_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = power_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test1304 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1304 tests QUASIGEOMETRIC_CDF, *_CDF_INV, *_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
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
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST1304\n";
  cout << "  For the Quasigeometric PDF:\n";
  cout << "  QUASIGEOMETRIC_CDF evaluates the CDF;\n";
  cout << "  QUASIGEOMETRIC_CDF_INV inverts the CDF.\n";
  cout << "  QUASIGEOMETRIC_PDF evaluates the PDF;\n";

  a = 0.4825;
  b = 0.5893;

  cout << "\n";
  cout << "  PDF parameter A = " << a   << "\n";
  cout << "  PDF parameter B = " << b   << "\n";

  if ( !quasigeometric_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST1304 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = quasigeometric_sample ( a, b, &seed );
    pdf = quasigeometric_pdf ( x, a, b );
    cdf = quasigeometric_cdf ( x, a, b );
    x2 = quasigeometric_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test1306 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1306 tests QUASIGEOMETRIC_MEAN, *_SAMPLE, *_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 January 2009
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  int j;
  double mean;
  int sample_num = 1000;
  int seed = 123456789;
  double variance;
  int *x;
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST1306\n";
  cout << "  For the Quasigeometric PDF:\n";
  cout << "  QUASIGEOMETRIC_MEAN computes the mean;\n";
  cout << "  QUASIGEOMETRIC_SAMPLE samples;\n";
  cout << "  QUASIGEOMETRIC_VARIANCE computes the variance.\n";

  a = 0.4825;
  b = 0.5893;

  cout << "\n";
  cout << "  PDF parameter A = " << a << "\n";
  cout << "  PDF parameter B = " << b << "\n";

  if ( !quasigeometric_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST1306 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = quasigeometric_mean ( a, b );
  variance = quasigeometric_variance ( a, b );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  x = new int[sample_num];

  for ( j = 0; j < sample_num; j++ )
  {
    x[j] = quasigeometric_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( sample_num, x );
  variance = i4vec_variance ( sample_num, x );
  xmax = i4vec_max ( sample_num, x );
  xmin = i4vec_min ( sample_num, x );

  cout << "\n";
  cout << "  Sample size =     " << sample_num  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  delete [] x;

  return;
}
//****************************************************************************80

void test131 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST131 tests RAYLEIGH_CDF, RAYLEIGH_CDF_INV, RAYLEIGH_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST131\n";
  cout << "  For the Rayleigh PDF:\n";
  cout << "  RAYLEIGH_CDF evaluates the CDF;\n";
  cout << "  RAYLEIGH_CDF_INV inverts the CDF.\n";
  cout << "  RAYLEIGH_PDF evaluates the PDF;\n";

  a = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !rayleigh_check ( a ) )
  {
    cout << "\n";
    cout << "TEST131 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = rayleigh_sample ( a, &seed );
    pdf = rayleigh_pdf ( x, a );
    cdf = rayleigh_cdf ( x, a );
    x2 = rayleigh_cdf_inv ( cdf, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test132 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST132 tests RAYLEIGH_MEAN, RAYLEIGH_SAMPLE, RAYLEIGH_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST132\n";
  cout << "  For the Rayleigh PDF:\n";
  cout << "  RAYLEIGH_MEAN computes the mean;\n";
  cout << "  RAYLEIGH_SAMPLE samples;\n";
  cout << "  RAYLEIGH_VARIANCE computes the variance.\n";

  a = 2.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !rayleigh_check ( a ) )
  {
    cout << "\n";
    cout << "TEST132 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = rayleigh_mean ( a );
  variance = rayleigh_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = rayleigh_sample ( a, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test133 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST133 tests RECIPROCAL_CDF, RECIPROCAL_CDF_INV, RECIPROCAL_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST133\n";
  cout << "  For the Reciprocal PDF:\n";
  cout << "  RECIPROCAL_CDF evaluates the CDF;\n";
  cout << "  RECIPROCAL_CDF_INV inverts the CDF.\n";
  cout << "  RECIPROCAL_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !reciprocal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST133 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = reciprocal_sample ( a, b, &seed );
    pdf = reciprocal_pdf ( x, a, b );
    cdf = reciprocal_cdf ( x, a, b );
    x2 = reciprocal_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test134 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST134 tests RECIPROCAL_MEAN, RECIPROCAL_SAMPLE, RECIPROCAL_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 February 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST134\n";
  cout << "  For the Reciprocal PDF:\n";
  cout << "  RECIPROCAL_MEAN computes the mean;\n";
  cout << "  RECIPROCAL_SAMPLE samples;\n";
  cout << "  RECIPROCAL_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 3.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !reciprocal_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST134 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = reciprocal_mean ( a, b );
  variance = reciprocal_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = reciprocal_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test1341 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1341 checks RIBESL against BESSEL_IX_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 March 2007
//
//  Author:
//
//    John Burkardt
//
{
# define NB_MAX 10

  double alpha;
  double alpha_frac;
  double b[NB_MAX];
  double fx;
  double fx2;
  int ize;
  int n_data;
  int nb;
  int ncalc;
  double x;

  cout << "\n";
  cout << "TEST1341:\n";
  cout << "  RIBESL computes values of Bessel functions\n";
  cout << "  of NONINTEGER order.\n";
  cout << "  BESSEL_IX_VALUES returns selected values of the\n";
  cout << "  Bessel function In for NONINTEGER order.\n";
  cout << "\n";
  cout << "      ALPHA         X             FX                        FX2\n";
  cout << "                                  (table)                   (RIBESL)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    bessel_ix_values ( &n_data, &alpha, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    ize = 1;

    nb = ( int ) alpha + 1;

    if ( NB_MAX < nb )
    {
      cout << "  [Skipping calculation, NB_MAX too small.]\n";
      continue;
    }

    alpha_frac = alpha - ( double ) ( ( int ) alpha );

    ncalc = ribesl ( x, alpha_frac, nb, ize, b );

    fx2 = b[nb-1];

    cout << "  "
         << setw(12) << alpha << "  "
         << setw(12) << x     << "  "
         << setw(24) << fx    << "  "
         << setw(24) << fx2   << "\n";
  }

  return;
# undef NB_MAX
}
//****************************************************************************80

void test1342 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1342 tests RUNS_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int m;
  int n;
  double pdf;
  double pdf_total;
  int r;

  cout << "\n";
  cout << "TEST1342\n";
  cout << "  For the RUNS PDF:\n";
  cout << "  RUNS_PDF evaluates the PDF;\n";
  cout << "\n";
  cout << "  M is the number of symbols of one kind,\n";
  cout << "  N is the number of symbols of the other kind,\n";
  cout << "  R is the number of runs (sequences of one symbol)\n";
  cout << "\n";
  cout << "         M         N         R      PDF\n";
  cout << "\n";

  m = 6;

  for ( n = 0; n <= 9; n++ )
  {
    cout << "\n";
    pdf_total = 0.0;

    for ( r = 1; r <= 2 * i4_min ( m, n ) + 2; r++ )
    {
      pdf = runs_pdf ( m, n, r );

      cout << "  " << setw(8) << m
           << "  " << setw(8) << n
           << "  " << setw(8) << r
           << "  " << setw(14) << pdf << "\n";

      pdf_total = pdf_total + pdf;
    }

    cout << "  " << setw(8) << m
         << "  " << "        "
         << "  " << "        "
         << "  " << setw(14) << pdf_total << "\n";

  }

  return;
}
//****************************************************************************80

void test1344 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1344 tests RUNS_MEAN, RUNS_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  int m;
  double mean;
  int n;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST1344\n";
  cout << "  For the RUNS PDF:\n";
  cout << "  RUNS_MEAN computes the mean;\n";
  cout << "  RUNS_VARIANCE computes the variance\n";

  m = 10;
  n = 5;

  cout << "\n";
  cout << "  PDF parameter M = " << m << "\n";
  cout << "  PDF parameter N = " << n << "\n";

  mean = runs_mean ( m, n );
  variance = runs_variance ( m, n );

  cout << "  PDF mean =        " << mean << "\n";
  cout << "  PDF variance =    " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = runs_sample ( m, n, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test135 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST135 tests SECH_CDF, SECH_CDF_INV, SECH_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST135\n";
  cout << "  For the Sech PDF:\n";
  cout << "  SECH_CDF evaluates the CDF;\n";
  cout << "  SECH_CDF_INV inverts the CDF.\n";
  cout << "  SECH_PDF evaluates the PDF;\n";

  a = 3.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !sech_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST135 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = sech_sample ( a, b, &seed );
    pdf = sech_pdf ( x, a, b );
    cdf = sech_cdf ( x, a, b );
    x2 = sech_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test136 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST136 tests SECH_MEAN, SECH_SAMPLE, SECH_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST136\n";
  cout << "  For the Sech PDF:\n";
  cout << "  SECH_MEAN computes the mean;\n";
  cout << "  SECH_SAMPLE samples;\n";
  cout << "  SECH_VARIANCE computes the variance;\n";

  a = 3.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !sech_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST136 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = sech_mean ( a, b );
  variance = sech_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = sech_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test137 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST137 tests SEMICIRCULAR_CDF, SEMICIRCULAR_CDF_INV, SEMICIRCULAR_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST137\n";
  cout << "  For the Semicircular PDF:\n";
  cout << "  SEMICIRCULAR_CDF evaluates the CDF;\n";
  cout << "  SEMICIRCULAR_CDF_INV inverts the CDF.\n";
  cout << "  SEMICIRCULAR_PDF evaluates the PDF;\n";

  a = 3.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !semicircular_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST137 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = semicircular_sample ( a, b, &seed );
    pdf = semicircular_pdf ( x, a, b );
    cdf = semicircular_cdf ( x, a, b );
    x2 = semicircular_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test138 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST138 tests SEMICIRCULAR_MEAN, SEMICIRCULAR_SAMPLE and SEMICIRCULAR_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST138\n";
  cout << "  For the Semicircular PDF:\n";
  cout << "  SEMICIRCULAR_MEAN computes the mean;\n";
  cout << "  SEMICIRCULAR_SAMPLE samples;\n";
  cout << "  SEMICIRCULAR_VARIANCE computes the variance;\n";

  a = 3.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !semicircular_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST138 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = semicircular_mean ( a, b );
  variance = semicircular_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = semicircular_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test139 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST139 tests STUDENT_CDF and STUDENT_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST139:\n";
  cout << "  STUDENT_CDF evaluates the cumulative density function\n";
  cout << "    for the Student's T PDF.\n";
  cout << "  STUDENT_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "   A      B      C     X       Exact F       STUDENT_CDF(A,B,C,X)\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    student_cdf_values ( &n_data, &c, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    a = 0.0;
    b = 1.0;

    fx2 = student_cdf ( x, a, b, c );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(8)  << c   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test140 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST140 tests STUDENT_CDF, STUDENT_PDF and STUDENT_SAMPLE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;

  cout << "\n";
  cout << "TEST140\n";
  cout << "  For the Student PDF:\n";
  cout << "  STUDENT_CDF evaluates the CDF;\n";
  cout << "  STUDENT_PDF evaluates the PDF;\n";
  cout << "  STUDENT_SAMPLE samples the PDF;\n";

  a = 0.5;
  b = 2.0;
  c = 6.0;

  cout << "\n";
  cout << "  PDF parameter A = " << a   << "\n";
  cout << "  PDF parameter B = " << b   << "\n";
  cout << "  PDF parameter C = " << c   << "\n";

  if ( !student_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST140 - Fatal error!\n";
    cout << "  The parameter values are illegal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = student_sample ( a, b, c, &seed );
    pdf = student_pdf ( x, a, b, c );
    cdf = student_cdf ( x, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test141 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST141 tests STUDENT_MEAN, STUDENT_SAMPLE and STUDENT_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST141\n";
  cout << "  For the Student PDF:\n";
  cout << "  STUDENT_MEAN evaluates the mean;\n";
  cout << "  STUDENT_SAMPLE samples the PDF;\n";
  cout << "  STUDENT_VARIANCE computes the variance;\n";

  a = 0.5;
  b = 2.0;
  c = 6.0;

  cout << "\n";
  cout << "  PDF parameter A = " << a   << "\n";
  cout << "  PDF parameter B = " << b   << "\n";
  cout << "  PDF parameter C = " << c   << "\n";

  if ( !student_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST141 - Fatal error!\n";
    cout << "  The parameter values are illegal.\n";
    return;
  }

  mean = student_mean ( a, b, c );
  variance = student_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = student_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";


  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test142 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST142 tests STUDENT_NONCENTRAL_CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double b;
  double cdf;
  int idf;
  double x;

  cout << "\n";
  cout << "TEST142\n";
  cout << "  For the Noncentral Student PDF:\n";
  cout << "  STUDENT_NONCENTRAL_CDF evaluates the CDF;\n";

  x = 0.50;
  idf = 10;
  b = 1.0;

  cdf = student_noncentral_cdf ( x, idf, b );

  cout << "\n";
  cout << "  PDF argument X =              " << x   << "\n";
  cout << "  PDF parameter IDF =           " << idf << "\n";
  cout << "  PDF parameter B =             " << b   << "\n";
  cout << "  CDF value =                   " << cdf << "\n";

  return;
}
//****************************************************************************80

void test1425 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1425 tests TFN, OWEN_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double h;
  int n_data;
  double t;
  double t2;

  cout << "\n";
  cout << "TEST1425\n";
  cout << "  TFN evaluates Owen's T function;\n";
  cout << "  OWEN_VALUES returns some exact values;\n";

  cout << "\n";
  cout << "      H             A           T(H,A)          Exact\n";
  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    owen_values ( &n_data, &h, &a, &t );

    if ( n_data <= 0 )
    {
      break;
    }

    t2 = tfn ( h, a );

    cout << "  "
         << setw(14) << h  << "  "
         << setw(14) << a  << "  "
         << setw(14) << t2  << "  "
         << setw(14) << t  << "\n";
  }

  return;
}
//****************************************************************************80

void test143 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST143 tests TRIANGLE_CDF, TRIANGLE_CDF_INV, TRIANGLE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST143\n";
  cout << "  For the Triangle PDF:\n";
  cout << "  TRIANGLE_CDF evaluates the CDF;\n";
  cout << "  TRIANGLE_CDF_INV inverts the CDF.\n";
  cout << "  TRIANGLE_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 3.0;
  c = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !triangle_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST143 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = triangle_sample ( a, b, c, &seed );
    pdf = triangle_pdf ( x, a, b, c );
    cdf = triangle_cdf ( x, a, b, c );
    x2 = triangle_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test144 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST144 tests TRIANGLE_MEAN, TRIANGLE_SAMPLE, TRIANGLE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST144\n";
  cout << "  For the Triangle PDF:\n";
  cout << "  TRIANGLE_MEAN computes the mean;\n";
  cout << "  TRIANGLE_SAMPLE samples;\n";
  cout << "  TRIANGLE_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 3.0;
  c = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =        " << a << "\n";
  cout << "  PDF parameter B =        " << b << "\n";
  cout << "  PDF parameter C =        " << c << "\n";

  if ( !triangle_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST144 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = triangle_mean ( a, b, c );
  variance = triangle_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = triangle_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test145 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST145 tests TRIANGULAR_CDF, TRIANGULAR_CDF_INV, TRIANGULAR_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST145\n";
  cout << "  For the Triangular PDF:\n";
  cout << "  TRIANGULAR_CDF evaluates the CDF;\n";
  cout << "  TRIANGULAR_CDF_INV inverts the CDF.\n";
  cout << "  TRIANGULAR_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !triangular_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST145 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = triangular_sample ( a, b, &seed );
    pdf = triangular_pdf ( x, a, b );
    cdf = triangular_cdf ( x, a, b );
    x2 = triangular_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test146 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST146 tests TRIANGULAR_MEAN, TRIANGULAR_SAMPLE, TRIANGULAR_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST146\n";
  cout << "  For the Triangular PDF:\n";
  cout << "  TRIANGULAR_MEAN computes the mean;\n";
  cout << "  TRIANGULAR_SAMPLE samples;\n";
  cout << "  TRIANGULAR_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !triangular_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST146 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = triangular_mean ( a, b );
  variance = triangular_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = triangular_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test147 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST147 tests UNIFORM_01_ORDER_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int n;
  int seed = 123456789;
  double *x;

  cout << "\n";
  cout << "TEST147\n";
  cout << "  For the Uniform 01 Order PDF:\n";
  cout << "  UNIFORM_ORDER_SAMPLE samples.\n";

  n = 10;
  x = uniform_01_order_sample ( n, &seed );

  r8vec_print ( n, x, "  Ordered sample:" );

  delete [] x;

  return;
}
//****************************************************************************80

void test148 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST148 tests UNIFORM_NSPHERE_SAMPLE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  int j;
  int n;
  int seed = 123456789;
  double *x;

  n = 3;

  cout << "\n";
  cout << "TEST148\n";
  cout << "  For the Uniform PDF on the N-Sphere:\n";
  cout << "  UNIFORM_NSPHERE_SAMPLE samples.\n";

  cout << "\n";
  cout << "  Dimension N of sphere = " << n << "\n";
  cout << "\n";
  cout << "  Points on the sphere:\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_nsphere_sample ( n, &seed );
    cout << "  " << setw(6) << i << "  ";
    for ( j = 0; j < n; j++ )
    {
      cout << setw(12) << x[j] << "  ";
    }
    cout << "\n";
    delete [] x;
  }

  return;
}
//****************************************************************************80

void test1485 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1485 tests UNIFORM_01_CDF, UNIFORM_01_CDF_INV, UNIFORM_01_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST1485\n";
  cout << "  For the Uniform 01 PDF:\n";
  cout << "  UNIFORM_01_CDF evaluates the CDF;\n";
  cout << "  UNIFORM_01_CDF_INV inverts the CDF.\n";
  cout << "  UNIFORM_01_PDF evaluates the PDF;\n";

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_01_sample ( &seed );
    pdf = uniform_01_pdf ( x );
    cdf = uniform_01_cdf ( x );
    x2 = uniform_01_cdf_inv ( cdf );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test1486 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1486 tests UNIFORM_01_MEAN, UNIFORM_01_SAMPLE, UNIFORM_01_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST1486\n";
  cout << "  For the Uniform 01 PDF:\n";
  cout << "  UNIFORM_01_MEAN computes the mean;\n";
  cout << "  UNIFORM_01_SAMPLE samples;\n";
  cout << "  UNIFORM_01_VARIANCE computes the variance.\n";

  mean     = uniform_01_mean ( );
  variance = uniform_01_variance ( );

  cout << "\n";
  cout << "  PDF mean =     " << mean << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_01_sample ( &seed );
  }

  mean     = r8vec_mean     ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax     = r8vec_max      ( SAMPLE_NUM, x );
  xmin     = r8vec_min      ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test149 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST149 tests UNIFORM_CDF, UNIFORM_CDF_INV, UNIFORM_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST149\n";
  cout << "  For the Uniform PDF:\n";
  cout << "  UNIFORM_CDF evaluates the CDF;\n";
  cout << "  UNIFORM_CDF_INV inverts the CDF.\n";
  cout << "  UNIFORM_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !uniform_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST149 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_sample ( a, b, &seed );
    pdf = uniform_pdf ( x, a, b );
    cdf = uniform_cdf ( x, a, b );
    x2 = uniform_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test150 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST150 tests UNIFORM_MEAN, UNIFORM_SAMPLE, UNIFORM_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST150\n";
  cout << "  For the Uniform PDF:\n";
  cout << "  UNIFORM_MEAN computes the mean;\n";
  cout << "  UNIFORM_SAMPLE samples;\n";
  cout << "  UNIFORM_VARIANCE computes the variance;\n";

  a = 1.0;
  b = 10.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !uniform_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST150 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = uniform_mean ( a, b );
  variance = uniform_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test151 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST151 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_CDF_INV, UNIFORM_DISCRETE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST151\n";
  cout << "  For the Uniform Discrete PDF:\n";
  cout << "  UNIFORM_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  UNIFORM_DISCRETE_CDF_INV inverts the CDF.\n";
  cout << "  UNIFORM_DISCRETE_PDF evaluates the PDF;\n";

  a = 1;
  b = 6;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !uniform_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST151 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = uniform_discrete_sample ( a, b, &seed );
    pdf = uniform_discrete_pdf ( x, a, b );
    cdf = uniform_discrete_cdf ( x, a, b );
    x2 = uniform_discrete_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test152 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST152 tests UNIFORM_DISCRETE_MEAN, UNIFORM_DISCRETE_SAMPLE, UNIFORM_DISCRETE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  int a;
  int b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST152\n";
  cout << "  For the Uniform Discrete PDF:\n";
  cout << "  UNIFORM_DISCRETE_MEAN computes the mean;\n";
  cout << "  UNIFORM_DISCRETE_SAMPLE samples;\n";
  cout << "  UNIFORM_DISCRETE_VARIANCE computes the variance;\n";

  a = 1;
  b = 6;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !uniform_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST152 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = uniform_discrete_mean ( a, b );
  variance = uniform_discrete_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = uniform_discrete_sample ( a, b, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test153 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST153 tests UNIFORM_DISCRETE_CDF, UNIFORM_DISCRETE_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  int a;
  int b;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST153\n";
  cout << "  For the Uniform Discrete PDF:\n";
  cout << "  UNIFORM_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  UNIFORM_DISCRETE_PDF evaluates the PDF;\n";

  a = 1;
  b = 6;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";
  cout << "  PDF parameter B =             " << b   << "\n";

  if ( !uniform_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST153 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 0; x <= 10; x++ )
  {
    pdf = uniform_discrete_pdf ( x, a, b );
    cdf = uniform_discrete_cdf ( x, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test154 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST154 tests VON_MISES_CDF, VON_MISES_CDF_INV, VON_MISES_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST154\n";
  cout << "  For the Von Mises PDF:\n";
  cout << "  VON_MISES_CDF evaluates the CDF;\n";
  cout << "  VON_MISES_CDF_INV inverts the CDF.\n";
  cout << "  VON_MISES_PDF evaluates the PDF;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !von_mises_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST154 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = von_mises_sample ( a, b, &seed );
    pdf = von_mises_pdf ( x, a, b );
    cdf = von_mises_cdf ( x, a, b );
    x2 = von_mises_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test155 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST155 tests VON_MISES_MEAN, VON_MISES_SAMPLE, VON_MISES_CIRCULAR_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST155\n";
  cout << "  For the Von Mises PDF:\n";
  cout << "  VON_MISES_MEAN computes the mean;\n";
  cout << "  VON_MISES_SAMPLE samples;\n";
  cout << "  VON_MISES_CIRCULAR_VARIANCE computes the circular variance;\n";

  a = 1.0;
  b = 2.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !von_mises_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST155 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = von_mises_mean ( a, b );
  variance = von_mises_circular_variance ( a, b );

  cout << "\n";
  cout << "  PDF mean =              " << mean     << "\n";
  cout << "  PDF circular variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = von_mises_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_circular_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =              " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =              " << mean     << "\n";
  cout << "  Sample circular variance = " << variance << "\n";
  cout << "  Sample maximum =           " << xmax     << "\n";
  cout << "  Sample minimum =           " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test1555 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST1555 tests VON_MISES_CDF, VON_MISES_CDF_VALUES.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double fx;
  double fx2;
  int n_data;
  double x;

  cout << "\n";
  cout << "TEST1555:\n";
  cout << "  VON_MISES_CDF evaluates the von Mises CDF.\n";
  cout << "  VON_MISES_CDF_VALUES returns some exact values.\n";
  cout << "\n";
  cout << "  A is the dominant angle;\n";
  cout << "  B is a measure of spread;\n";
  cout << "  X is the angle;\n";
  cout << "\n";
  cout << "      A     B         X   Exact F     Computed F\n";

  cout << "\n";

  n_data = 0;

  for ( ; ; )
  {
    von_mises_cdf_values ( &n_data, &a, &b, &x, &fx );

    if ( n_data == 0 )
    {
      break;
    }

    fx2 = von_mises_cdf ( x, a, b );

    cout << "  "
         << setw(8)  << a   << "  "
         << setw(8)  << b   << "  "
         << setw(8)  << x   << "  "
         << setw(16) << fx  << "  "
         << setw(16) << fx2 << "\n";
  }

  return;
}
//****************************************************************************80

void test156 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST156 tests WEIBULL_CDF, WEIBULL_CDF_INV, WEIBULL_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double c;
  double cdf;
  int i;
  double pdf;
  int seed = 123456789;
  double x;
  double x2;

  cout << "\n";
  cout << "TEST156\n";
  cout << "  For the Weibull PDF:\n";
  cout << "  WEIBULL_CDF evaluates the CDF;\n";
  cout << "  WEIBULL_CDF_INV inverts the CDF.\n";
  cout << "  WEIBULL_PDF evaluates the PDF;\n";

  a = 2.0;
  b = 3.0;
  c = 4.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !weibull_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST156 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = weibull_sample ( a, b, c, &seed );
    pdf = weibull_pdf ( x, a, b, c );
    cdf = weibull_cdf ( x, a, b, c );
    x2 = weibull_cdf_inv ( cdf, a, b, c );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test157 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST157 tests WEIBULL_MEAN, WEIBULL_SAMPLE, WEIBULL_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  double c;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST157\n";
  cout << "  For the Weibull PDF:\n";
  cout << "  WEIBULL_MEAN computes the mean;\n";
  cout << "  WEIBULL_SAMPLE samples;\n";
  cout << "  WEIBULL_VARIANCE computes the variance.\n";

  a = 2.0;
  b = 3.0;
  c = 4.0;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";
  cout << "  PDF parameter C =      " << c << "\n";

  if ( !weibull_check ( a, b, c ) )
  {
    cout << "\n";
    cout << "TEST157 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = weibull_mean ( a, b, c );
  variance = weibull_variance ( a, b, c );

  cout << "\n";
  cout << "  PDF mean =     " << mean     << "\n";
  cout << "  PDF variance = " << variance << "\n";

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = weibull_sample ( a, b, c, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;

# undef SAMPLE_NUM
}
//****************************************************************************80

void test158 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST158 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_CDF_INV, WEIBULL_DISCRETE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
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
  double pdf;
  int seed = 123456789;
  int x;
  int x2;

  cout << "\n";
  cout << "TEST158\n";
  cout << "  For the Weibull Discrete PDF:\n";
  cout << "  WEIBULL_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  WEIBULL_DISCRETE_CDF_INV inverts the CDF.\n";
  cout << "  WEIBULL_DISCRETE_PDF evaluates the PDF;\n";

  a = 0.5;
  b = 1.5;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !weibull_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST158 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF            CDF_INV\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    x = weibull_discrete_sample ( a, b, &seed );
    pdf = weibull_discrete_pdf ( x, a, b );
    cdf = weibull_discrete_cdf ( x, a, b );
    x2 = weibull_discrete_cdf_inv ( cdf, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "  "
         << setw(12) << x2  << "\n";
  }

  return;
}
//****************************************************************************80

void test159 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST159 tests WEIBULL_DISCRETE_CDF, WEIBULL_DISCRETE_PDF;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double b;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST159\n";
  cout << "  For the Weibull Discrete PDF:\n";
  cout << "  WEIBULL_DISCRETE_CDF evaluates the CDF;\n";
  cout << "  WEIBULL_DISCRETE_PDF evaluates the PDF;\n";

  a = 0.5;
  b = 1.5;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !weibull_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST159 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 0; x <= 10; x++ )
  {
    pdf = weibull_discrete_pdf ( x, a, b );
    cdf = weibull_discrete_cdf ( x, a, b );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test160 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST160 tests WEIBULL_DISCRETE_MEAN, WEIBULL_DISCRETE_SAMPLE, WEIBULL_DISCRETE_VARIANCE;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  double b;
  int i;
  double mean;
  int seed = 123456789;
  double variance;
  double x[SAMPLE_NUM];
  double xmax;
  double xmin;

  cout << "\n";
  cout << "TEST160\n";
  cout << "  For the Weibull Discrete PDF:\n";
  cout << "  WEIBULL_DISCRETE_SAMPLE samples;\n";

  a = 0.5;
  b = 1.5;

  cout << "\n";
  cout << "  PDF parameter A =      " << a << "\n";
  cout << "  PDF parameter B =      " << b << "\n";

  if ( !weibull_discrete_check ( a, b ) )
  {
    cout << "\n";
    cout << "TEST160 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  for ( i = 0; i < SAMPLE_NUM; i++ )
  {
    x[i] = weibull_discrete_sample ( a, b, &seed );
  }

  mean = r8vec_mean ( SAMPLE_NUM, x );
  variance = r8vec_variance ( SAMPLE_NUM, x );
  xmax = r8vec_max ( SAMPLE_NUM, x );
  xmin = r8vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
//****************************************************************************80

void test161 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST161 tests ZIPF_CDF, ZIPF_PDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
  double a;
  double cdf;
  double pdf;
  int x;

  cout << "\n";
  cout << "TEST161\n";
  cout << "  For the ZIPF PDF:\n";
  cout << "  ZIPF_CDF evaluates the CDF;\n";
  cout << "  ZIPF_PDF evaluates the PDF;\n";

  a = 2.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a   << "\n";

  if ( !zipf_check ( a ) )
  {
    cout << "\n";
    cout << "TEST161 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  cout << "\n";
  cout << "       X            PDF           CDF\n";
  cout << "\n";

  for ( x = 1; x <= 20; x++ )
  {
    pdf = zipf_pdf ( x, a );
    cdf = zipf_cdf ( x, a );

    cout << "  "
         << setw(12) << x   << "  "
         << setw(12) << pdf << "  "
         << setw(12) << cdf << "\n";
  }

  return;
}
//****************************************************************************80

void test162 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST162 tests ZIPF_MEAN, ZIPF_SAMPLE, ZIPF_VARIANCE.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2007
//
//  Author:
//
//    John Burkardt
//
{
# define SAMPLE_NUM 1000

  double a;
  int i;
  int j;
  double mean;
  int seed = 123456789;
  double variance;
  int x[SAMPLE_NUM];
  int xmax;
  int xmin;

  cout << "\n";
  cout << "TEST162\n";
  cout << "  For the Zipf PDF:\n";
  cout << "  ZIPF_MEAN computes the mean;\n";
  cout << "  ZIPF_SAMPLE samples;\n";
  cout << "  ZIPF_VARIANCE computes the variance.\n";

  a = 4.0E+00;

  cout << "\n";
  cout << "  PDF parameter A =             " << a        << "\n";

  if ( !zipf_check ( a ) )
  {
    cout << "\n";
    cout << "TEST162 - Fatal error!\n";
    cout << "  The parameters are not legal.\n";
    return;
  }

  mean = zipf_mean ( a );
  variance = zipf_variance ( a );

  cout << "  PDF mean =                    " << mean     << "\n";
  cout << "  PDF variance =                " << variance << "\n";

  for ( j = 0; j < SAMPLE_NUM; j++ )
  {
    x[j] = zipf_sample ( a, &seed );
  }

  mean = i4vec_mean ( SAMPLE_NUM, x );
  variance = i4vec_variance ( SAMPLE_NUM, x );
  xmax = i4vec_max ( SAMPLE_NUM, x );
  xmin = i4vec_min ( SAMPLE_NUM, x );

  cout << "\n";
  cout << "  Sample size =     " << SAMPLE_NUM  << "\n";
  cout << "  Sample mean =     " << mean     << "\n";
  cout << "  Sample variance = " << variance << "\n";
  cout << "  Sample maximum =  " << xmax     << "\n";
  cout << "  Sample minimum =  " << xmin     << "\n";

  return;
# undef SAMPLE_NUM
}
