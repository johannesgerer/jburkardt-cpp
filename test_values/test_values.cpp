# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <complex>

using namespace std;

# include "test_values.hpp"

//****************************************************************************80

void abram0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM0_VALUES returns some values of the Abramowitz0 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      ABRAM0(X) = integral ( 0 <= T < +oo ) exp ( -T * T - X / T ) dT
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.87377726306985360531E+00,
     0.84721859650456925922E+00,
     0.77288934483988301615E+00,
     0.59684345853450151603E+00,
     0.29871735283675888392E+00,
     0.15004596450516388138E+00,
     0.11114662419157955096E+00,
     0.83909567153151897766E-01,
     0.56552321717943417515E-01,
     0.49876496603033790206E-01,
     0.44100889219762791328E-01,
     0.19738535180254062496E-01,
     0.86193088287161479900E-02,
     0.40224788162540127227E-02,
     0.19718658458164884826E-02,
     0.10045868340133538505E-02,
     0.15726917263304498649E-03,
     0.10352666912350263437E-04,
     0.91229759190956745069E-06,
     0.25628287737952698742E-09 };
  static double x_vec[N_MAX] = {
     0.0019531250E+00,
     0.0078125000E+00,
     0.0312500000E+00,
     0.1250000000E+00,
     0.5000000000E+00,
     1.0000000000E+00,
     1.2500000000E+00,
     1.5000000000E+00,
     1.8750000000E+00,
     2.0000000000E+00,
     2.1250000000E+00,
     3.0000000000E+00,
     4.0000000000E+00,
     5.0000000000E+00,
     6.0000000000E+00,
     7.0000000000E+00,
     10.0000000000E+00,
     15.0000000000E+00,
     20.0000000000E+00,
     40.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void abram1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM1_VALUES returns some values of the Abramowitz1 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      ABRAM1(x) = integral ( 0 <= t < oo ) t * exp ( -t^2 - x / t ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.49828219848799921792E+00,
     0.49324391773047288556E+00,
     0.47431612784691234649E+00,
     0.41095983258760410149E+00,
     0.25317617388227035867E+00,
     0.14656338138597777543E+00,
     0.11421547056018366587E+00,
     0.90026307383483764795E-01,
     0.64088214170742303375E-01,
     0.57446614314166191085E-01,
     0.51581624564800730959E-01,
     0.25263719555776416016E-01,
     0.11930803330196594536E-01,
     0.59270542280915272465E-02,
     0.30609215358017829567E-02,
     0.16307382136979552833E-02,
     0.28371851916959455295E-03,
     0.21122150121323238154E-04,
     0.20344578892601627337E-05,
     0.71116517236209642290E-09 };

  static double x_vec[N_MAX] = {
     0.0019531250E+00,
     0.0078125000E+00,
     0.0312500000E+00,
     0.1250000000E+00,
     0.5000000000E+00,
     1.0000000000E+00,
     1.2500000000E+00,
     1.5000000000E+00,
     1.8750000000E+00,
     2.0000000000E+00,
     2.1250000000E+00,
     3.0000000000E+00,
     4.0000000000E+00,
     5.0000000000E+00,
     6.0000000000E+00,
     7.0000000000E+00,
     10.0000000000E+00,
     15.0000000000E+00,
     20.0000000000E+00,
     40.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void abram2_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ABRAM2_VALUES returns some values of the Abramowitz2 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      ABRAM2(x) = Integral ( 0 <= t < +oo ) t^2 * exp( -t^2 - x / t ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.44213858162107913430E+00,
     0.43923379545684026308E+00,
     0.42789857297092602234E+00,
     0.38652825661854504406E+00,
     0.26538204413231368110E+00,
     0.16848734838334595000E+00,
     0.13609200032513227112E+00,
     0.11070330027727917352E+00,
     0.82126019995530382267E-01,
     0.74538781999594581763E-01,
     0.67732034377612811390E-01,
     0.35641808698811851022E-01,
     0.17956589956618269083E-01,
     0.94058737143575370625E-02,
     0.50809356204299213556E-02,
     0.28149565414209719359E-02,
     0.53808696422559303431E-03,
     0.44821756380146327259E-04,
     0.46890678427324100410E-05,
     0.20161544850996420504E-08 };

  static double x_vec[N_MAX] = {
     0.0019531250E+00,
     0.0078125000E+00,
     0.0312500000E+00,
     0.1250000000E+00,
     0.5000000000E+00,
     1.0000000000E+00,
     1.2500000000E+00,
     1.5000000000E+00,
     1.8750000000E+00,
     2.0000000000E+00,
     2.1250000000E+00,
     3.0000000000E+00,
     4.0000000000E+00,
     5.0000000000E+00,
     6.0000000000E+00,
     7.0000000000E+00,
     10.0000000000E+00,
     15.0000000000E+00,
     20.0000000000E+00,
     40.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void agm_values ( int &n_data, double &a, double &b, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    AGM_VALUES returns some values of the AGM.
//
//  Discussion:
//
//    The AGM is defined for nonnegative A and B.
//
//    The AGM of numbers A and B is defined by setting
//
//      A(0) = A,
//      B(0) = B
//
//      A(N+1) = ( A(N) + B(N) ) / 2
//      B(N+1) = sqrt ( A(N) * B(N) )
//
//    The two sequences both converge to AGM(A,B).
//
//    In Mathematica, the AGM can be evaluated by
//
//      ArithmeticGeometricMean [ a, b ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the argument ofs the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double a_vec[N_MAX] = {
        22.0,
        83.0,
        42.0,
        26.0,
         4.0,
         6.0,
        40.0,
        80.0,
        90.0,
         9.0,
        53.0,
         1.0,
         1.0,
         1.0,
         1.5 };
  static double b_vec[N_MAX] = {
        96.0,
        56.0,
         7.0,
        11.0,
        63.0,
        45.0,
        75.0,
         0.0,
        35.0,
         1.0,
        53.0,
         2.0,
         4.0,
         8.0,
         8.0 };
  static double fx_vec[N_MAX] = {
        52.274641198704240049,
        68.836530059858524345,
        20.659301196734009322,
        17.696854873743648823,
        23.867049721753300163,
        20.717015982805991662,
        56.127842255616681863,
         0.000000000000000000,
        59.269565081229636528,
        3.9362355036495554780,
        53.000000000000000000,
        1.4567910310469068692,
        2.2430285802876025701,
        3.6157561775973627487,
        4.0816924080221632670 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_ai_values ( int &n_data, double &x, double &ai )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_VALUES returns some values of the Airy Ai(x) function.
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryAi[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &AI, the value of the Airy AI function.
//
{
# define N_MAX 11

  static double ai_vec[N_MAX] = {
     0.3550280538878172E+00,
     0.3292031299435381E+00,
     0.3037031542863820E+00,
     0.2788064819550049E+00,
     0.2547423542956763E+00,
     0.2316936064808335E+00,
     0.2098000616663795E+00,
     0.1891624003981501E+00,
     0.1698463174443649E+00,
     0.1518868036405444E+00,
     0.1352924163128814E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    ai = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    ai = ai_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_ai_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_INT_VALUES returns some values of the integral of the Airy function.
//
//  Discussion:
//
//    The function is defined by:
//
//      AIRY_AI_INT(x) = Integral ( 0 <= t <= x ) Ai(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     -0.75228838916610124300E+00,
     -0.57348350185854889466E+00,
     -0.76569840313421291743E+00,
     -0.65181015505382467421E+00,
     -0.55881974894471876922E+00,
     -0.56902352870716815309E+00,
     -0.47800749642926168100E+00,
     -0.46567398346706861416E+00,
     -0.96783140945618013679E-01,
     -0.34683049857035607494E-03,
      0.34658366917927930790E-03,
      0.27657581846051227124E-02,
      0.14595330491185717833E+00,
      0.23631734191710977960E+00,
      0.33289264538612212697E+00,
      0.33318759129779422976E+00,
      0.33332945170523851439E+00,
      0.33333331724248357420E+00,
      0.33333333329916901594E+00,
      0.33333333333329380187E+00 };

  static double x_vec[N_MAX] = {
     -12.0000000000E+00,
     -11.0000000000E+00,
     -10.0000000000E+00,
      -9.5000000000E+00,
      -9.0000000000E+00,
      -6.5000000000E+00,
      -4.0000000000E+00,
      -1.0000000000E+00,
      -0.2500000000E+00,
      -0.0009765625E+00,
       0.0009765625E+00,
       0.0078125000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       4.0000000000E+00,
       4.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_ai_prime_values ( int &n_data, double &x, double &aip )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_AI_PRIME_VALUES returns some values of the Airy function Ai'(x).
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryAiPrime[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &AIP, the derivative of the Airy AI function.
//
{
# define N_MAX 11

  static double aip_vec[N_MAX] = {
     -0.2588194037928068E+00,
     -0.2571304219075862E+00,
     -0.2524054702856195E+00,
     -0.2451463642190548E+00,
     -0.2358320344192082E+00,
     -0.2249105326646839E+00,
     -0.2127932593891585E+00,
     -0.1998511915822805E+00,
     -0.1864128638072717E+00,
     -0.1727638434616347E+00,
     -0.1591474412967932E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    aip = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    aip = aip_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_bi_values ( int &n_data, double &x, double &bi )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_VALUES returns some values of the Airy Bi(x) function.
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryBi[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &BI, the value of the Airy BI function.
//
{
# define N_MAX 11

  static double bi_vec[N_MAX] = {
     0.6149266274460007E+00,
     0.6598616901941892E+00,
     0.7054642029186612E+00,
     0.7524855850873156E+00,
     0.8017730000135972E+00,
     0.8542770431031555E+00,
     0.9110633416949405E+00,
     0.9733286558781659E+00,
     0.1042422171231561E+01,
     0.1119872813134447E+01,
     0.1207423594952871E+01 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    bi = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    bi = bi_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_bi_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_INT_VALUES returns some values of the integral of the Airy function.
//
//  Discussion:
//
//    The function is defined by:
//
//      AIRY_BI_INT(x) = Integral ( 0 <= t <= x ) Bi(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.17660819031554631869E-01,
     -0.15040424806140020451E-01,
      0.14756446293227661920E-01,
     -0.11847304264848446271E+00,
     -0.64916741266165856037E-01,
      0.97260832464381044540E-01,
      0.50760058495287539119E-01,
     -0.37300500963429492179E+00,
     -0.13962988442666578531E+00,
     -0.12001735266723296160E-02,
      0.12018836117890354598E-02,
      0.36533846550952011043E+00,
      0.87276911673800812196E+00,
      0.48219475263803429675E+02,
      0.44006525804904178439E+06,
      0.17608153976228301458E+07,
      0.73779211705220007228E+07,
      0.14780980310740671617E+09,
      0.97037614223613433849E+11,
      0.11632737638809878460E+15 };

  static double x_vec[N_MAX] = {
     -12.0000000000E+00,
     -10.0000000000E+00,
      -8.0000000000E+00,
      -7.5000000000E+00,
      -7.0000000000E+00,
      -6.5000000000E+00,
      -4.0000000000E+00,
      -1.0000000000E+00,
      -0.2500000000E+00,
      -0.0019531250E+00,
       0.0019531250E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       4.0000000000E+00,
       8.0000000000E+00,
       8.5000000000E+00,
       9.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00,
      14.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_bi_prime_values ( int &n_data, double &x, double &bip )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_BI_PRIME_VALUES returns some values of the Airy function Bi'(x).
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryBiPrime[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &BIP, the derivative of the Airy BI function.
//
{
# define N_MAX 11

  static double bip_vec[N_MAX] = {
     0.4482883573538264E+00,
     0.4515126311496465E+00,
     0.4617892843621509E+00,
     0.4800490287524480E+00,
     0.5072816760506224E+00,
     0.5445725641405923E+00,
     0.5931444786342857E+00,
     0.6544059191721400E+00,
     0.7300069016152518E+00,
     0.8219038903072090E+00,
     0.9324359333927756E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    bip = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    bip = bip_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_cai_values ( int &n_data, complex <double> &x,
  complex <double> &cai )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_CAI_VALUES returns some values of the Airy Ai(x) for complex argument.
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryAi[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, complex <double> &X, the argument of the function.
//
//    Output, complex <double> &CAI, the value of the Airy AI function.
//
{
# define N_MAX 10

  static complex <double> cai_vec[N_MAX] = {
   complex <double> ( 0.1352924163128814, + 0.0000000000000000  ),
   complex <double> ( 0.1433824486882056, - 0.1092193342707378  ),
   complex <double> ( 0.2215404472324631, - 0.2588711788891803  ),
   complex <double> ( 0.4763929771766866, - 0.3036484220291284  ),
   complex <double> ( 0.5983692170633874, - 0.08154602160771214 ),
   complex <double> ( 0.5355608832923521, + 0.00000000000000000 ),
   complex <double> ( 0.5983692170633874, + 0.08154602160771214 ),
   complex <double> ( 0.4763929771766866, + 0.3036484220291284  ),
   complex <double> ( 0.2215404472324631, + 0.2588711788891803  ),
   complex <double> ( 0.1433824486882056, + 0.1092193342707378  ) };

  static complex <double> x_vec[N_MAX] = {
   complex <double> (  1.0000000000000000, + 0.0000000000000000 ),
   complex <double> (  0.8090169943749474, + 0.5877852522924731 ),
   complex <double> (  0.3090169943749474, + 0.9510565162951536 ),
   complex <double> ( -0.3090169943749474, + 0.9510565162951536 ),
   complex <double> ( -0.8090169943749474, + 0.5877852522924731 ),
   complex <double> ( -1.0000000000000000, + 0.0000000000000000 ),
   complex <double> ( -0.8090169943749474, - 0.5877852522924731 ),
   complex <double> ( -0.3090169943749474, - 0.9510565162951536 ),
   complex <double> (  0.3090169943749474, - 0.9510565162951536 ),
   complex <double> (  0.8090169943749474, - 0.5877852522924731 ) };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    cai = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    cai = cai_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_cbi_values ( int &n_data, complex <double> &x,
  complex <double> &cbi )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_CBI_VALUES returns some values of the Airy Bi(x) for complex argument.
//
//  Discussion:
//
//    The Airy functions Ai(X) and Bi(X) are a pair of linearly independent
//    solutions of the differential equation:
//
//      W'' - X * W = 0;
//
//    In Mathematica, the function can be evaluated by:
//
//      AiryAi[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, complex <double> &X, the argument of the function.
//
//    Output, complex <double> &CBI, the value of the Airy BI function.
//
{
# define N_MAX 10

  static complex <double> cbi_vec[N_MAX] = {
    complex <double> ( 1.207423594952871,  + 0.0000000000000000 ),
    complex <double> ( 0.9127160108293936, + 0.3800456133135556 ),
    complex <double> ( 0.6824453575635721, + 0.3343047153635002 ),
    complex <double> ( 0.5726265660086474, + 0.3988641086982559 ),
    complex <double> ( 0.2511841251049547, + 0.3401447690712719 ),
    complex <double> ( 0.1039973894969446, + 0.0000000000000000 ),
    complex <double> ( 0.2511841251049547, - 0.3401447690712719 ),
    complex <double> ( 0.5726265660086474, - 0.3988641086982559 ),
    complex <double> ( 0.6824453575635721, - 0.3343047153635002 ),
    complex <double> ( 0.9127160108293936, - 0.3800456133135556 ) };

  static complex <double> x_vec[N_MAX] = {
    complex <double> (  1.0000000000000000, + 0.0000000000000000 ),
    complex <double> (  0.8090169943749474, + 0.5877852522924731 ),
    complex <double> (  0.3090169943749474, + 0.9510565162951536 ),
    complex <double> ( -0.3090169943749474, + 0.9510565162951536 ),
    complex <double> ( -0.8090169943749474, + 0.5877852522924731 ),
    complex <double> ( -1.0000000000000000, + 0.0000000000000000 ),
    complex <double> ( -0.8090169943749474, - 0.5877852522924731 ),
    complex <double> ( -0.3090169943749474, - 0.9510565162951536 ),
    complex <double> (  0.3090169943749474, - 0.9510565162951536 ),
    complex <double> (  0.8090169943749474, - 0.5877852522924731 ) };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    cbi = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    cbi = cbi_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_gi_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_GI_VALUES returns some values of the Airy Gi function.
//
//  Discussion:
//
//    The function is defined by:
//
//      AIRY_GI(x) = Integral ( 0 <= t < +oo ) sin ( x*t+t^3/3) dt / pi
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.20468308070040542435E+00,
      0.18374662832557904078E+00,
     -0.11667221729601528265E+00,
      0.31466934902729557596E+00,
     -0.37089040722426257729E+00,
     -0.25293059772424019694E+00,
      0.28967410658692701936E+00,
     -0.34644836492634090590E+00,
      0.28076035913873049496E+00,
      0.21814994508094865815E+00,
      0.20526679000810503329E+00,
      0.22123695363784773258E+00,
      0.23521843981043793760E+00,
      0.82834303363768729338E-01,
      0.45757385490989281893E-01,
      0.44150012014605159922E-01,
      0.39951133719508907541E-01,
      0.35467706833949671483E-01,
      0.31896005100679587981E-01,
      0.26556892713512410405E-01 };

  static double x_vec[N_MAX] = {
      -0.0019531250E+00,
      -0.1250000000E+00,
      -1.0000000000E+00,
      -4.0000000000E+00,
      -8.0000000000E+00,
      -8.2500000000E+00,
      -9.0000000000E+00,
     -10.0000000000E+00,
     -11.0000000000E+00,
     -13.0000000000E+00,
       0.0019531250E+00,
       0.1250000000E+00,
       1.0000000000E+00,
       4.0000000000E+00,
       7.0000000000E+00,
       7.2500000000E+00,
       8.0000000000E+00,
       9.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void airy_hi_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    AIRY_HI_VALUES returns some values of the Airy Hi function.
//
//  Discussion:
//
//    The function is defined by:
//
//      AIRY_HI(x) = Integral ( 0 <= t < +oo ) exp(x*t-t^3/3) dt / pi
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.40936798278458884024E+00,
     0.37495291608048868619E+00,
     0.22066960679295989454E+00,
     0.77565356679703713590E-01,
     0.39638826473124717315E-01,
     0.38450072575004151871E-01,
     0.35273216868317898556E-01,
     0.31768535282502272742E-01,
     0.28894408288051391369E-01,
     0.24463284011678541180E-01,
     0.41053540139998941517E+00,
     0.44993502381204990817E+00,
     0.97220515514243332184E+00,
     0.83764237105104371193E+02,
     0.80327744952044756016E+05,
     0.15514138847749108298E+06,
     0.11995859641733262114E+07,
     0.21472868855967642259E+08,
     0.45564115351632913590E+09,
     0.32980722582904761929E+12 };

  static double x_vec[N_MAX] = {
      -0.0019531250E+00,
      -0.1250000000E+00,
      -1.0000000000E+00,
      -4.0000000000E+00,
      -8.0000000000E+00,
      -8.2500000000E+00,
      -9.0000000000E+00,
     -10.0000000000E+00,
     -11.0000000000E+00,
     -13.0000000000E+00,
       0.0019531250E+00,
       0.1250000000E+00,
       1.0000000000E+00,
       4.0000000000E+00,
       7.0000000000E+00,
       7.2500000000E+00,
       8.0000000000E+00,
       9.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arccos_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCCOS_VALUES returns some values of the arc cosine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcCos[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
    1.6709637479564564156,
    1.5707963267948966192,
    1.4706289056333368229,
    1.3694384060045658278,
    1.2661036727794991113,
    1.1592794807274085998,
    1.0471975511965977462,
    0.92729521800161223243,
    0.79539883018414355549,
    0.64350110879328438680,
    0.45102681179626243254,
    0.00000000000000000000 };

  static double x_vec[N_MAX] = {
    -0.1,
     0.0,
     0.1,
     0.2,
     0.3,
     0.4,
     0.5,
     0.6,
     0.7,
     0.8,
     0.9,
     1.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arccosh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCCOSH_VALUES returns some values of the hyperbolic arc cosine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcCosh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double fx_vec[N_MAX] = {
    0.0000000000000000000,
    0.14130376948564857735,
    0.44356825438511518913,
    0.62236250371477866781,
    0.75643291085695958624,
    0.86701472649056510395,
    0.96242365011920689500,
    1.3169578969248167086,
    1.7627471740390860505,
    1.8115262724608531070,
    2.0634370688955605467,
    2.2924316695611776878,
    2.9932228461263808979,
    5.2982923656104845907,
    7.6009022095419886114 };

  static double x_vec[N_MAX] = {
       1.0,
       1.01,
       1.1,
       1.2,
       1.3,
       1.4,
       1.5,
       2.0,
       3.0,
       3.1415926535897932385,
       4.0,
       5.0,
      10.0,
     100.0,
    1000.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arcsin_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSIN_VALUES returns some values of the arc sine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcSin[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
    -0.10016742116155979635,
     0.00000000000000000000,
     0.10016742116155979635,
     0.20135792079033079146,
     0.30469265401539750797,
     0.41151684606748801938,
     0.52359877559829887308,
     0.64350110879328438680,
     0.77539749661075306374,
     0.92729521800161223243,
     1.1197695149986341867,
     1.5707963267948966192 };

  static double x_vec[N_MAX] = {
    -0.1,
     0.0,
     0.1,
     0.2,
     0.3,
     0.4,
     0.5,
     0.6,
     0.7,
     0.8,
     0.9,
     1.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arcsinh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCSINH_VALUES returns some values of the hyperbolic arc sine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcSinh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
    -2.3124383412727526203,
    -0.88137358701954302523,
     0.00000000000000000000,
     0.099834078899207563327,
     0.19869011034924140647,
     0.29567304756342243910,
     0.39003531977071527608,
     0.48121182505960344750,
     0.56882489873224753010,
     0.65266656608235578681,
     0.73266825604541086415,
     0.80886693565278246251,
     0.88137358701954302523,
     1.4436354751788103425,
     1.8184464592320668235,
     2.0947125472611012942,
     2.3124383412727526203,
     2.9982229502979697388,
     5.2983423656105887574,
     7.6009027095419886115 };

  static double x_vec[N_MAX] = {
       -5.0,
       -1.0,
        0.0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        2.0,
        3.0,
        4.0,
        5.0,
       10.0,
      100.0,
     1000.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arctan_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTAN_VALUES returns some values of the arc tangent function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcTan[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double f_vec[N_MAX] = {
    0.00000000000000000000,
    0.24497866312686415417,
    0.32175055439664219340,
    0.46364760900080611621,
    0.78539816339744830962,
    1.1071487177940905030,
    1.2490457723982544258,
    1.3258176636680324651,
    1.3734007669450158609,
    1.4711276743037345919,
    1.5208379310729538578 };

  static double x_vec[N_MAX] = {
    0.00000000000000000000,
    0.25000000000000000000,
    0.33333333333333333333,
    0.50000000000000000000,
    1.0000000000000000000,
    2.0000000000000000000,
    3.0000000000000000000,
    4.0000000000000000000,
    5.0000000000000000000,
    10.000000000000000000,
    20.000000000000000000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arctan_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTAN_INT_VALUES returns some values of the inverse tangent integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      ARCTAN_INT(x) = Integral ( 0 <= t <= x ) arctan ( t ) / t dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double f_vec[N_MAX] = {
      0.19531241721588483191E-02,
     -0.39062433772980711281E-02,
      0.78124470192576499535E-02,
      0.15624576181996527280E-01,
     -0.31246610349485401551E-01,
      0.62472911335014397321E-01,
      0.12478419717389654039E+00,
     -0.24830175098230686908E+00,
      0.48722235829452235711E+00,
      0.91596559417721901505E+00,
      0.12749694484943800618E+01,
     -0.15760154034463234224E+01,
      0.24258878412859089996E+01,
      0.33911633326292997361E+01,
      0.44176450919422186583E+01,
     -0.47556713749547247774E+01,
      0.50961912150934111303E+01,
      0.53759175735714876256E+01,
     -0.61649904785027487422E+01,
      0.72437843013083534973E+01 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
      -0.0039062500E+00,
       0.0078125000E+00,
       0.0156250000E+00,
      -0.0312500000E+00,
       0.0625000000E+00,
       0.1250000000E+00,
      -0.2500000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
      -2.0000000000E+00,
       4.0000000000E+00,
       8.0000000000E+00,
      16.0000000000E+00,
     -20.0000000000E+00,
      25.0000000000E+00,
      30.0000000000E+00,
     -50.0000000000E+00,
     100.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arctan2_values ( int &n_data, double &x, double &y, double &fxy )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTAN2_VALUES: arc tangent function of two arguments.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcTan[x,y]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, &Y, the arguments of the function.
//
//    Output, double &FXY, the value of the function.
//
{
# define N_MAX 19

  static double f_vec[N_MAX] = {
   -1.5707963267948966192,
   -1.0471975511965977462,
   -0.52359877559829887308,
    0.00000000000000000000,
    0.52359877559829887308,
    1.0471975511965977462,
    1.5707963267948966192,
    2.0943951023931954923,
    2.6179938779914943654,
    3.1415926535897932385,
   -2.6179938779914943654,
   -2.0943951023931954923,
   -1.5707963267948966192,
   -1.0471975511965977462,
   -0.52359877559829887308,
    0.00000000000000000000,
    0.52359877559829887308,
    1.0471975511965977462,
    1.5707963267948966192 };

  static double x_vec[N_MAX] = {
    0.00000000000000000000,
    0.50000000000000000000,
    0.86602540378443864676,
    1.00000000000000000000,
    0.86602540378443864676,
    0.50000000000000000000,
    0.00000000000000000000,
   -0.50000000000000000000,
   -0.86602540378443864676,
   -1.00000000000000000000,
   -0.86602540378443864676,
   -0.50000000000000000000,
    0.00000000000000000000,
    0.50000000000000000000,
    0.86602540378443864676,
    1.00000000000000000000,
    0.86602540378443864676,
    0.50000000000000000000,
    0.00000000000000000000 };

  static double y_vec[N_MAX] = {
   -1.00000000000000000000,
   -0.86602540378443864676,
   -0.50000000000000000000,
    0.00000000000000000000,
    0.50000000000000000000,
    0.86602540378443864676,
    1.00000000000000000000,
    0.86602540378443864676,
    0.50000000000000000000,
    0.00000000000000000000,
   -0.50000000000000000000,
   -0.86602540378443864676,
   -1.00000000000000000000,
   -0.86602540378443864676,
   -0.50000000000000000000,
    0.00000000000000000000,
    0.50000000000000000000,
    0.86602540378443864676,
    1.00000000000000000000
  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    y = 0.0;
    fxy = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    y = y_vec[n_data-1];
    fxy = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void arctanh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ARCTANH_VALUES returns some values of the hyperbolic arc tangent function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ArcTanh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double fx_vec[N_MAX] = {
    -0.54930614433405484570,
     0.00000000000000000000,
     0.0010000003333335333335,
     0.10033534773107558064,
     0.20273255405408219099,
     0.30951960420311171547,
     0.42364893019360180686,
     0.54930614433405484570,
     0.69314718055994530942,
     0.86730052769405319443,
     1.0986122886681096914,
     1.4722194895832202300,
     2.6466524123622461977,
     3.8002011672502000318,
     7.2543286192620472067 };

  static double x_vec[N_MAX] = {
    -0.5,
     0.0,
     0.001,
     0.1,
     0.2,
     0.3,
     0.4,
     0.5,
     0.6,
     0.7,
     0.8,
     0.9,
     0.99,
     0.999,
     0.999999 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bei0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BEI0_VALUES returns some values of the Kelvin BEI function of order NU = 0.
//
//  Discussion:
//
//    The function is defined by:
//
//      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
//
//    where J(NU,X) is the J Bessel function.
//
//    In Mathematica, BEI(NU,X) can be defined by:
//
//      Im [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
    0.0000000000000000,
    0.06249321838219946,
    0.2495660400366597,
    0.5575600623030867,
    0.9722916273066612,
    1.457182044159804,
   1.937586785266043,
   2.283249966853915,
   2.292690322699300,
   1.686017203632139,
   0.1160343815502004 };
  static double x_vec[N_MAX] = {
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bei1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BEI1_VALUES returns some values of the Kelvin BEI function of order NU = 1.
//
//  Discussion:
//
//    The function is defined by:
//
//      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
//
//    where J(NU,X) is the J Bessel function.
//
//    In Mathematica, BEI(NU,X) can be defined by:
//
//      Im [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
    0.0000000000000000,
    0.1711951797170153,
    0.3075566313755366,
    0.3678649890020899,
    0.2997754370020335,
    0.03866844396595048,
   -0.4874541770160708,
   -1.344042373111174,
   -2.563821688561078,
   -4.105685408400878,
   -5.797907901792625 };
  static double x_vec[N_MAX] = {
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bell_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    BELL_VALUES returns some values of the Bell numbers.
//
//  Discussion:
//
//    The Bell number B(N) is the number of restricted growth functions on N.
//
//    Note that the Stirling numbers of the second kind, S^m_n, count the
//    number of partitions of N objects into M classes, and so it is
//    true that
//
//      B(N) = S^1_N + S^2_N + ... + S^N_N.
//
//    The Bell numbers were named for Eric Temple Bell.
//
//    In Mathematica, the function can be evaluated by
//
//      Sum[StirlingS2[n,m],{m,1,n}]
//
//  Definition:
//
//    The Bell number B(N) is defined as the number of partitions (of
//    any size) of a set of N distinguishable objects.
//
//    A partition of a set is a division of the objects of the set into
//    subsets.
//
//  Examples:
//
//    There are 15 partitions of a set of 4 objects:
//
//      (1234),
//      (123) (4),
//      (124) (3),
//      (12) (34),
//      (12) (3) (4),
//      (134) (2),
//      (13) (24),
//      (13) (2) (4),
//      (14) (23),
//      (1) (234),
//      (1) (23) (4),
//      (14) (2) (3),
//      (1) (24) (3),
//      (1) (2) (34),
//      (1) (2) (3) (4).
//
//    and so B(4) = 15.
//
//  First values:
//
//     N         B(N)
//     0           1
//     1           1
//     2           2
//     3           5
//     4          15
//     5          52
//     6         203
//     7         877
//     8        4140
//     9       21147
//    10      115975
//
//  Recursion:
//
//    B(I) = sum ( 1 <= J <=I ) Binomial ( I-1, J-1 ) * B(I-J)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Bell number.
//
//    Output, int &C, the value of the Bell number.
//
{
# define N_MAX 11

  static int c_vec[N_MAX] = {
    1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ber0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BER0_VALUES returns some values of the Kelvin BER function of order NU = 0.
//
//  Discussion:
//
//    The function is defined by:
//
//      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
//
//    where J(NU,X) is the J Bessel function.
//
//    In Mathematica, BER(NU,X) can be defined by:
//
//      Re [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
    1.0000000000000000,
    0.9990234639908383,
    0.9843817812130869,
    0.9210721835462558,
    0.7517341827138082,
    0.3999684171295313,
   -0.2213802495986939,
   -1.193598179589928,
   -2.563416557258580,
   -4.299086551599756,
   -6.230082478666358 };
  static double x_vec[N_MAX] = {
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ber1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BER1_VALUES returns some values of the Kelvin BER function of order NU = 1.
//
//  Discussion:
//
//    The function is defined by:
//
//      BER(NU,X) + i * BEI(NU,X) = exp(NU*Pi*I) * J(NU,X*exp(-PI*I/4))
//
//    where J(NU,X) is the J Bessel function.
//
//    In Mathematica, BER(NU,X) can be defined by:
//
//      Re [ Exp [ NU * Pi * I ] * BesselJ [ NU, X * Exp[ -Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
     0.0000000000000000,
    -0.1822431237551121,
    -0.3958682610197114,
    -0.6648654179597691,
    -0.9970776519264285,
    -1.373096897645111,
   -1.732644221128481,
   -1.959644131289749,
   -1.869248459031899,
   -1.202821631480086,
    0.3597766667766728 };
  static double x_vec[N_MAX] = {
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bernoulli_number_values ( int &n_data, int &n, double &c )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_NUMBER_VALUES returns some values of the Bernoulli numbers.
//
//  Discussion:
//
//    The Bernoulli numbers are rational.
//
//    If we define the sum of the M-th powers of the first N integers as:
//
//      SIGMA(M,N) = sum ( 0 <= I <= N ) I**M
//
//    and let C(I,J) be the combinatorial coefficient:
//
//      C(I,J) = I! / ( ( I - J )! * J! )
//
//    then the Bernoulli numbers B(J) satisfy:
//
//      SIGMA(M,N) = 1/(M+1) * sum ( 0 <= J <= M ) C(M+1,J) B(J) * (N+1)**(M+1-J)
//
//    In Mathematica, the function can be evaluated by:
//
//      BernoulliB[n]
//
//  First values:
//
//   B0  1                   =         1.00000000000
//   B1 -1/2                 =        -0.50000000000
//   B2  1/6                 =         1.66666666666
//   B3  0                   =         0
//   B4 -1/30                =        -0.03333333333
//   B5  0                   =         0
//   B6  1/42                =         0.02380952380
//   B7  0                   =         0
//   B8 -1/30                =        -0.03333333333
//   B9  0                   =         0
//  B10  5/66                =         0.07575757575
//  B11  0                   =         0
//  B12 -691/2730            =        -0.25311355311
//  B13  0                   =         0
//  B14  7/6                 =         1.16666666666
//  B15  0                   =         0
//  B16 -3617/510            =        -7.09215686274
//  B17  0                   =         0
//  B18  43867/798           =        54.97117794486
//  B19  0                   =         0
//  B20 -174611/330          =      -529.12424242424
//  B21  0                   =         0
//  B22  854,513/138         =      6192.123
//  B23  0                   =         0
//  B24 -236364091/2730      =    -86580.257
//  B25  0                   =         0
//  B26  8553103/6           =   1425517.16666
//  B27  0                   =         0
//  B28 -23749461029/870     = -27298231.0678
//  B29  0                   =         0
//  B30  8615841276005/14322 = 601580873.901
//
//  Recursion:
//
//    With C(N+1,K) denoting the standard binomial coefficient,
//
//    B(0) = 1.0
//    B(N) = - ( sum ( 0 <= K < N ) C(N+1,K) * B(K) ) / C(N+1,N)
//
//  Special Values:
//
//    Except for B(1), all Bernoulli numbers of odd index are 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Bernoulli number.
//
//    Output, double &C, the value of the Bernoulli number.
//
{
# define N_MAX 10

  static double c_vec[N_MAX] = {
      0.1000000000000000E+01,
     -0.5000000000000000E+00,
      0.1666666666666667E+00,
      0.0000000000000000E+00,
     -0.3333333333333333E-01,
     -0.2380952380952380E-01,
     -0.3333333333333333E-01,
      0.7575757575757575E-01,
     -0.5291242424242424E+03,
      0.6015808739006424E+09 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,  4, 6,  8, 10, 20, 30 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0.0E+00;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bernoulli_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BERNOULLI_POLY_VALUES returns some values of the Bernoulli polynomials.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BernoulliB[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Bernoulli polynomial.
//
//    Output, double &X, the argument of the Bernoulli polynomial.
//
//    Output, double &FX, the value of the Bernoulli polynomial.
//
{
# define N_MAX 27

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
     -0.3000000000000000E+00,
      0.6666666666666667E-02,
      0.4800000000000000E-01,
     -0.7733333333333333E-02,
     -0.2368000000000000E-01,
      0.6913523809523810E-02,
      0.2490880000000000E-01,
     -0.1014997333333333E-01,
     -0.4527820800000000E-01,
      0.2332631815757576E-01,
     -0.3125000000000000E+00,
     -0.1142400000000000E+00,
     -0.0176800000000000E+00,
      0.0156800000000000E+00,
      0.0147400000000000E+00,
      0.0000000000000000E+00,
     -0.1524000000000000E-01,
     -0.2368000000000000E-01,
     -0.2282000000000000E-01,
     -0.1376000000000000E-01,
      0.0000000000000000E+01,
      0.1376000000000000E-01,
      0.2282000000000000E-01,
      0.2368000000000000E-01,
      0.1524000000000000E-01,
      0.0000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0,
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
     -0.5E+00,
     -0.4E+00,
     -0.3E+00,
     -0.2E+00,
     -0.1E+00,
      0.0E+00,
      0.1E+00,
      0.2E+00,
      0.3E+00,
      0.4E+00,
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bernstein_poly_values ( int &n_data, int &n, int &k, double &x, double &b )

//****************************************************************************80
//
//  Purpose:
//
//    BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
//
//  Discussion:
//
//    The Bernstein polynomials are assumed to be based on [0,1].
//
//    The formula for the Bernstein polynomials is
//
//      B(N,I)(X) = [N!/(I//(N-I)!)] * (1-X)^(N-I) * X^I
//
//    In Mathematica, the function can be evaluated by:
//
//      Binomial[n,i] * (1-x)^(n-i) * x^i
//
//  First values:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =               X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)   * X
//    B(2,2)(X) =               X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
//  Special values:
//
//    B(N,I)(X) has a unique maximum value at X = I/N.
//
//    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
//
//    B(N,I)(1/2) = C(N,K) / 2^N
//
//    For a fixed X and N, the polynomials add up to 1:
//
//      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the degree of the polynomial.
//
//    Output, int &K, the index of the polynomial.
//
//    Output, double &X, the argument of the polynomial.
//
//    Output, double &B, the value of the polynomial B(N,K)(X).
//
{
# define N_MAX 15

  static double b_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.7500000000000000E+00,
     0.2500000000000000E+00,
     0.5625000000000000E+00,
     0.3750000000000000E+00,
     0.6250000000000000E-01,
     0.4218750000000000E+00,
     0.4218750000000000E+00,
     0.1406250000000000E+00,
     0.1562500000000000E-01,
     0.3164062500000000E+00,
     0.4218750000000000E+00,
     0.2109375000000000E+00,
     0.4687500000000000E-01,
     0.3906250000000000E-02 };

  static int k_vec[N_MAX] = {
    0,
    0, 1,
    0, 1, 2,
    0, 1, 2, 3,
    0, 1, 2, 3, 4 };

  static int n_vec[N_MAX] = {
    0,
    1, 1,
    2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4, 4 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    k = 0;
    x = 0.0;
    b = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    k = k_vec[n_data-1];
    x = x_vec[n_data-1];
    b = b_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_i0_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_INT_VALUES returns some values of the Bessel I0 integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      I0_INT(x) = Integral ( 0 <= t <= x ) I0(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.19531256208818052282E-02,
     -0.39062549670565734544E-02,
      0.62520348032546565850E-01,
      0.12516285581366971819E+00,
     -0.51051480879740303760E+00,
      0.10865210970235898158E+01,
      0.27750019054282535299E+01,
     -0.13775208868039716639E+02,
      0.46424372058106108576E+03,
      0.64111867658021584522E+07,
     -0.10414860803175857953E+08,
      0.44758598913855743089E+08,
     -0.11852985311558287888E+09,
      0.31430078220715992752E+09,
     -0.83440212900794309620E+09,
      0.22175367579074298261E+10,
      0.58991731842803636487E+10,
     -0.41857073244691522147E+11,
      0.79553885818472357663E+12,
      0.15089715082719201025E+17 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
      -0.0039062500E+00,
       0.0625000000E+00,
       0.1250000000E+00,
      -0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
      -4.0000000000E+00,
       8.0000000000E+00,
      18.0000000000E+00,
     -18.5000000000E+00,
      20.0000000000E+00,
     -21.0000000000E+00,
      22.0000000000E+00,
     -23.0000000000E+00,
      24.0000000000E+00,
      25.0000000000E+00,
     -27.0000000000E+00,
      30.0000000000E+00,
      40.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_i0_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_SPHERICAL_VALUES returns some values of the Spherical Bessel function i0.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselI[1/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    1.001667500198440E+00,
    1.006680012705470E+00,
    1.026880814507039E+00,
    1.061089303580402E+00,
    1.110132477734529E+00,
    1.175201193643801E+00,
    1.257884462843477E+00,
    1.360215358179667E+00,
    1.484729970750144E+00,
    1.634541271164267E+00,
    1.813430203923509E+00,
    2.025956895698133E+00,
    2.277595505698373E+00,
    2.574897010920645E+00,
    2.925685126512827E+00,
    3.339291642469967E+00,
    3.826838748926716E+00,
    4.401577467270101E+00,
    5.079293155726485E+00,
    5.878791279137455E+00,
    6.822479299281938E+00 };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_i0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I0_VALUES returns some values of the I0 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    The modified Bessel function I0(Z) corresponds to N = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[0,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1010025027795146E+01,
     0.1040401782229341E+01,
     0.1092045364317340E+01,
     0.1166514922869803E+01,
     0.1266065877752008E+01,
     0.1393725584134064E+01,
     0.1553395099731217E+01,
     0.1749980639738909E+01,
     0.1989559356618051E+01,
     0.2279585302336067E+01,
     0.3289839144050123E+01,
     0.4880792585865024E+01,
     0.7378203432225480E+01,
     0.1130192195213633E+02,
     0.1748117185560928E+02,
     0.2723987182360445E+02,
     0.6723440697647798E+02,
     0.4275641157218048E+03,
     0.2815716628466254E+04 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.20E+00,
     0.40E+00,
     0.60E+00,
     0.80E+00,
     0.10E+01,
     0.12E+01,
     0.14E+01,
     0.16E+01,
     0.18E+01,
     0.20E+01,
     0.25E+01,
     0.30E+01,
     0.35E+01,
     0.40E+01,
     0.45E+01,
     0.50E+01,
     0.60E+01,
     0.80E+01,
     0.10E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_i1_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1_SPHERICAL_VALUES returns some values of the Spherical Bessel function i1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselI[3/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    0.03336667857363341E+00,
    0.06693371456802954E+00,
    0.1354788933285401E+00,
    0.2072931911031093E+00,
    0.2841280857128948E+00,
    0.3678794411714423E+00,
    0.4606425870674146E+00,
    0.5647736480096238E+00,
    0.6829590627779635E+00,
    0.8182955028627777E+00,
    0.9743827435800610E+00,
    1.155432469636406E+00,
    1.366396525527973E+00,
    1.613118767572064E+00,
    1.902515460838681E+00,
    2.242790117769266E+00,
    2.643689828630357E+00,
    3.116811526884873E+00,
    3.675968313148932E+00,
    4.337627987747642E+00,
    5.121438384183637E+00 };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_i1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_I1_VALUES returns some values of the I1 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.1005008340281251E+00,
     0.2040267557335706E+00,
     0.3137040256049221E+00,
     0.4328648026206398E+00,
     0.5651591039924850E+00,
     0.7146779415526431E+00,
     0.8860919814143274E+00,
     0.1084810635129880E+01,
     0.1317167230391899E+01,
     0.1590636854637329E+01,
     0.2516716245288698E+01,
     0.3953370217402609E+01,
     0.6205834922258365E+01,
     0.9759465153704450E+01,
     0.1538922275373592E+02,
     0.2433564214245053E+02,
     0.6134193677764024E+02,
     0.3998731367825601E+03,
     0.2670988303701255E+04 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.20E+00,
     0.40E+00,
     0.60E+00,
     0.80E+00,
     0.10E+01,
     0.12E+01,
     0.14E+01,
     0.16E+01,
     0.18E+01,
     0.20E+01,
     0.25E+01,
     0.30E+01,
     0.35E+01,
     0.40E+01,
     0.45E+01,
     0.50E+01,
     0.60E+01,
     0.80E+01,
     0.10E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_in_values ( int &n_data, int &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_IN_VALUES returns some values of the In Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
     0.5016687513894678E-02,
     0.1357476697670383E+00,
     0.6889484476987382E+00,
     0.1276466147819164E+01,
     0.2245212440929951E+01,
     0.1750561496662424E+02,
     0.2281518967726004E+04,
     0.3931278522104076E+08,
     0.2216842492433190E-01,
     0.2127399592398527E+00,
     0.1033115016915114E+02,
     0.1758380716610853E+04,
     0.2677764138883941E+21,
     0.2714631559569719E-03,
     0.9825679323131702E-02,
     0.2157974547322546E+01,
     0.7771882864032600E+03,
     0.2278548307911282E+21,
     0.2752948039836874E-09,
     0.3016963879350684E-06,
     0.4580044419176051E-02,
     0.2189170616372337E+02,
     0.1071597159477637E+21,
     0.3966835985819020E-24,
     0.4310560576109548E-18,
     0.5024239357971806E-10,
     0.1250799735644948E-03,
     0.5442008402752998E+19 };

  static int nu_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  2,  2,  2,
     3,  3,  3,  3,
     3,  5,  5,  5,
     5,  5, 10, 10,
    10, 10, 10, 20,
    20, 20, 20, 20 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_ix_values ( int &n_data, double &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_IX_VALUES returns some values of the Ix Bessel function.
//
//  Discussion:
//
//    This set of data considers the less common case in which the
//    index of the Bessel function In is actually not an integer.
//    We may suggest this case by occasionally replacing the symbol
//    "In" by "Ix".
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselI[n,x]
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
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
    0.3592084175833614E+00,
    0.9376748882454876E+00,
    2.046236863089055E+00,
    3.053093538196718E+00,
    4.614822903407601E+00,
    26.47754749755907E+00,
    2778.784603874571E+00,
    4.327974627242893E+07,
    0.2935253263474798E+00,
    1.099473188633110E+00,
    21.18444226479414E+00,
    2500.906154942118E+00,
    2.866653715931464E+20,
    0.05709890920304825E+00,
    0.3970270801393905E+00,
    13.76688213868258E+00,
    2028.512757391936E+00,
    2.753157630035402E+20,
    0.4139416015642352E+00,
    1.340196758982897E+00,
    22.85715510364670E+00,
    2593.006763432002E+00,
    2.886630075077766E+20,
    0.03590910483251082E+00,
    0.2931108636266483E+00,
    11.99397010023068E+00,
    1894.575731562383E+00,
    2.716911375760483E+20 };

  static double nu_vec[N_MAX] = {
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_j0_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_INT_VALUES returns some values of the Bessel J0 integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      J0_INT(x) = Integral ( 0 <= t <= x ) J0(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.97656242238978822427E-03,
      0.39062450329491108875E-02,
     -0.62479657927917933620E-01,
      0.12483733492120479139E+00,
     -0.48968050664604505505E+00,
      0.91973041008976023931E+00,
     -0.14257702931970265690E+01,
      0.10247341594606064818E+01,
     -0.12107468348304501655E+01,
      0.11008652032736190799E+01,
     -0.10060334829904124192E+01,
      0.81330572662485953519E+00,
     -0.10583788214211277585E+01,
      0.87101492116545875169E+00,
     -0.88424908882547488420E+00,
      0.11257761503599914603E+01,
     -0.90141212258183461184E+00,
      0.91441344369647797803E+00,
     -0.94482281938334394886E+00,
      0.92266255696016607257E+00 };

  static double x_vec[N_MAX] = {
       0.0009765625E+00,
       0.0039062500E+00,
      -0.0625000000E+00,
       0.1250000000E+00,
      -0.5000000000E+00,
       1.0000000000E+00,
      -2.0000000000E+00,
       4.0000000000E+00,
      -8.0000000000E+00,
      16.0000000000E+00,
     -16.5000000000E+00,
      18.0000000000E+00,
     -20.0000000000E+00,
      25.0000000000E+00,
     -30.0000000000E+00,
      40.0000000000E+00,
     -50.0000000000E+00,
      75.0000000000E+00,
     -80.0000000000E+00,
     100.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_j0_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_SPHERICAL_VALUES returns some values of the Spherical Bessel function j0.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselJ[1/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
      0.9983341664682815E+00,
      0.9933466539753061E+00,
      0.9735458557716262E+00,
      0.9410707889917256E+00,
      0.8966951136244035E+00,
      0.8414709848078965E+00,
      0.7766992383060220E+00,
      0.7038926642774716E+00,
      0.6247335019009407E+00,
      0.5410264615989973E+00,
      0.4546487134128408E+00,
      0.3674983653725410E+00,
      0.2814429918963129E+00,
      0.1982697583928709E+00,
      0.1196386250556803E+00,
      0.4704000268662241E-01,
     -0.1824191982111872E-01,
     -0.7515914765495039E-01,
     -0.1229223453596812E+00,
     -0.1610152344586103E+00,
     -0.1892006238269821E+00 };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_j0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J0_VALUES returns some values of the J0 Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselJ[0,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     -0.1775967713143383E+00,
     -0.3971498098638474E+00,
     -0.2600519549019334E+00,
      0.2238907791412357E+00,
      0.7651976865579666E+00,
      0.1000000000000000E+01,
      0.7651976865579666E+00,
      0.2238907791412357E+00,
     -0.2600519549019334E+00,
     -0.3971498098638474E+00,
     -0.1775967713143383E+00,
      0.1506452572509969E+00,
      0.3000792705195556E+00,
      0.1716508071375539E+00,
     -0.9033361118287613E-01,
     -0.2459357644513483E+00,
     -0.1711903004071961E+00,
      0.4768931079683354E-01,
      0.2069261023770678E+00,
      0.1710734761104587E+00,
     -0.1422447282678077E-01 };

  static double x_vec[N_MAX] = {
     -5.0E+00,
     -4.0E+00,
     -3.0E+00,
     -2.0E+00,
     -1.0E+00,
      0.0E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
     11.0E+00,
     12.0E+00,
     13.0E+00,
     14.0E+00,
     15.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_j1_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J1_SPHERICAL_VALUES returns some values of the Spherical Bessel function j1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselJ[3/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.3330001190255757E-01,
     0.6640038067032223E-01,
     0.1312121544218529E+00,
     0.1928919568034122E+00,
     0.2499855053465475E+00,
     0.3011686789397568E+00,
     0.3452845698577903E+00,
     0.3813753724123076E+00,
     0.4087081401263934E+00,
     0.4267936423844913E+00,
     0.4353977749799916E+00,
     0.4345452193763121E+00,
     0.4245152947656493E+00,
     0.4058301968314685E+00,
     0.3792360591872637E+00,
     0.3456774997623560E+00,
     0.3062665174917607E+00,
     0.2622467779189737E+00,
     0.2149544641595738E+00,
     0.1657769677515280E+00,
     0.1161107492591575E+00  };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_j1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_J1_VALUES returns some values of the J1 Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselJ[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
      0.3275791375914652E+00,
      0.6604332802354914E-01,
     -0.3390589585259365E+00,
     -0.5767248077568734E+00,
     -0.4400505857449335E+00,
      0.0000000000000000E+00,
      0.4400505857449335E+00,
      0.5767248077568734E+00,
      0.3390589585259365E+00,
     -0.6604332802354914E-01,
     -0.3275791375914652E+00,
     -0.2766838581275656E+00,
     -0.4682823482345833E-02,
      0.2346363468539146E+00,
      0.2453117865733253E+00,
      0.4347274616886144E-01,
     -0.1767852989567215E+00,
     -0.2234471044906276E+00,
     -0.7031805212177837E-01,
      0.1333751546987933E+00,
      0.2051040386135228E+00 };

  static double x_vec[N_MAX] = {
     -5.0E+00,
     -4.0E+00,
     -3.0E+00,
     -2.0E+00,
     -1.0E+00,
      0.0E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
     11.0E+00,
     12.0E+00,
     13.0E+00,
     14.0E+00,
     15.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_jn_values ( int &n_data, int &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_JN_VALUES returns some values of the Jn Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselJ[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 April 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.1149034849319005E+00,
      0.3528340286156377E+00,
      0.4656511627775222E-01,
      0.2546303136851206E+00,
     -0.5971280079425882E-01,
      0.2497577302112344E-03,
      0.7039629755871685E-02,
      0.2611405461201701E+00,
     -0.2340615281867936E+00,
     -0.8140024769656964E-01,
      0.2630615123687453E-09,
      0.2515386282716737E-06,
      0.1467802647310474E-02,
      0.2074861066333589E+00,
     -0.1138478491494694E+00,
      0.3873503008524658E-24,
      0.3918972805090754E-18,
      0.2770330052128942E-10,
      0.1151336924781340E-04,
     -0.1167043527595797E+00 };

  static int nu_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  5,  5,  5,
     5,  5, 10, 10,
    10, 10, 10, 20,
    20, 20, 20, 20 };

  static double x_vec[N_MAX] = {
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_jx_values ( int &n_data, double &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_JX_VALUES returns some values of the Jx Bessel function.
//
//  Discussion:
//
//    This set of data considers the less common case in which the
//    index of the Bessel function Jn is actually not an integer.
//    We may suggest this case by occasionally replacing the symbol
//    "Jn" by "Jx".
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselJ[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
       0.3544507442114011E+00,
        0.6713967071418031E+00,
        0.5130161365618278E+00,
        0.3020049060623657E+00,
        0.06500818287737578E+00,
       -0.3421679847981618E+00,
       -0.1372637357550505E+00,
        0.1628807638550299E+00,
        0.2402978391234270E+00,
        0.4912937786871623E+00,
       -0.1696513061447408E+00,
        0.1979824927558931E+00,
       -0.1094768729883180E+00,
        0.04949681022847794E+00,
        0.2239245314689158E+00,
        0.2403772011113174E+00,
        0.1966584835818184E+00,
        0.02303721950962553E+00,
        0.3314145508558904E+00,
        0.5461734240402840E+00,
       -0.2616584152094124E+00,
        0.1296035513791289E+00,
       -0.1117432171933552E+00,
        0.03142623570527935E+00,
        0.1717922192746527E+00,
        0.3126634069544786E+00,
        0.1340289119304364E+00,
        0.06235967135106445E+00 };
  static double nu_vec[N_MAX] = {
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_k0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K0_VALUES returns some values of the K0 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    The modified Bessel function K0(Z) corresponds to N = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselK[0,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.2427069024702017E+01,
     0.1752703855528146E+01,
     0.1114529134524434E+01,
     0.7775220919047293E+00,
     0.5653471052658957E+00,
     0.4210244382407083E+00,
     0.3185082202865936E+00,
     0.2436550611815419E+00,
     0.1879547519693323E+00,
     0.1459314004898280E+00,
     0.1138938727495334E+00,
     0.6234755320036619E-01,
     0.3473950438627925E-01,
     0.1959889717036849E-01,
     0.1115967608585302E-01,
     0.6399857243233975E-02,
     0.3691098334042594E-02,
     0.1243994328013123E-02,
     0.1464707052228154E-03,
     0.1778006231616765E-04 };

  static double x_vec[N_MAX] = {
      0.1E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00,
      5.0E+00,
      6.0E+00,
      8.0E+00,
     10.0E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_k0_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K0_INT_VALUES returns some values of the Bessel K0 integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      K0_INT(x) = Integral ( 0 <= t <= x ) K0(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.78587929563466784589E-02,
     0.26019991617330578111E-01,
     0.24311842237541167904E+00,
     0.39999633750480508861E+00,
     0.92710252093114907345E+00,
     0.12425098486237782662E+01,
     0.14736757343168286825E+01,
     0.15606495706051741364E+01,
     0.15673873907283660493E+01,
     0.15696345532693743714E+01,
     0.15701153443250786355E+01,
     0.15706574852894436220E+01,
     0.15707793116159788598E+01,
     0.15707942066465767196E+01,
     0.15707962315469192247E+01,
     0.15707963262340149876E+01,
     0.15707963267948756308E+01,
     0.15707963267948966192E+01,
     0.15707963267948966192E+01,
     0.15707963267948966192E+01  };

  static double x_vec[N_MAX] = {
       0.0009765625E+00,
       0.0039062500E+00,
       0.0625000000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
       4.0000000000E+00,
       5.0000000000E+00,
       6.0000000000E+00,
       6.5000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00,
      80.0000000000E+00,
     100.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_k1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_K1_VALUES returns some values of the K1 Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    The modified Bessel function K1(Z) corresponds to N = 1.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselK[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.9853844780870606E+01,
     0.4775972543220472E+01,
     0.2184354424732687E+01,
     0.1302834939763502E+01,
     0.8617816344721803E+00,
     0.6019072301972346E+00,
     0.4345923910607150E+00,
     0.3208359022298758E+00,
     0.2406339113576119E+00,
     0.1826230998017470E+00,
     0.1398658818165224E+00,
     0.7389081634774706E-01,
     0.4015643112819418E-01,
     0.2223939292592383E-01,
     0.1248349888726843E-01,
     0.7078094908968090E-02,
     0.4044613445452164E-02,
     0.1343919717735509E-02,
     0.1553692118050011E-03,
     0.1864877345382558E-04 };

  static double x_vec[N_MAX] = {
      0.1E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00,
      5.0E+00,
      6.0E+00,
      8.0E+00,
     10.0E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_kn_values ( int &n_data, int &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_KN_VALUES returns some values of the Kn Bessel function.
//
//  Discussion:
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 * W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselK[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
     0.4951242928773287E+02,
     0.1624838898635177E+01,
     0.2537597545660559E+00,
     0.1214602062785638E+00,
     0.6151045847174204E-01,
     0.5308943712223460E-02,
     0.2150981700693277E-04,
     0.6329543612292228E-09,
     0.7101262824737945E+01,
     0.6473853909486342E+00,
     0.8291768415230932E-02,
     0.2725270025659869E-04,
     0.3727936773826211E-22,
     0.3609605896012407E+03,
     0.9431049100596467E+01,
     0.3270627371203186E-01,
     0.5754184998531228E-04,
     0.4367182254100986E-22,
     0.1807132899010295E+09,
     0.1624824039795591E+06,
     0.9758562829177810E+01,
     0.1614255300390670E-02,
     0.9150988209987996E-22,
     0.6294369360424535E+23,
     0.5770856852700241E+17,
     0.4827000520621485E+09,
     0.1787442782077055E+03,
     0.1706148379722035E-20 };

  static int nu_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  2,  2,  2,
     3,  3,  3,  3,
     3,  5,  5,  5,
     5,  5, 10, 10,
    10, 10, 10, 20,
    20, 20, 20, 20 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_kx_values ( int &n_data, double &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_KX_VALUES returns some values of the Kx Bessel function.
//
//  Discussion:
//
//    This set of data considers the less common case in which the
//    index of the Bessel function Kn is actually not an integer.
//    We may suggest this case by occasionally replacing the symbol
//    "Kn" by "Kx".
//
//    The modified Bessel functions In(Z) and Kn(Z) are solutions of
//    the differential equation
//
//      Z^2 W'' + Z * W' - ( Z^2 + N^2 ) * W = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselK[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
       2.294489339798475E+00,
       0.4610685044478946E+00,
       0.1199377719680614E+00,
       0.06506594315400999E+00,
       0.03602598513176459E+00,
       0.003776613374642883E+00,
       0.00001799347809370518E+00,
       5.776373974707445E-10,
       0.9221370088957891E+00,
       0.1799066579520922E+00,
       0.004531936049571459E+00,
       0.00001979282590307570E+00,
       3.486992497366216E-23,
       3.227479531135262E+00,
       0.3897977588961997E+00,
       0.006495775004385758E+00,
       0.00002393132586462789E+00,
       3.627839645299048E-23,
       0.7311451879202114E+00,
       0.1567475478393932E+00,
       0.004257389528177461E+00,
       0.00001915541065869563E+00,
       3.463337593569306E-23,
       4.731184839919541E+00,
       0.4976876225514758E+00,
       0.007300864610941163E+00,
       0.00002546421294106458E+00,
       3.675275677913656E-23 };

  static double nu_vec[N_MAX] = {
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_y0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_VALUES returns some values of the Y0 Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselY[0,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     -0.1534238651350367E+01,
      0.8825696421567696E-01,
      0.5103756726497451E+00,
      0.3768500100127904E+00,
     -0.1694073932506499E-01,
     -0.3085176252490338E+00,
     -0.2881946839815792E+00,
     -0.2594974396720926E-01,
      0.2235214893875662E+00,
      0.2499366982850247E+00,
      0.5567116728359939E-01,
     -0.1688473238920795E+00,
     -0.2252373126343614E+00,
     -0.7820786452787591E-01,
      0.1271925685821837E+00,
      0.2054642960389183E+00 };


  static double x_vec[N_MAX] = {
      0.1E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
     11.0E+00,
     12.0E+00,
     13.0E+00,
     14.0E+00,
     15.0E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_y0_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_INT_VALUES returns some values of the Bessel Y0 integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      Y0_INT(x) = Integral ( 0 <= t <= x ) Y0(t) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     -0.91442642860172110926E-02,
     -0.29682047390397591290E-01,
     -0.25391431276585388961E+00,
     -0.56179545591464028187E+00,
     -0.63706937660742309754E+00,
     -0.28219285008510084123E+00,
      0.38366964785312561103E+00,
     -0.12595061285798929390E+00,
      0.24129031832266684828E+00,
      0.17138069757627037938E+00,
      0.18958142627134083732E+00,
      0.17203846136449706946E+00,
     -0.16821597677215029611E+00,
     -0.93607927351428988679E-01,
      0.88229711948036648408E-01,
     -0.89324662736274161841E-02,
     -0.54814071000063488284E-01,
     -0.94958246003466381588E-01,
     -0.19598064853404969850E-01,
     -0.83084772357154773468E-02 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0078125000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
       4.0000000000E+00,
       6.0000000000E+00,
      10.0000000000E+00,
      16.0000000000E+00,
      16.2500000000E+00,
      17.0000000000E+00,
      20.0000000000E+00,
      25.0000000000E+00,
      30.0000000000E+00,
      40.0000000000E+00,
      50.0000000000E+00,
      70.0000000000E+00,
     100.0000000000E+00,
     125.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_y0_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y0_SPHERICAL_VALUES returns some values of the Spherical Bessel function y0.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselY[1/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     -0.9950041652780258E+01,
     -0.4900332889206208E+01,
     -0.2302652485007213E+01,
     -0.1375559358182797E+01,
     -0.8708833866839568E+00,
     -0.5403023058681397E+00,
     -0.3019647953972280E+00,
     -0.1214051020716007E+00,
      0.1824970143830545E-01,
      0.1262233859406039E+00,
      0.2080734182735712E+00,
      0.2675005078433390E+00,
      0.3072473814755190E+00,
      0.3295725974495951E+00,
      0.3365079788102351E+00,
      0.3299974988668152E+00,
      0.3119671174358603E+00,
      0.2843524095821944E+00,
      0.2490995600928186E+00,
      0.2081493978722149E+00,
      0.1634109052159030E+00 };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_y1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y1_VALUES returns some values of the Y1 Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselY[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     -0.6458951094702027E+01,
     -0.7812128213002887E+00,
     -0.1070324315409375E+00,
      0.3246744247918000E+00,
      0.3979257105571000E+00,
      0.1478631433912268E+00,
     -0.1750103443003983E+00,
     -0.3026672370241849E+00,
     -0.1580604617312475E+00,
      0.1043145751967159E+00,
      0.2490154242069539E+00,
      0.1637055374149429E+00,
     -0.5709921826089652E-01,
     -0.2100814084206935E+00,
     -0.1666448418561723E+00,
      0.2107362803687351E-01 };

  static double x_vec[N_MAX] = {
      0.1E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
     11.0E+00,
     12.0E+00,
     13.0E+00,
     14.0E+00,
     15.0E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_y1_spherical_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_Y1_SPHERICAL_VALUES returns some values of the Spherical Bessel function y1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi/(2*x)] * BesselY[3/2,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     -0.1004987506942709E+03,
     -0.2549501110000635E+02,
     -0.6730177068289658E+01,
     -0.3233669719296388E+01,
     -0.1985299346979349E+01,
     -0.1381773290676036E+01,
     -0.1028336567803712E+01,
     -0.7906105943286149E+00,
     -0.6133274385019998E+00,
     -0.4709023582986618E+00,
     -0.3506120042760553E+00,
     -0.2459072254437506E+00,
     -0.1534232496148467E+00,
     -0.7151106706610352E-01,
      0.5427959479750482E-03,
      0.6295916360231598E-01,
      0.1157316440198251E+00,
      0.1587922092967723E+00,
      0.1921166676076864E+00,
      0.2157913917934037E+00,
      0.2300533501309578E+00 };

  static double x_vec[N_MAX] = {
     0.1E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00,
     3.2E+00,
     3.4E+00,
     3.6E+00,
     3.8E+00,
     4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_yn_values ( int &n_data, int &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_YN_VALUES returns some values of the Yn Bessel function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselY[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
//
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     -0.1650682606816254E+01,
     -0.6174081041906827E+00,
      0.3676628826055245E+00,
     -0.5868082442208615E-02,
      0.9579316872759649E-01,
     -0.2604058666258122E+03,
     -0.9935989128481975E+01,
     -0.4536948224911019E+00,
      0.1354030476893623E+00,
     -0.7854841391308165E-01,
     -0.1216180142786892E+09,
     -0.1291845422080393E+06,
     -0.2512911009561010E+02,
     -0.3598141521834027E+00,
      0.5723897182053514E-02,
     -0.4113970314835505E+23,
     -0.4081651388998367E+17,
     -0.5933965296914321E+09,
     -0.1597483848269626E+04,
      0.1644263394811578E-01 };

  static int nu_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  5,  5,  5,
     5,  5, 10, 10,
    10, 10, 10, 20,
    20, 20, 20, 20 };

  static double x_vec[N_MAX] = {
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bessel_yx_values ( int &n_data, double &nu, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BESSEL_YX_VALUES returns some values of the Yx Bessel function.
//
//  Discussion:
//
//    This set of data considers the less common case in which the
//    index of the Bessel function Yn is actually not an integer.
//    We may suggest this case by occasionally replacing the symbol
//    "Yn" by "Yx".
//
//    In Mathematica, the function can be evaluated by:
//
//      BesselY[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &NU, the order of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
       -1.748560416961876E+00,
       -0.4310988680183761E+00,
        0.2347857104062485E+00,
        0.4042783022390569E+00,
        0.4560488207946332E+00,
       -0.1012177091851084E+00,
        0.2117088663313982E+00,
       -0.07280690478506185E+00,
       -1.102495575160179E+00,
       -0.3956232813587035E+00,
        0.3219244429611401E+00,
        0.1584346223881903E+00,
        0.02742813676191382E+00,
       -2.876387857462161E+00,
       -0.8282206324443037E+00,
        0.2943723749617925E+00,
       -0.1641784796149411E+00,
        0.1105304445562544E+00,
       -0.9319659251969881E+00,
       -0.2609445010948933E+00,
        0.2492796362185881E+00,
        0.2174410301416733E+00,
       -0.01578576650557229E+00,
       -4.023453301501028E+00,
       -0.9588998694752389E+00,
        0.2264260361047367E+00,
       -0.2193617736566760E+00,
        0.09413988344515077E+00 };

  static double nu_vec[N_MAX] = {
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    0.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    1.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    2.50E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    1.25E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00,
    2.75E+00 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      1.0E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      5.0E+00,
     10.0E+00,
     20.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00,
      1.0E+00,
      2.0E+00,
      5.0E+00,
     10.0E+00,
     50.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nu = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    nu = nu_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void beta_cdf_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_CDF_VALUES returns some values of the Beta CDF.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
//                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    The incomplete Beta function is also sometimes called the
//    "modified" Beta function, or the "normalized" Beta function
//    or the Beta CDF (cumulative density function.
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//    The function can also be evaluated by using the Statistics package:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 42

  static double a_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.5E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     30.0E+00,
     30.0E+00,
     40.0E+00,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01 };

  static double b_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.0E+00,
      0.5E+00,
      5.0E+00,
      5.0E+00,
     10.0E+00,
      5.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01 };

  static double fx_vec[N_MAX] = {
     0.6376856085851985E-01,
     0.2048327646991335E+00,
     0.1000000000000000E+01,
     0.0000000000000000E+00,
     0.5012562893380045E-02,
     0.5131670194948620E-01,
     0.2928932188134525E+00,
     0.5000000000000000E+00,
     0.2800000000000000E-01,
     0.1040000000000000E+00,
     0.2160000000000000E+00,
     0.3520000000000000E+00,
     0.5000000000000000E+00,
     0.6480000000000000E+00,
     0.7840000000000000E+00,
     0.8960000000000000E+00,
     0.9720000000000000E+00,
     0.4361908850559777E+00,
     0.1516409096347099E+00,
     0.8978271484375000E-01,
     0.1000000000000000E+01,
     0.5000000000000000E+00,
     0.4598773297575791E+00,
     0.2146816102371739E+00,
     0.9507364826957875E+00,
     0.5000000000000000E+00,
     0.8979413687105918E+00,
     0.2241297491808366E+00,
     0.7586405487192086E+00,
     0.7001783247477069E+00,
     0.5131670194948620E-01,
     0.1055728090000841E+00,
     0.1633399734659245E+00,
     0.2254033307585166E+00,
     0.3600000000000000E+00,
     0.4880000000000000E+00,
     0.5904000000000000E+00,
     0.6723200000000000E+00,
     0.2160000000000000E+00,
     0.8370000000000000E-01,
     0.3078000000000000E-01,
     0.1093500000000000E-01 };

  static double x_vec[N_MAX] = {
     0.01E+00,
     0.10E+00,
     1.00E+00,
     0.00E+00,
     0.01E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     0.50E+00,
     0.90E+00,
     0.50E+00,
     1.00E+00,
     0.50E+00,
     0.80E+00,
     0.60E+00,
     0.80E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.70E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void beta_inc_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_INC_VALUES returns some values of the incomplete Beta function.
//
//  Discussion:
//
//    The incomplete Beta function may be written
//
//      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
//                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
//
//    Thus,
//
//      BETA_INC(A,B,0.0) = 0.0;
//      BETA_INC(A,B,1.0) = 1.0
//
//    The incomplete Beta function is also sometimes called the
//    "modified" Beta function, or the "normalized" Beta function
//    or the Beta CDF (cumulative density function.
//
//    In Mathematica, the function can be evaluated by:
//
//      BETA[X,A,B] / BETA[A,B]
//
//    The function can also be evaluated by using the Statistics package:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = BetaDistribution [ a, b ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Karl Pearson,
//    Tables of the Incomplete Beta Function,
//    Cambridge University Press, 1968.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 42

  static double a_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.5E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     20.0E+00,
     30.0E+00,
     30.0E+00,
     40.0E+00,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01 };

  static double b_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      2.0E+00,
      5.0E+00,
      0.5E+00,
      5.0E+00,
      5.0E+00,
     10.0E+00,
      5.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
     20.0E+00,
     10.0E+00,
     10.0E+00,
     20.0E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01 };

  static double fx_vec[N_MAX] = {
     0.6376856085851985E-01,
     0.2048327646991335E+00,
     0.1000000000000000E+01,
     0.0000000000000000E+00,
     0.5012562893380045E-02,
     0.5131670194948620E-01,
     0.2928932188134525E+00,
     0.5000000000000000E+00,
     0.2800000000000000E-01,
     0.1040000000000000E+00,
     0.2160000000000000E+00,
     0.3520000000000000E+00,
     0.5000000000000000E+00,
     0.6480000000000000E+00,
     0.7840000000000000E+00,
     0.8960000000000000E+00,
     0.9720000000000000E+00,
     0.4361908850559777E+00,
     0.1516409096347099E+00,
     0.8978271484375000E-01,
     0.1000000000000000E+01,
     0.5000000000000000E+00,
     0.4598773297575791E+00,
     0.2146816102371739E+00,
     0.9507364826957875E+00,
     0.5000000000000000E+00,
     0.8979413687105918E+00,
     0.2241297491808366E+00,
     0.7586405487192086E+00,
     0.7001783247477069E+00,
     0.5131670194948620E-01,
     0.1055728090000841E+00,
     0.1633399734659245E+00,
     0.2254033307585166E+00,
     0.3600000000000000E+00,
     0.4880000000000000E+00,
     0.5904000000000000E+00,
     0.6723200000000000E+00,
     0.2160000000000000E+00,
     0.8370000000000000E-01,
     0.3078000000000000E-01,
     0.1093500000000000E-01 };

  static double x_vec[N_MAX] = {
     0.01E+00,
     0.10E+00,
     1.00E+00,
     0.00E+00,
     0.01E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     0.50E+00,
     0.90E+00,
     0.50E+00,
     1.00E+00,
     0.50E+00,
     0.80E+00,
     0.60E+00,
     0.80E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.70E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void beta_log_values ( int &n_data, double &x, double &y, double &fxy )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_LOG_VALUES returns some values of the logarithm of the Beta function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[Beta[x]]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, &Y, the arguments of the function.
//
//    Output, double &FXY, the value of the function.
//
{
# define N_MAX 17

  static double fxy_vec[N_MAX] = {
      0.1609437912434100E+01,
      0.9162907318741551E+00,
      0.5108256237659907E+00,
      0.2231435513142098E+00,
      0.1609437912434100E+01,
      0.9162907318741551E+00,
      0.0000000000000000E+00,
     -0.1791759469228055E+01,
     -0.3401197381662155E+01,
     -0.4941642422609304E+01,
     -0.6445719819385578E+01,
     -0.3737669618283368E+01,
     -0.5123963979403259E+01,
     -0.6222576268071369E+01,
     -0.7138866999945524E+01,
     -0.7927324360309794E+01,
     -0.9393661429103221E+01 };

  static double x_vec[N_MAX] = {
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     7.0E+00 };

  static double y_vec[N_MAX] = {
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     0.2E+00,
     0.4E+00,
     1.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     6.0E+00,
     7.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    y = 0.0;
    fxy = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    y = y_vec[n_data-1];
    fxy = fxy_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void beta_noncentral_cdf_values ( int &n_data, double &a, double &b,
  double &lambda, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
//
//  Discussion:
//
//    The values presented here are taken from the reference, where they
//    were given to a limited number of decimal places.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    R Chattamvelli, R Shanmugam,
//    Algorithm AS 310:
//    Computing the Non-central Beta Distribution Function,
//    Applied Statistics,
//    Volume 46, Number 1, 1997, pages 146-156.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the shape parameters.
//
//    Output, double &LAMBDA, the noncentrality parameter.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 25

  static double a_vec[N_MAX] = {
        5.0,
        5.0,
        5.0,
       10.0,
       10.0,
       10.0,
       20.0,
       20.0,
       20.0,
       10.0,
       10.0,
       15.0,
       20.0,
       20.0,
       20.0,
       30.0,
       30.0,
       10.0,
       10.0,
       10.0,
       15.0,
       10.0,
       12.0,
       30.0,
       35.0 };
  static double b_vec[N_MAX] = {
        5.0,
        5.0,
        5.0,
       10.0,
       10.0,
       10.0,
       20.0,
       20.0,
       20.0,
       20.0,
       10.0,
        5.0,
       10.0,
       30.0,
       50.0,
       20.0,
       40.0,
        5.0,
       10.0,
       30.0,
       20.0,
        5.0,
       17.0,
       30.0,
       30.0 };
  static double fx_vec[N_MAX] = {
       0.4563021,
       0.1041337,
       0.6022353,
       0.9187770,
       0.6008106,
       0.0902850,
       0.9998655,
       0.9925997,
       0.9641112,
       0.9376626573,
       0.7306817858,
       0.1604256918,
       0.1867485313,
       0.6559386874,
       0.9796881486,
       0.1162386423,
       0.9930430054,
       0.0506899273,
       0.1030959706,
       0.9978417832,
       0.2555552369,
       0.0668307064,
       0.0113601067,
       0.7813366615,
       0.8867126477 };
  static double lambda_vec[N_MAX] = {
        54.0,
       140.0,
       170.0,
        54.0,
       140.0,
       250.0,
        54.0,
       140.0,
       250.0,
       150.0,
       120.0,
        80.0,
       110.0,
        65.0,
       130.0,
        80.0,
       130.0,
        20.0,
        54.0,
        80.0,
       120.0,
        55.0,
        64.0,
       140.0,
        20.0 };
  static double x_vec[N_MAX] = {
       0.8640,
       0.9000,
       0.9560,
       0.8686,
       0.9000,
       0.9000,
       0.8787,
       0.9000,
       0.9220,
       0.868,
       0.900,
       0.880,
       0.850,
       0.660,
       0.720,
       0.720,
       0.800,
       0.644,
       0.700,
       0.780,
       0.760,
       0.795,
       0.560,
       0.800,
       0.670 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    lambda = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    lambda = lambda_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void beta_values ( int &n_data, double &x, double &y, double &fxy )

//****************************************************************************80
//
//  Purpose:
//
//    BETA_VALUES returns some values of the Beta function.
//
//  Discussion:
//
//    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
//
//    Both X and Y must be greater than 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      Beta[X,Y]
//
//  Properties:
//
//    Beta(X,Y) = Beta(Y,X).
//    Beta(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
//    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, &Y, the arguments of the function.
//
//    Output, double &FXY, the value of the function.
//
{
# define N_MAX 17

  static double b_vec[N_MAX] = {
     0.5000000000000000E+01,
     0.2500000000000000E+01,
     0.1666666666666667E+01,
     0.1250000000000000E+01,
     0.5000000000000000E+01,
     0.2500000000000000E+01,
     0.1000000000000000E+01,
     0.1666666666666667E+00,
     0.3333333333333333E-01,
     0.7142857142857143E-02,
     0.1587301587301587E-02,
     0.2380952380952381E-01,
     0.5952380952380952E-02,
     0.1984126984126984E-02,
     0.7936507936507937E-03,
     0.3607503607503608E-03,
     0.8325008325008325E-04 };

  static double x_vec[N_MAX] = {
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     6.0E+00,
     7.0E+00 };

  static double y_vec[N_MAX] = {
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     0.2E+00,
     0.4E+00,
     1.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     5.0E+00,
     6.0E+00,
     7.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    y = 0.0;
    fxy = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    y = y_vec[n_data-1];
    fxy = b_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void binomial_values ( int &n_data, int &a, int &b, int &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_VALUES returns some values of the binomial coefficients.
//
//  Discussion:
//
//    The formula for the binomial coefficient is
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//    In Mathematica, the function can be evaluated by:
//
//      Binomial[n,k]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &A, &B, the arguments of the function.
//
//    Output, int &FX, the value of the function.
//
{
# define N_MAX 20

  static int a_vec[N_MAX] = {
     1,  6,  6,  6, 15,
    15, 15, 15, 15, 15,
    15, 25, 25, 25, 25,
    25, 25, 25, 25, 25  };

  static int b_vec[N_MAX] = {
     0,  1,  3,  5,  1,
     3,  5,  7,  9, 11,
    13,  1,  3,  5,  7,
     9, 11, 13, 15, 17 };

  static int fx_vec[N_MAX] = {
           1,
           6,
          20,
           6,
          15,
         455,
        3003,
        6435,
        5005,
        1365,
         105,
          25,
        2300,
       53130,
      480700,
     2042975,
     4457400,
     5200300,
     3268760,
     1081575 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0;
    b = 0;
    fx = 0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void binomial_cdf_values ( int &n_data, int &a, double &b, int &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    BINOMIAL_CDF_VALUES returns some values of the binomial CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = BinomialDistribution [ n, p ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &A, a parameter of the function.
//
//    Output, double &B, a parameter of the function.
//
//    Output, int &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 17

  static int a_vec[N_MAX] = {
     2,  2,  2,  2,
     2,  4,  4,  4,
     4, 10, 10, 10,
    10, 10, 10, 10,
    10 };

  static double b_vec[N_MAX] = {
     0.05E+00,
     0.05E+00,
     0.05E+00,
     0.50E+00,
     0.50E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.05E+00,
     0.10E+00,
     0.15E+00,
     0.20E+00,
     0.25E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00  };

  static double fx_vec[N_MAX] = {
     0.9025000000000000E+00,
     0.9975000000000000E+00,
     0.1000000000000000E+01,
     0.2500000000000000E+00,
     0.7500000000000000E+00,
     0.3164062500000000E+00,
     0.7382812500000000E+00,
     0.9492187500000000E+00,
     0.9960937500000000E+00,
     0.9999363101685547E+00,
     0.9983650626000000E+00,
     0.9901259090013672E+00,
     0.9672065024000000E+00,
     0.9218730926513672E+00,
     0.8497316674000000E+00,
     0.6331032576000000E+00,
     0.3769531250000000E+00 };

  static int x_vec[N_MAX] = {
     0, 1, 2, 0,
     1, 0, 1, 2,
     3, 4, 4, 4,
     4, 4, 4, 4,
     4 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0;
    b = 0.0;
    x = 0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void bivariate_normal_cdf_values ( int &n_data, double &x, double &y,
  double &r, double &fxy )

//****************************************************************************80
//
//  Purpose:
//
//    BIVARIATE_NORMAL_CDF_VALUES returns some values of the bivariate normal CDF.
//
//  Discussion:
//
//    FXY is the probability that two variables A and B, which are
//    related by a bivariate normal distribution with correlation R,
//    respectively satisfy A <= X and B <= Y.
//
//    Mathematica can evaluate the bivariate normal CDF via the commands:
//
//      <<MultivariateStatistics`
//      cdf = CDF[MultinormalDistribution[{0,0}{{1,r},{r,1}}],{x,y}]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 November 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    National Bureau of Standards,
//    Tables of the Bivariate Normal Distribution and Related Functions,
//    NBS, Applied Mathematics Series, Number 50, 1959.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, &Y, the parameters of the function.
//
//    Output, double &R, the correlation value.
//
//    Output, double &FXY, the value of the function.
//
{
#  define N_MAX 41

  static double fxy_vec[N_MAX] = {
  0.02260327218569867E+00,
  0.1548729518584100E+00,
  0.4687428083352184E+00,
  0.7452035868929476E+00,
  0.8318608306874188E+00,
  0.8410314261134202E+00,
  0.1377019384919464E+00,
  0.1621749501739030E+00,
  0.1827411243233119E+00,
  0.2010067421506235E+00,
  0.2177751155265290E+00,
  0.2335088436446962E+00,
  0.2485057781834286E+00,
  0.2629747825154868E+00,
  0.2770729823404738E+00,
  0.2909261168683812E+00,
  0.3046406378726738E+00,
  0.3183113449213638E+00,
  0.3320262544108028E+00,
  0.3458686754647614E+00,
  0.3599150462310668E+00,
  0.3742210899871168E+00,
  0.3887706405282320E+00,
  0.4032765198361344E+00,
  0.4162100291953678E+00,
  0.6508271498838664E+00,
  0.8318608306874188E+00,
  0.0000000000000000,
  0.1666666666539970,
  0.2500000000000000,
  0.3333333333328906,
  0.5000000000000000,
  0.7452035868929476,
  0.1548729518584100,
  0.1548729518584100,
  0.06251409470431653,
  0.7452035868929476,
  0.1548729518584100,
  0.1548729518584100,
  0.06251409470431653,
  0.6337020457912916 };
  static double r_vec[N_MAX] = {
     0.500,  0.500,  0.500,  0.500,  0.500,
     0.500, -0.900, -0.800, -0.700, -0.600,
    -0.500, -0.400, -0.300, -0.200, -0.100,
     0.000,  0.100,  0.200,  0.300,  0.400,
     0.500,  0.600,  0.700,  0.800,  0.900,
     0.673,  0.500, -1.000, -0.500,  0.000,
     0.500,  1.000,  0.500,  0.500,  0.500,
     0.500,  0.500,  0.500,  0.500,  0.500,
     0.500 };
  static double x_vec[N_MAX] = {
    -2.0, -1.0,  0.0,  1.0,  2.0,
     3.0, -0.2, -0.2, -0.2, -0.2,
    -0.2, -0.2, -0.2, -0.2, -0.2,
    -0.2, -0.2, -0.2, -0.2, -0.2,
    -0.2, -0.2, -0.2, -0.2, -0.2,
     1.0,  2.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  1.0,  1.0, -1.0,
    -1.0,  1.0,  1.0, -1.0, -1.0,
     0.7071067811865475 };
  static double y_vec[N_MAX] = {
     1.0,  1.0,  1.0,  1.0,  1.0,
     1.0,  0.5,  0.5,  0.5,  0.5,
     0.5,  0.5,  0.5,  0.5,  0.5,
     0.5,  0.5,  0.5,  0.5,  0.5,
     0.5,  0.5,  0.5,  0.5,  0.5,
     0.5,  1.0,  0.0,  0.0,  0.0,
     0.0,  0.0,  1.0, -1.0,  1.0,
    -1.0,  1.0, -1.0,  1.0, -1.0,
     0.7071067811865475 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    r = 0.0;
    x = 0.0;
    y = 0.0;
    fxy = 0.0;
  }
  else
  {
    r = r_vec[n_data-1];
    x = x_vec[n_data-1];
    y = y_vec[n_data-1];
    fxy = fxy_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void catalan_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    CATALAN_VALUES returns some values of the Catalan numbers.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Binomial[2*n,n] / ( n + 1 )
//
//  First values:
//
//     C(0)     1
//     C(1)     1
//     C(2)     2
//     C(3)     5
//     C(4)    14
//     C(5)    42
//     C(6)   132
//     C(7)   429
//     C(8)  1430
//     C(9)  4862
//    C(10) 16796
//
//  Formula:
//
//    C(N) = (2*N)! / ( (N+1) * (N!) * (N!) )
//         = 1 / (N+1) * COMB ( 2N, N )
//         = 1 / (2N+1) * COMB ( 2N+1, N+1).
//
//  Recursion:
//
//    C(N) = 2 * (2*N-1) * C(N-1) / (N+1)
//    C(N) = sum ( 1 <= I <= N-1 ) C(I) * C(N-I)
//
//  Discussion:
//
//    The Catalan number C(N) counts:
//
//    1) the number of binary trees on N vertices;
//    2) the number of ordered trees on N+1 vertices;
//    3) the number of full binary trees on 2N+1 vertices;
//    4) the number of well formed sequences of 2N parentheses;
//    5) the number of ways 2N ballots can be counted, in order,
//       with N positive and N negative, so that the running sum
//       is never negative;
//    6) the number of standard tableaus in a 2 by N rectangular Ferrers diagram;
//    7) the number of monotone functions from [1..N} to [1..N} which
//       satisfy f(i) <= i for all i;
//    8) the number of ways to triangulate a polygon with N+2 vertices.
//
//  Example:
//
//    N = 3
//
//    ()()()
//    ()(())
//    (()())
//    (())()
//    ((()))
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Catalan number.
//
//    Output, int &C, the value of the Catalan number.
//
{
# define N_MAX 11

  static int c_vec[N_MAX] = {
    1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,  4, 5,  6,  7,  8,  9,  10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cauchy_cdf_values ( int &n_data, double &mu, double &sigma, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CAUCHY_CDF_VALUES returns some values of the Cauchy CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = CauchyDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &SIGMA, the variance of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.8524163823495667E+00,
     0.9220208696226307E+00,
     0.9474315432887466E+00,
     0.6475836176504333E+00,
     0.6024163823495667E+00,
     0.5779791303773693E+00,
     0.5628329581890012E+00,
     0.6475836176504333E+00,
     0.5000000000000000E+00,
     0.3524163823495667E+00,
     0.2500000000000000E+00 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cbrt_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CBRT_VALUES returns some values of the cube root function.
//
//  Discussion:
//
//    CBRT(X) = real number Y such that Y * Y * Y = X.
//
//    In Mathematica, the function can be evaluated by:
//
//      Sign[x] * ( Abs[x] )^(1/3)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output double FX, the value of the function.
//
{
# define N_MAX 14

  static double fx_vec[N_MAX] = {
    0.0000000000000000E+00,
   -0.0020082988563383484484E+00,
    0.44814047465571647087E+00,
   -0.46415888336127788924E+00,
    0.73680629972807732116E+00,
   -1.0000000000000000000E+00,
    1.2599210498948731648E+00,
   -1.4422495703074083823E+00,
    1.4645918875615232630E+00,
   -2.6684016487219448673E+00,
    3.0723168256858472933E+00,
   -4.1408177494228532500E+00,
    4.5947008922070398061E+00,
   -497.93385921817447440E+00 };

  static double x_vec[N_MAX] = {
     0.0000000000000000E+00,
    -0.8100000073710001E-08,
     0.9000000000000000E-01,
    -0.1000000000000000E+00,
     0.4000000000000000E+00,
    -0.1000000000000000E+01,
     0.2000000000000000E+01,
    -0.3000000000000000E+01,
     0.3141592653589793E+01,
    -0.1900000000000000E+02,
     0.2900000000000000E+02,
    -0.7100000000000000E+02,
     0.9700000000000000E+02,
    -0.1234567890000000E+09 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_t_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_T_POLY_VALUES returns values of Chebyshev polynomials T(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevT[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.8000000000000000E+00,
      0.2800000000000000E+00,
     -0.3520000000000000E+00,
     -0.8432000000000000E+00,
     -0.9971200000000000E+00,
     -0.7521920000000000E+00,
     -0.2063872000000000E+00,
      0.4219724800000000E+00,
      0.8815431680000000E+00,
      0.9884965888000000E+00,
      0.7000513740800000E+00,
      0.1315856097280000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_u_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_U_POLY_VALUES returns values of Chebyshev polynomials U(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ChebyshevU[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.1600000000000000E+01,
      0.1560000000000000E+01,
      0.8960000000000000E+00,
     -0.1264000000000000E+00,
     -0.1098240000000000E+01,
     -0.1630784000000000E+01,
     -0.1511014400000000E+01,
     -0.7868390400000000E+00,
      0.2520719360000000E+00,
      0.1190154137600000E+01,
      0.1652174684160000E+01,
      0.1453325357056000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_v_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_V_POLY_VALUES returns values of Chebyshev polynomials V(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      u = Sqrt[(x+1)/2],
//      ChebyshevT[2*n+1,u] / u
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     1.0000000000000000E+00, 
     0.6000000000000000E+00, 
    -0.0400000000000000E+00, 
    -0.6640000000000000E+00, 
    -1.0224000000000000E+00, 
    -0.9718400000000000E+00, 
    -0.5325440000000000E+00, 
     0.1197696000000000E+00, 
     0.7241753600000000E+00, 
     1.0389109760000000E+00, 
     0.9380822016000000E+00, 
     0.4620205465600000E+00, 
    -0.1988493271040000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cheby_w_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHEBY_W_POLY_VALUES returns values of Chebyshev polynomials W(n,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      u = Sqrt[(x+1)/2],
//      ChebyshevU[2*n,u]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     1.000000000000000E+00, 
     2.600000000000000E+00, 
     3.160000000000000E+00, 
     2.456000000000000E+00, 
     0.769600000000000E+00, 
    -1.224640000000000E+00, 
    -2.729024000000000E+00, 
    -3.141798400000000E+00, 
    -2.297853440000000E+00, 
    -0.534767104000000E+00, 
     1.442226073600000E+00, 
     2.842328821760000E+00, 
     3.105500041216000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12 };

  static double x_vec[N_MAX] = {
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00,
     0.8E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void chi_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_VALUES returns some values of the hyperbolic cosine integral function.
//
//  Discussion:
//
//    The hyperbolic cosine integral is defined by
//
//      CHI(X) = gamma + log ( x )
//        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
//
//    where gamma is Euler's constant.
//
//    In Mathematica, the function can be evaluated by:
//
//      CoshIntegral[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
    -0.05277684495649362,
     0.1577508933739787,
     0.3455691756953907,
     0.5183999848333915,
     0.6813138871854339,
     0.8378669409802082,
     1.141841924170595,
     1.445494075789644,
     1.759505807660965,
     2.092577214062032,
     2.452666922646915,
     3.524425488354165,
     4.960392094765610,
     6.959191927647393,
     9.813547558823186,
    13.96581164859243 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void chi_square_cdf_values ( int &n_data, int &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ChiSquareDistribution [ df ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
//
# define N_MAX 21

  static int a_vec[N_MAX] = {
     1,  2,  1,  2,
     1,  2,  3,  4,
     1,  2,  3,  4,
     5,  3,  3,  3,
     3,  3, 10, 10,
    10 };

  static double fx_vec[N_MAX] = {
     0.7965567455405796E-01,
     0.4987520807317687E-02,
     0.1124629160182849E+00,
     0.9950166250831946E-02,
     0.4729107431344619E+00,
     0.1812692469220181E+00,
     0.5975750516063926E-01,
     0.1752309630642177E-01,
     0.6826894921370859E+00,
     0.3934693402873666E+00,
     0.1987480430987992E+00,
     0.9020401043104986E-01,
     0.3743422675270363E-01,
     0.4275932955291202E+00,
     0.6083748237289110E+00,
     0.7385358700508894E+00,
     0.8282028557032669E+00,
     0.8883897749052874E+00,
     0.1721156299558408E-03,
     0.3659846827343712E-02,
     0.1857593622214067E-01 };

  static double x_vec[N_MAX] = {
     0.01E+00,
     0.01E+00,
     0.02E+00,
     0.02E+00,
     0.40E+00,
     0.40E+00,
     0.40E+00,
     0.40E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     3.00E+00,
     4.00E+00,
     5.00E+00,
     6.00E+00,
     1.00E+00,
     2.00E+00,
     3.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void chi_square_noncentral_cdf_values ( int &n_data, int &df, double &lambda,
  double &x, double &cdf )

//****************************************************************************80
//
//  Purpose:
//
//    CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralChiSquareDistribution [ df, lambda ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &DF, the number of degrees of freedom.
//
//    Output, double &LAMBDA, the noncentrality parameter.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &CDF, the noncentral chi CDF.
//
{
# define N_MAX 28

  static double cdf_vec[N_MAX] = {
     0.8399444269398261E+00,
     0.6959060300435139E+00,
     0.5350879697078847E+00,
     0.7647841496310313E+00,
     0.6206436532195436E+00,
     0.4691667375373180E+00,
     0.3070884345937569E+00,
     0.2203818092990903E+00,
     0.1500251895581519E+00,
     0.3071163194335791E-02,
     0.1763982670131894E-02,
     0.9816792594625022E-03,
     0.1651753140866208E-01,
     0.2023419573950451E-03,
     0.4984476352854074E-06,
     0.1513252400654827E-01,
     0.2090414910614367E-02,
     0.2465021206048452E-03,
     0.2636835050342939E-01,
     0.1857983220079215E-01,
     0.1305736595486640E-01,
     0.5838039534819351E-01,
     0.4249784402463712E-01,
     0.3082137716021596E-01,
     0.1057878223400849E+00,
     0.7940842984598509E-01,
     0.5932010895599639E-01,
     0.2110395656918684E+00 };

  static int df_vec[N_MAX] = {
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
      1,   2,   3,
     60,  80, 100,
      1,   2,   3,
     10,  10,  10,
     10,  10,  10,
     10,  10,  10,
      8 };

  static double lambda_vec[N_MAX] = {
     0.5E+00,
     0.5E+00,
     0.5E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
    20.0E+00,
    20.0E+00,
    20.0E+00,
    30.0E+00,
    30.0E+00,
    30.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     0.5E+00 };

  static double x_vec[N_MAX] = {
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
      3.000E+00,
     60.000E+00,
     60.000E+00,
     60.000E+00,
      0.050E+00,
      0.050E+00,
      0.050E+00,
      4.000E+00,
      4.000E+00,
      4.000E+00,
      5.000E+00,
      5.000E+00,
      5.000E+00,
      6.000E+00,
      6.000E+00,
      6.000E+00,
      5.000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    lambda = 0.0;
    df = 0;
    cdf = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    lambda = lambda_vec[n_data-1];
    df = df_vec[n_data-1];
    cdf = cdf_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ci_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CI_VALUES returns some values of the cosine integral function.
//
//  Discussion:
//
//    The cosine integral is defined by
//
//      CI(X) = - integral ( X <= T < +oo ) ( cos ( T ) ) / T  dT
//
//    In Mathematica, the function can be evaluated by:
//
//      CosIntegral[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     -0.1777840788066129E+00,
     -0.2227070695927976E-01,
      0.1005147070088978E+00,
      0.1982786159524672E+00,
      0.2760678304677729E+00,
      0.3374039229009681E+00,
      0.4204591828942405E+00,
      0.4620065850946773E+00,
      0.4717325169318778E+00,
      0.4568111294183369E+00,
      0.4229808287748650E+00,
      0.2858711963653835E+00,
      0.1196297860080003E+00,
     -0.3212854851248112E-01,
     -0.1409816978869304E+00,
     -0.1934911221017388E+00 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cin_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CIN_VALUES returns some values of the alternate cosine integral function.
//
//  Discussion:
//
//    The alternate cosine integral is defined by
//
//      CIN(X) = gamma + log(X) + integral ( 0 <= T <= X ) ( cos ( T ) - 1 ) / T  dT
//
//    In Mathematica, the function can be evaluated by:
//
//      EulerGamma + Log[x] - CosIntegral[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.6185256314820045E-01,
     0.8866074809482194E-01,
     0.1200260139539026E+00,
     0.1557934976348559E+00,
     0.1957873187759337E+00,
     0.2398117420005647E+00,
     0.3390780388012470E+00,
     0.4516813164280685E+00,
     0.5754867772153906E+00,
     0.7081912003853150E+00,
     0.8473820166866132E+00,
     0.1207635200410304E+01,
     0.1556198167561642E+01,
     0.1862107181909382E+01,
     0.2104491723908354E+01,
     0.2274784183779546E+01 };

  static double x_vec[N_MAX] = {
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.5E+00,
     3.0E+00,
     3.5E+00,
     4.0E+00,
     4.5E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cinh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CINH_VALUES returns some values of the alternate hyperbolic cosine integral.
//
//  Discussion:
//
//    The alternate hyperbolic cosine integral is defined by
//
//      CINH(X) = integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
//
//    In Mathematica, the function can be evaluated by:
//
//      Integrate [ ( Cosh[t] - 1 ) / t, { t, 0, x } ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 17

  static double fx_vec[N_MAX] = {
     0.00000000000000000,
     0.06315467070191883,
     0.09136085223843649,
     0.1250284547325902,
     0.1643278712460683,
     0.2094587379417273,
     0.2606512760786754,
     0.3823047024751071,
     0.5318061742668980,
     0.7122865135136963,
     0.9275748842583805,
     1.182304077185436,
     2.030919091578478,
     3.284564141195967,
     5.129213294250493,
     7.850037532801762,
    11.88451858691463 };
  static double x_vec[N_MAX] = {
     0.0,
     0.5,
     0.6,
     0.7,
     0.8,
     0.9,
     1.0,
     1.2,
     1.4,
     1.6,
     1.8,
     2.0,
     2.5,
     3.0,
     3.5,
     4.0,
     4.5 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void clausen_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CLAUSEN_VALUES returns some values of the Clausen's integral.
//
//  Discussion:
//
//    The function is defined by:
//
//      CLAUSEN(x) = integral ( 0 <= t <= x ) -ln ( 2 * sin ( t / 2 ) ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.14137352886760576684E-01,
      0.13955467081981281934E+00,
     -0.38495732156574238507E+00,
      0.84831187770367927099E+00,
      0.10139591323607685043E+01,
     -0.93921859275409211003E+00,
      0.72714605086327924743E+00,
      0.43359820323553277936E+00,
     -0.98026209391301421161E-01,
     -0.56814394442986978080E+00,
     -0.70969701784448921625E+00,
      0.99282013254695671871E+00,
     -0.98127747477447367875E+00,
     -0.64078266570172320959E+00,
      0.86027963733231192456E+00,
      0.39071647608680211043E+00,
      0.47574793926539191502E+00,
      0.10105014481412878253E+01,
      0.96332089044363075154E+00,
     -0.61782699481929311757E+00 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
      -0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
      -1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
      -3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
      -5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
     -10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
     -30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void clebsch_gordan_values ( int &n_data, double &j1, double &j2, double &j3,
  double &m1, double &m2, double &m3, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    CLEBSCH_GORDAN_VALUES returns some values of the Clebsch-Gordan function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ClebschGordan[{j1,m1},{j2,m2},{j3,m3}]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &J1, &J2, &J3, &M1, &M2, &M3, the arguments
//    of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.7071067811865475,
     1.000000000000000,
     0.5773502691896258,
    -0.2581988897471611,
    -0.6324555320336759,
    -0.7745966692414834,
     0.4082482904638630,
     0.8164965809277260,
     0.5345224838248488,
     0.2672612419124244,
     0.8944271909999159,
     0.3380617018914066 };
  static double j1_vec[N_MAX] = {
    0.5,
    0.5,
    0.5,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    2.0,
    2.0,
    1.5,
    1.5 };
  static double j2_vec[N_MAX] = {
    0.5,
    0.5,
    1.0,
    1.5,
    1.5,
    1.5,
    1.0,
    1.0,
    2.0,
    2.0,
    2.0,
    2.0 };
  static double j3_vec[N_MAX] = {
    1.0,
    1.0,
    1.5,
    1.5,
    1.5,
    1.5,
    2.0,
    2.0,
    2.0,
    2.0,
    2.5,
    3.5 };
  static double m1_vec[N_MAX] = {
    0.5,
    0.5,
   -0.5,
    0.0,
   -1.0,
    0.0,
    1.0,
    0.0,
    2.0,
    1.0,
    0.5,
    1.5 };
  static double m2_vec[N_MAX] = {
   -0.5,
    0.5,
    1.0,
    0.5,
    1.5,
    1.5,
   -1.0,
    0.0,
   -2.0,
   -1.0,
    1.0,
   -1.0 };
  static double m3_vec[N_MAX] = {
   -0.5,
    0.5,
    1.0,
    0.5,
    1.5,
    1.5,
   -1.0,
    0.0,
   -2.0,
   -1.0,
    1.0,
   -1.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    j1 = 0.0;
    j2 = 0.0;
    j3 = 0.0;
    m1 = 0.0;
    m2 = 0.0;
    m3 = 0.0;
    fx = 0.0;
  }
  else
  {
    j1 = j1_vec[n_data-1];
    j2 = j2_vec[n_data-1];
    j3 = j3_vec[n_data-1];
    m1 = m1_vec[n_data-1];
    m2 = m2_vec[n_data-1];
    m3 = m3_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void collatz_count_values ( int &n_data, int &n, int &count )

//****************************************************************************80
//
//  Purpose:
//
//    COLLATZ_COUNT_VALUES returns some values of the Collatz count function.
//
//  Discussion:
//
//    The rules for generation of the Collatz sequence are recursive.
//    If T is the current entry of the sequence, (T is
//    assumed to be a positive integer), then the next
//    entry, U is determined as follows:
//
//      if T is 1 (or less)
//        terminate the sequence;
//      else if T is even
//        U = T/2.
//      else (if T is odd and not 1)
//        U = 3*T+1;
//
//    The Collatz count is the length of the Collatz sequence for a given
//    starting value.  By convention, we include the initial value in the
//    count, so the minimum value of the count is 1.
//
//     N  Sequence                                                 Count
//
//     1                                                               1
//     2   1                                                           2
//     3  10,  5, 16,  8,  4,  2,  1                                   8
//     4   2   1                                                       3
//     5  16,  8,  4,  2,  1                                           6
//     6   3, 10,  5, 16,  8,  4,  2,  1                               9
//     7  22, 11, 34, 17, 52, 26, 13, 40, 20, 10, 5, 16, 8, 4, 2, 1   17
//     8   4,  2,  1                                                   4
//     9  28, 14,  7, ...                                             20
//    10   5, 16,  8,  4,  2,  1                                       7
//    11  34, 17, 52, 26, 13, 40, 20, 10,  5, 16, 8, 4, 2, 1          15
//    12   6,  3, 10,  5, 16,  8,  4,  2,  1                          10
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 March 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Eric Weisstein,
//    "The Collatz Problem",
//    CRC Concise Encyclopedia of Mathematics,
//    CRC 1998.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the initial value of a Collatz sequence.
//
//    Output, int &COUNT, the length of the Collatz sequence starting
//    with N.
//
{
# define N_MAX 20

  static int count_vec[N_MAX] = {
      1,   2,   8,   3,   6,   9,   17,   4,  20,   7,
    112,  25,  26,  27,  17,  28,  111,  18,  83,  29 };
  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     27,  50, 100, 200, 300, 400, 500, 600, 700, 800 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    count = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    count = count_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cos_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COS_VALUES returns some values of the cosine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Cos[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     1.0000000000000000000,
     0.96592582628906828675,
     0.87758256189037271612,
     0.86602540378443864676,
     0.70710678118654752440,
     0.54030230586813971740,
     0.50000000000000000000,
     0.00000000000000000000,
    -0.41614683654714238700,
    -0.98999249660044545727,
    -1.0000000000000000000,
    -0.65364362086361191464,
     0.28366218546322626447 };

  static double x_vec[N_MAX] = {
    0.0000000000000000000,
    0.26179938779914943654,
    0.50000000000000000000,
    0.52359877559829887308,
    0.78539816339744830962,
    1.0000000000000000000,
    1.0471975511965977462,
    1.5707963267948966192,
    2.0000000000000000000,
    3.0000000000000000000,
    3.1415926535897932385,
    4.0000000000000000000,
    5.0000000000000000000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cos_degree_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COS_DEGREE_VALUES: values of the cosine function for degree arguments.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Cos[x Degree]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
     0.99619469809174553230,
     1.0000000000000000000,
     0.99984769515639123916,
     0.99939082701909573001,
     0.99862953475457387378,
     0.99756405025982424761,
     0.99619469809174553230,
     0.98480775301220805937,
     0.96592582628906828675,
     0.86602540378443864676,
     0.70710678118654752440,
     0.50000000000000000000,
     0.25881904510252076235,
     0.087155742747658173558,
     0.069756473744125300776,
     0.052335956242943832722,
     0.034899496702500971646,
     0.017452406437283512819,
     0.000000000000000000000,
    -0.017452406437283512819,
    -0.25881904510252076235,
    -1.0000000000000000000 };
  static double x_vec[N_MAX] = {
     -5.0,
      0.0,
      1.0,
      2.0,
      3.0,
      4.0,
      5.0,
     10.0,
     15.0,
     30.0,
     45.0,
     60.0,
     75.0,
     85.0,
     86.0,
     87.0,
     88.0,
     89.0,
     90.0,
     91.0,
    105.0,
    180.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cos_power_int_values ( int &n_data, double &a, double &b, int &n,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COS_POWER_INT_VALUES returns some values of the sine power integral.
//
//  Discussion:
//
//    The function has the form
//
//      COS_POWER_INT(A,B,N) = integral ( A <= T <= B ) ( cos(T) )^N dt
//
//    In Mathematica, the function can be evaluated by:
//
//      Integrate [ ( Cos[x] )^n, { x, a, b } ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the limits of integration.
//
//    Output, int &N, the power.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double a_vec[N_MAX] = {
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00 };

  static double b_vec[N_MAX] = {
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793,
      3.141592653589793 };

  static double fx_vec[N_MAX] = {
     3.141592653589793, 
     0.0, 
     1.570796326794897, 
     0.0, 
     1.178097245096172, 
     0.0, 
     0.9817477042468104, 
     0.0, 
     0.8590292412159591, 
     0.0, 
     0.7731263170943632 };

  static int n_vec[N_MAX] = {
     0,
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cosh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COSH_VALUES returns some values of the hyperbolic cosine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Cosh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
      74.209948524787844444,
       1.5430806348152437785,
       1.0000000000000000000,
       1.0050041680558035990,
       1.0200667556190758463,
       1.0453385141288604850,
       1.0810723718384548093,
       1.1276259652063807852,
       1.1854652182422677038,
       1.2551690056309430182,
       1.3374349463048445980,
       1.4330863854487743878,
       1.5430806348152437785,
       3.7621956910836314596,
      10.067661995777765842,
      27.308232836016486629,
      74.209948524787844444,
   11013.232920103323140 };

  static double x_vec[N_MAX] = {
   -5.0,
   -1.0,
    0.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
   10.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cot_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    COT_VALUES returns some values of the cotangent function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Cot[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double fx_vec[N_MAX] = {
   11.972209353628661620,
    3.7320508075688772935,
    1.8304877217124519193,
    1.7320508075688772935,
    1.0000000000000000000,
    0.64209261593433070301,
    0.57735026918962576451,
    0.26794919243112270647,
    0.00000000000000000000,
    0.13165249758739585347,
    0.065543462815238228565,
   -0.45765755436028576375,
   -7.0152525514345334694,
    0.86369115445061661395,
   -0.29581291553274554043 };

  static double x_vec[N_MAX] = {
   0.083333333333333333333,
   0.26179938779914943654,
   0.50000000000000000000,
   0.52359877559829887308,
   0.78539816339744830962,
   1.0000000000000000000,
   1.0471975511965977462,
   1.3089969389957471827,
   1.5707963267948966192,
   1.4398966328953219010,
   1.5053464798451092601,
   2.0000000000000000000,
   3.0000000000000000000,
   4.0000000000000000000,
   5.0000000000000000000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void cp_values ( int &n_data, double &tc, double &p, double &cp )

//****************************************************************************80
//
//  Purpose:
//
//    CP_VALUES returns some values of the specific heat at constant pressure.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, pages 229-237.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &CP, the specific heat at constant pressure,
//    in KJ/(kg K).
//
{
# define N_MAX 24

  static double cp_vec[N_MAX] = {
     4.228E+00,
     2.042E+00,
     1.975E+00,
     2.013E+00,
     2.040E+00,
     2.070E+00,
     2.135E+00,
     2.203E+00,
     2.378E+00,
     2.541E+00,
     2.792E+00,
     2.931E+00,
     4.226E+00,
     4.223E+00,
     4.202E+00,
     4.177E+00,
     4.130E+00,
     4.089E+00,
     4.053E+00,
     4.021E+00,
     3.909E+00,
     3.844E+00,
     3.786E+00,
     2.890E+00  };

  static double p_vec[N_MAX] = {
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        5.0E+00,
       10.0E+00,
       50.0E+00,
      100.0E+00,
      200.0E+00,
      300.0E+00,
      400.0E+00,
      500.0E+00,
     1000.0E+00,
     1500.0E+00,
     2000.0E+00,
     5000.0E+00 };

  static double tc_vec[N_MAX] = {
        0.0E+00,
      100.0E+00,
      200.0E+00,
      300.0E+00,
      350.0E+00,
      400.0E+00,
      500.0E+00,
      600.0E+00,
      850.0E+00,
     1100.0E+00,
     1600.0E+00,
     2000.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    cp = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    cp = cp_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void dawson_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DAWSON_VALUES returns some values of Dawson's integral.
//
//  Discussion:
//
//    The definition of Dawson's integral is
//
//      D(X) = exp ( -X * X ) * integral ( 0 <= Y <= X ) exp ( Y * Y ) dY
//
//    Dawson's integral has a maximum at roughly
//
//      X = 0.9241388730
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[Pi] * Exp[-x^2] * I * Erf[I*x] / 2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 1998.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.9933599239785286E-01,
     0.1947510333680280E+00,
     0.2826316650213119E+00,
     0.3599434819348881E+00,
     0.4244363835020223E+00,
     0.4747632036629779E+00,
     0.5105040575592318E+00,
     0.5321017070563654E+00,
     0.5407243187262987E+00,
     0.5380795069127684E+00,
     0.5262066799705525E+00,
     0.5072734964077396E+00,
     0.4833975173848241E+00,
     0.4565072375268973E+00,
     0.4282490710853986E+00,
     0.3999398943230814E+00,
     0.3725593489740788E+00,
     0.3467727691148722E+00,
     0.3229743193228178E+00,
     0.3013403889237920E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void debye1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE1_VALUES returns some values of Debye's function of order 1.
//
//  Discussion:
//
//    The function is defined by:
//
//      DEBYE1(x) = 1 / x * integral ( 0 <= t <= x ) t / ( exp ( t ) - 1 ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.99951182471380889183E+00,
      0.99221462647120597836E+00,
      0.96918395997895308324E+00,
      0.88192715679060552968E+00,
      0.77750463411224827642E+00,
      0.68614531078940204342E+00,
      0.60694728460981007205E+00,
      0.53878956907785587703E+00,
      0.48043521957304283829E+00,
      0.38814802129793784501E+00,
      0.36930802829242526815E+00,
      0.32087619770014612104E+00,
      0.29423996623154246701E+00,
      0.27126046678502189985E+00,
      0.20523930310221503723E+00,
      0.16444346567994602563E+00,
      0.10966194482735821276E+00,
      0.82246701178200016086E-01,
      0.54831135561510852445E-01,
      0.32898681336964528729E-01 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void debye2_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE2_VALUES returns some values of Debye's function of order 2.
//
//  Discussion:
//
//    The function is defined by:
//
//      DEBYE2(x) = 2 / x^2 * integral ( 0 <= t <= x ) t^2 / ( exp ( t ) - 1 ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.99934911727904599738E+00,
      0.98962402299599181205E+00,
      0.95898426200345986743E+00,
      0.84372119334725358934E+00,
      0.70787847562782928288E+00,
      0.59149637225671282917E+00,
      0.49308264399053185014E+00,
      0.41079413579749669069E+00,
      0.34261396060786351671E+00,
      0.24055368752127897660E+00,
      0.22082770061202308232E+00,
      0.17232915939014138975E+00,
      0.14724346738730182894E+00,
      0.12666919046715789982E+00,
      0.74268805954862854626E-01,
      0.47971498020121871622E-01,
      0.21369201683658373846E-01,
      0.12020564476446432799E-01,
      0.53424751249537071952E-02,
      0.19232910450553508562E-02 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void debye3_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE3_VALUES returns some values of Debye's function of order 3.
//
//  Discussion:
//
//    The function is defined by:
//
//      DEBYE3(x) = 3 / x^3 * integral ( 0 <= t <= x ) t^3 / ( exp ( t ) - 1 ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.99926776885985461940E+00,
      0.98833007755734698212E+00,
      0.95390610472023510237E+00,
      0.82496296897623372315E+00,
      0.67441556407781468010E+00,
      0.54710665141286285468E+00,
      0.44112847372762418113E+00,
      0.35413603481042394211E+00,
      0.28357982814342246206E+00,
      0.18173691382177474795E+00,
      0.16277924385112436877E+00,
      0.11759741179993396450E+00,
      0.95240802723158889887E-01,
      0.77581324733763020269E-01,
      0.36560295673194845002E-01,
      0.19295765690345489563E-01,
      0.57712632276188798621E-02,
      0.24352200674805479827E-02,
      0.72154882216335666096E-03,
      0.15585454565440389896E-03 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void debye4_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DEBYE4_VALUES returns some values of Debye's function of order 4.
//
//  Discussion:
//
//    The function is defined by:
//
//      DEBYE4(x) = 4 / x^4 * integral ( 0 <= t <= x ) t^4 / ( exp ( t ) - 1 ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.99921896192761576256E+00,
      0.98755425280996071022E+00,
      0.95086788606389739976E+00,
      0.81384569172034042516E+00,
      0.65487406888673697092E+00,
      0.52162830964878715188E+00,
      0.41189273671788528876E+00,
      0.32295434858707304628E+00,
      0.25187863642883314410E+00,
      0.15185461258672022043E+00,
      0.13372661145921413299E+00,
      0.91471377664481164749E-01,
      0.71227828197462523663E-01,
      0.55676547822738862783E-01,
      0.21967566525574960096E-01,
      0.96736755602711590082E-02,
      0.19646978158351837850E-02,
      0.62214648623965450200E-03,
      0.12289514092077854510E-03,
      0.15927210319002161231E-04 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void dielectric_values ( int &n_data, double &tc, double &p, double &eps )

//****************************************************************************80
//
//  Purpose:
//
//    DIELECTRIC_VALUES returns some values of the static dielectric constant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, page 266.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &EPS, the dielectric constant, dimensionless.
//
{
# define N_MAX 15

  static double eps_vec[N_MAX] = {
      88.29E+00,
      90.07E+00,
      92.02E+00,
      95.14E+00,
     100.77E+00,
      78.85E+00,
      70.27E+00,
      62.60E+00,
      55.78E+00,
      44.31E+00,
      35.11E+00,
      20.40E+00,
       1.17E+00,
       1.11E+00,
       1.08E+00 };

  static double p_vec[N_MAX] = {
      100.0E+00,
      500.0E+00,
     1000.0E+00,
     2000.0E+00,
     5000.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00,
      100.0E+00 };

  static double tc_vec[N_MAX] = {
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
      25.0E+00,
      50.0E+00,
      75.0E+00,
     100.0E+00,
     150.0E+00,
     200.0E+00,
     300.0E+00,
     400.0E+00,
     500.0E+00,
     600.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    eps = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    eps = eps_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void dedekind_sum_values ( int &n_data, int &p, int &q, int &n, int &d )

//****************************************************************************80
//
//  Purpose:
//
//    DEDEKIND_SUM_VALUES returns some values of the Dedekind sum.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hans Rademacher, Emil Grosswald,
//    Dedekind Sums,
//    Mathematics Association of America, 1972,
//    LC: QA241.R2.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA
//    by 1, and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &P, &Q, the arguments of the function.
//
//    Output, int &N, &D, the numerator and denominator of the
//    function value.
//
{
# define N_MAX 95

  static int d_vec[N_MAX] = {
     1,  1, 18,  8,  5, 18, 14, 16, 27,  5,
    22, 72, 13, 14, 90, 32, 17, 27, 38, 40,
     1, 18,  1, 14, 27, 22, 13, 18, 17, 38,
     1,  1,  8,  1, 14, 16,  1, 22, 13, 14,
    32, 17, 38,  8,  1, 18,  5, 14, 27, 22,
    13, 90,  1, 38,  1,  1, 18,  8, 18, 14,
    16, 27, 22, 72,  1, 14, 32, 17, 27, 38,
     1,  5, 14, 22, 13, 17, 38,  1,  1, 18,
     8,  1, 18, 16, 27,  1, 22, 72, 13, 18,
    32, 17, 27, 38,  8 };
  static int n_vec[N_MAX] = {
     0,  0,  1,  1,  1,  5,  5,  7, 14,  3,
    15, 55, 11, 13, 91, 35, 20, 34, 51, 57,
     0, -1,  0,  1,  4,  5,  4,  7,  8, 21,
     0,  0, -1,  0, -1,  1,  0,  3,  1,  3,
     5,  5,  9,  3,  0,  1, -1,  1, -4,  3,
    -1, 19,  0, 11,  0,  0, -1,  1, -5, -1,
    -1,  4, -5, -1,  0,  3, -5,  1,  2, 11,
     0,  1, -5,  5, -4,  5, -9,  0,  0,  1,
    -1,  0,  5, -7, -4,  0, -3,  1,  4, -7,
    -3,  1, -2,  3,  3 };
  static int p_vec[N_MAX] = {
     1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
     2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
     3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
     3,  3,  3,  3,  4,  4,  4,  4,  4,  4,
     4,  4,  4,  4,  5,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     6,  6,  6,  6,  6,  6,  6,  7,  7,  7,
     7,  7,  7,  7,  7,  7,  7,  7,  7,  7,
     7,  7,  7,  7,  7 };
  static int q_vec[N_MAX] = {
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     1,  3,  5,  7,  9, 11, 13, 15, 17, 19,
     1,  2,  4,  5,  7,  8, 10, 11, 13, 14,
    16, 17, 19, 20,  1,  3,  5,  7,  9, 11,
    13, 15, 17, 19,  1,  2,  3,  4,  6,  7,
     8,  9, 11, 12, 13, 14, 16, 17, 18, 19,
     1,  5,  7, 11, 13, 17, 19,  1,  2,  3,
     4,  5,  6,  8,  9, 10, 11, 12, 13, 15,
    16, 17, 18, 19, 20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    p = 0;
    q = 0;
    n = 0;
    d = 0;
  }
  else
  {
    p = p_vec[n_data-1];
    q = q_vec[n_data-1];
    n = n_vec[n_data-1];
    d = d_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void dilogarithm_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    DILOGARITHM_VALUES returns some values of the dilogarithm function.
//
//  Discussion:
//
//    The dilogarithm is defined as
//
//      Li_2(X) = - integral ( 0 <= T <= X ) ln ( 1 - T ) / T dT
//
//    The dilogarithm is also known as Spence's integral.
//
//    In Abramowitz and Stegun form of the function is different,
//    and is equivalent to evaluated Li_2(1-X).
//
//    The dilogarithm is the special case, with N = 2, of the
//    polylogarithm Li_N(X).
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyLog[2,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Wolfram Media / Cambridge University Press, 1999.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.5063929246449603E-01,
     0.1026177910993911E+00,
     0.1560350339454831E+00,
     0.2110037754397048E+00,
     0.2676526390827326E+00,
     0.3261295100754761E+00,
     0.3866059411605865E+00,
     0.4492829744712817E+00,
     0.5143989891542119E+00,
     0.5822405264650125E+00,
     0.6531576315069018E+00,
     0.7275863077163334E+00,
     0.8060826895177240E+00,
     0.8893776242860387E+00,
     0.9784693929303061E+00,
     0.1074794600008248E+01,
     0.1180581123830255E+01,
     0.1299714723004959E+01,
     0.1440633796970039E+01,
     0.1644934066848226E+01 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.05E+00,
     0.10E+00,
     0.15E+00,
     0.20E+00,
     0.25E+00,
     0.30E+00,
     0.35E+00,
     0.40E+00,
     0.45E+00,
     0.50E+00,
     0.55E+00,
     0.60E+00,
     0.65E+00,
     0.70E+00,
     0.75E+00,
     0.80E+00,
     0.85E+00,
     0.90E+00,
     0.95E+00,
     0.10E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void e1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    E1_VALUES returns some values of the exponential integral function E1(X).
//
//  Definition:
//
//    The exponential integral E1(X) is defined by the formula:
//
//      E1(X) = integral ( 1 <= T <= +oo ) exp ( -X*T ) / T dT
//
//    In Mathematica, the function can be evaluated by:
//
//      ExpIntegralE[1,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.5597735947761608E+00,
     0.4543795031894021E+00,
     0.3737688432335091E+00,
     0.3105965785455430E+00,
     0.2601839393259996E+00,
     0.2193839343955203E+00,
     0.1859909045360402E+00,
     0.1584084368514626E+00,
     0.1354509578491291E+00,
     0.1162193125713579E+00,
     0.1000195824066327E+00,
     0.8630833369753979E-01,
     0.7465464440125305E-01,
     0.6471312936386886E-01,
     0.5620437817453485E-01,
     0.4890051070806112E-01 };

  static double x_vec[N_MAX] = {
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ei_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EI_VALUES returns some values of the exponential integral function EI(X).
//
//  Definition:
//
//    The exponential integral EI(X) has the formula:
//
//      EI(X) = - integral ( -X <= T < +oo ) exp ( -T ) / T dT
//
//    In Mathematica, the function can be evaluated by:
//
//      ExpIntegralEi[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.4542199048631736E+00,
     0.7698812899373594E+00,
     0.1064907194624291E+01,
     0.1347396548212326E+01,
     0.1622811713696867E+01,
     0.1895117816355937E+01,
     0.2167378279563403E+01,
     0.2442092285192652E+01,
     0.2721398880232024E+01,
     0.3007207464150646E+01,
     0.3301285449129798E+01,
     0.3605319949019469E+01,
     0.3920963201354904E+01,
     0.4249867557487934E+01,
     0.4593713686953585E+01,
     0.4954234356001890E+01 };

  static double x_vec[N_MAX] = {
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void elliptic_ea_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EA_VALUES returns values of the complete elliptic integral E(ALPHA).
//
//  Discussion:
//
//    This is one form of what is sometimes called the complete elliptic
//    integral of the second kind.
//
//    The function is defined by the formula:
//
//      E(ALPHA) = integral ( 0 <= T <= PI/2 )
//        sqrt ( 1 - sin ( ALPHA )^2 * sin ( T )^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      EllipticE[(Sin[Pi*alpha/180])^2]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function, measured
//    in degrees.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 19

  static double fx_vec[N_MAX] = {
     1.570796326794897E+00,
     1.567809073977622E+00,
     1.558887196601596E+00,
     1.544150496914673E+00,
     1.523799205259774E+00,
     1.498114928422116E+00,
     1.467462209339427E+00,
     1.432290969306756E+00,
     1.393140248523812E+00,
     1.350643881047676E+00,
     1.305539094297794E+00,
     1.258679624779997E+00,
     1.211056027568459E+00,
     1.163827964493139E+00,
     1.118377737969864E+00,
     1.076405113076403E+00,
     1.040114395706010E+00,
     1.012663506234396E+00,
     1.000000000000000E+00 };

  static double x_vec[N_MAX] = {
      0.0E+00,
      5.0E+00,
     10.0E+00,
     15.0E+00,
     20.0E+00,
     25.0E+00,
     30.0E+00,
     35.0E+00,
     40.0E+00,
     45.0E+00,
     50.0E+00,
     55.0E+00,
     60.0E+00,
     65.0E+00,
     70.0E+00,
     75.0E+00,
     80.0E+00,
     85.0E+00,
     90.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void elliptic_em_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_EM_VALUES returns values of the complete elliptic integral E(M).
//
//  Discussion:
//
//    This is one form of what is sometimes called the complete elliptic
//    integral of the second kind.
//
//    The function is defined by the formula:
//
//      E(M) = integral ( 0 <= T <= PI/2 )
//        sqrt ( 1 - M * sin ( T )^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      EllipticE[m]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     1.570796326794897E+00,
     1.550973351780472E+00,
     1.530757636897763E+00,
     1.510121832092819E+00,
     1.489035058095853E+00,
     1.467462209339427E+00,
     1.445363064412665E+00,
     1.422691133490879E+00,
     1.399392138897432E+00,
     1.375401971871116E+00,
     1.350643881047676E+00,
     1.325024497958230E+00,
     1.298428035046913E+00,
     1.270707479650149E+00,
     1.241670567945823E+00,
     1.211056027568459E+00,
     1.178489924327839E+00,
     1.143395791883166E+00,
     1.104774732704073E+00,
     1.060473727766278E+00,
     1.000000000000000E+00 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.05E+00,
     0.10E+00,
     0.15E+00,
     0.20E+00,
     0.25E+00,
     0.30E+00,
     0.35E+00,
     0.40E+00,
     0.45E+00,
     0.50E+00,
     0.55E+00,
     0.60E+00,
     0.65E+00,
     0.70E+00,
     0.75E+00,
     0.80E+00,
     0.85E+00,
     0.90E+00,
     0.95E+00,
     1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void elliptic_ka_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_KA_VALUES returns values of the complete elliptic integral K(ALPHA).
//
//  Discussion:
//
//    This is one form of what is sometimes called the complete elliptic integral
//    of the first kind.
//
//    The function is defined by the formula:
//
//      K(ALPHA) = integral ( 0 <= T <= PI/2 )
//        dT / sqrt ( 1 - sin ( ALPHA )^2 * sin ( T )^2 )
//
//    In Mathematica, the function can be evaluated by:
//
//      EllipticK[(Sin[alpha*Pi/180])^2]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function, measured
//    in degrees.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
     0.1570796326794897E+01,
     0.1573792130924768E+01,
     0.1582842804338351E+01,
     0.1598142002112540E+01,
     0.1620025899124204E+01,
     0.1648995218478530E+01,
     0.1685750354812596E+01,
     0.1731245175657058E+01,
     0.1786769134885021E+01,
     0.1854074677301372E+01,
     0.1935581096004722E+01,
     0.2034715312185791E+01,
     0.2156515647499643E+01,
     0.2308786798167196E+01,
     0.2504550079001634E+01,
     0.2768063145368768E+01,
     0.3153385251887839E+01,
     0.3831741999784146E+01 };

  static double x_vec[N_MAX] = {
      0.0E+00,
      5.0E+00,
     10.0E+00,
     15.0E+00,
     20.0E+00,
     25.0E+00,
     30.0E+00,
     35.0E+00,
     40.0E+00,
     45.0E+00,
     50.0E+00,
     55.0E+00,
     60.0E+00,
     65.0E+00,
     70.0E+00,
     75.0E+00,
     80.0E+00,
     85.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void elliptic_km_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ELLIPTIC_KM_VALUES returns values of the complete elliptic integral K(M).
//
//  Discussion:
//
//    This is one form of what is sometimes called the complete elliptic
//    integral of the first kind.
//
//    The function is defined by the formula:
//
//      K(M) = integral ( 0 <= T <= PI/2 )
//        dT / sqrt ( 1 - M * sin ( T )^2 )
//
//    In Mathematica, the function can be evaluated by:
//
//      EllipticK[m]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     1.570796326794897E+00,
     1.591003453790792E+00,
     1.612441348720219E+00,
     1.635256732264580E+00,
     1.659623598610528E+00,
     1.685750354812596E+00,
     1.713889448178791E+00,
     1.744350597225613E+00,
     1.777519371491253E+00,
     1.813883936816983E+00,
     1.854074677301372E+00,
     1.898924910271554E+00,
     1.949567749806026E+00,
     2.007598398424376E+00,
     2.075363135292469E+00,
     2.156515647499643E+00,
     2.257205326820854E+00,
     2.389016486325580E+00,
     2.578092113348173E+00,
     2.908337248444552E+00 };

  static double x_vec[N_MAX] = {
    0.00E+00,
    0.05E+00,
    0.10E+00,
    0.15E+00,
    0.20E+00,
    0.25E+00,
    0.30E+00,
    0.35E+00,
    0.40E+00,
    0.45E+00,
    0.50E+00,
    0.55E+00,
    0.60E+00,
    0.65E+00,
    0.70E+00,
    0.75E+00,
    0.80E+00,
    0.85E+00,
    0.90E+00,
    0.95E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void erf_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERF_VALUES returns some values of the ERF or "error" function.
//
//  Discussion:
//
//    The error function is defined by:
//
//      ERF(X) = ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      Erf[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.1124629160182849E+00,
     0.2227025892104785E+00,
     0.3286267594591274E+00,
     0.4283923550466685E+00,
     0.5204998778130465E+00,
     0.6038560908479259E+00,
     0.6778011938374185E+00,
     0.7421009647076605E+00,
     0.7969082124228321E+00,
     0.8427007929497149E+00,
     0.8802050695740817E+00,
     0.9103139782296354E+00,
     0.9340079449406524E+00,
     0.9522851197626488E+00,
     0.9661051464753107E+00,
     0.9763483833446440E+00,
     0.9837904585907746E+00,
     0.9890905016357307E+00,
     0.9927904292352575E+00,
     0.9953222650189527E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.1E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.6E+00,
     0.7E+00,
     0.8E+00,
     0.9E+00,
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void erfc_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    ERFC_VALUES returns some values of the ERFC or "complementary error" function.
//
//  Discussion:
//
//    The complementary error function is defined by:
//
//      ERFC(X) = 1 - ( 2 / sqrt ( PI ) * integral ( 0 <= T <= X ) exp ( - T^2 ) dT
//
//    In Mathematica, the function can be evaluated by:
//
//      Erfc[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    1.000000000000000E+00,
    0.7772974107895215E+00,
    0.5716076449533315E+00,
    0.3961439091520741E+00,
    0.2578990352923395E+00,
    0.1572992070502851E+00,
    0.08968602177036462E+00,
    0.04771488023735119E+00,
    0.02365161665535599E+00,
    0.01090949836426929E+00,
    0.004677734981047266E+00,
    0.001862846297981891E+00,
    0.0006885138966450786E+00,
    0.0002360344165293492E+00,
    0.00007501319466545902E+00,
    0.00002209049699858544E+00,
    6.025761151762095E-06,
    1.521993362862285E-06,
    3.558629930076853E-07,
    7.700392745696413E-08,
    1.541725790028002E-08 };

  static double x_vec[N_MAX] = {
    0.0E+00,
    0.2E+00,
    0.4E+00,
    0.6E+00,
    0.8E+00,
    1.0E+00,
    1.2E+00,
    1.4E+00,
    1.6E+00,
    1.8E+00,
    2.0E+00,
    2.2E+00,
    2.4E+00,
    2.6E+00,
    2.8E+00,
    3.0E+00,
    3.2E+00,
    3.4E+00,
    3.6E+00,
    3.8E+00,
    4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void euler_number_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_NUMBER_VALUES returns some values of the Euler numbers.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      EulerE[n]
//
//    These numbers rapidly get too big to store in an ordinary integer!
//
//    The terms of odd index are 0.
//
//    E(N) = -C(N,N-2) * E(N-2) - C(N,N-4) * E(N-4) - ... - C(N,0) * E(0).
//
//  First terms:
//
//    E0  = 1
//    E1  = 0;
//    E2  = -1
//    E3  = 0;
//    E4  = 5
//    E5  = 0;
//    E6  = -61
//    E7  = 0;
//    E8  = 1385
//    E9  = 0;
//    E10 = -50521
//    E11 = 0;
//    E12 = 2702765
//    E13 = 0;
//    E14 = -199360981
//    E15 = 0;
//    E16 = 19391512145
//    E17 = 0;
//    E18 = -2404879675441
//    E19 = 0;
//    E20 = 370371188237525
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Euler number.
//
//    Output, int &C, the value of the Euler number.
//
{
# define N_MAX 8

  static int c_vec[N_MAX] = {
    1, 0, -1, 5, 61, 1385, -50521, 2702765 };

  static int n_vec[N_MAX] = {
     0, 1, 2, 4, 6, 8, 10, 12 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void euler_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EULER_POLY_VALUES returns some values of the Euler polynomials.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      EulerE[n,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the Euler polynomial.
//
//    Output, double &X, the argument of the Euler polynomial.
//
//    Output, double &FX, the value of the Euler polynomial.
//
{
# define N_MAX 27

  static double fx_vec[N_MAX] = {
      0.100000000000E+01,
     -0.300000000000E+00,
     -0.160000000000E+00,
      0.198000000000E+00,
      0.185600000000E+00,
     -0.403680000000E+00,
     -0.560896000000E+00,
      0.171878880000E+01,
      0.318043136000E+01,
     -0.125394670080E+02,
     -0.289999384576E+02,
     -0.625000000000E-01,
     -0.174240000000E+00,
     -0.297680000000E+00,
     -0.404320000000E+00,
     -0.475260000000E+00,
     -0.500000000000E+00,
     -0.475240000000E+00,
     -0.403680000000E+00,
     -0.292820000000E+00,
     -0.153760000000E+00,
      0.000000000000E+00,
      0.153760000000E+00,
      0.292820000000E+00,
      0.403680000000E+00,
      0.475240000000E+00,
      0.500000000000E+00 };

  static int n_vec[N_MAX] = {
     0,
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5 };

  static double x_vec[N_MAX] = {
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
      0.2E+00,
     -0.5E+00,
     -0.4E+00,
     -0.3E+00,
     -0.2E+00,
     -0.1E+00,
      0.0E+00,
      0.1E+00,
      0.2E+00,
      0.3E+00,
      0.4E+00,
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void exp_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXP_VALUES returns some values of the exponential function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Exp[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 March 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 24

  static double fx_vec[N_MAX] = {
    0.000045399929762484851536E+00,
    0.0067379469990854670966E+00,
    0.36787944117144232160E+00,
    1.0000000000000000000E+00,
    1.0000000100000000500E+00,
    1.0001000050001666708E+00,
    1.0010005001667083417E+00,
    1.0100501670841680575E+00,
    1.1051709180756476248E+00,
    1.2214027581601698339E+00,
    1.3498588075760031040E+00,
    1.4918246976412703178E+00,
    1.6487212707001281468E+00,
    1.8221188003905089749E+00,
    2.0137527074704765216E+00,
    2.2255409284924676046E+00,
    2.4596031111569496638E+00,
    2.7182818284590452354E+00,
    7.3890560989306502272E+00,
    23.140692632779269006E+00,
    148.41315910257660342E+00,
    22026.465794806716517E+00,
    4.8516519540979027797E+08,
    2.3538526683701998541E+17 };

  static double x_vec[N_MAX] = {
     -10.0E+00,
      -5.0E+00,
      -1.0E+00,
       0.0E+00,
       0.00000001E+00,
       0.0001E+00,
       0.001E+00,
       0.01E+00,
       0.1E+00,
       0.2E+00,
       0.3E+00,
       0.4E+00,
       0.5E+00,
       0.6E+00,
       0.7E+00,
       0.8E+00,
       0.9E+00,
       1.0E+00,
       2.0E+00,
       3.1415926535897932385E+00,
       5.0E+00,
      10.0E+00,
      20.0E+00,
      40.0E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void exp3_int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXP3_INT_VALUES returns some values of the EXP3 integral function.
//
//  Discussion:
//
//    The function is defined by:
//
//      EXP3_INT(x) = Integral ( 0 <= t <= x ) exp ( -t^3 ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.19531249963620212007E-02,
      0.78124990686775522671E-02,
      0.31249761583499728667E-01,
      0.12493899888803079984E+00,
      0.48491714311363971332E+00,
      0.80751118213967145286E+00,
      0.86889265412623270696E+00,
      0.88861722235357162648E+00,
      0.89286018500218176869E+00,
      0.89295351429387631138E+00,
      0.89297479112737843939E+00,
      0.89297880579798112220E+00,
      0.89297950317496621294E+00,
      0.89297951152951902903E+00,
      0.89297951156918122102E+00,
      0.89297951156924734716E+00,
      0.89297951156924917298E+00,
      0.89297951156924921121E+00,
      0.89297951156924921122E+00,
      0.89297951156924921122E+00 };

  static double x_vec[N_MAX] = {
      0.0019531250E+00,
      0.0078125000E+00,
      0.0312500000E+00,
      0.1250000000E+00,
      0.5000000000E+00,
      1.0000000000E+00,
      1.2500000000E+00,
      1.5000000000E+00,
      1.8750000000E+00,
      2.0000000000E+00,
      2.1250000000E+00,
      2.2500000000E+00,
      2.5000000000E+00,
      2.7500000000E+00,
      3.0000000000E+00,
      3.1250000000E+00,
      3.2500000000E+00,
      3.5000000000E+00,
      3.7500000000E+00,
      4.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void exponential_cdf_values ( int &n_data, double &lambda, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ExponentialDistribution [ lambda ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &LAMBDA, the parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 9

  static double fx_vec[N_MAX] = {
     0.3934693402873666E+00,
     0.6321205588285577E+00,
     0.7768698398515702E+00,
     0.8646647167633873E+00,
     0.8646647167633873E+00,
     0.9816843611112658E+00,
     0.9975212478233336E+00,
     0.9996645373720975E+00,
     0.9999546000702375E+00 };

  static double lambda_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    lambda = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    lambda = lambda_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void extreme_values_cdf_values ( int &n_data, double &alpha, double &beta,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    EXTREME_VALUES_CDF_VALUES returns some values of the Extreme Values CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = ExtremeValuesDistribution [ alpha, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &ALPHA, the first parameter of the distribution.
//
//    Output, double &BETA, the second parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double alpha_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double beta_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double fx_vec[N_MAX] = {
     0.3678794411714423E+00,
     0.8734230184931166E+00,
     0.9818510730616665E+00,
     0.9975243173927525E+00,
     0.5452392118926051E+00,
     0.4884435800065159E+00,
     0.4589560693076638E+00,
     0.4409910259429826E+00,
     0.5452392118926051E+00,
     0.3678794411714423E+00,
     0.1922956455479649E+00,
     0.6598803584531254E-01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    alpha = 0.0;
    beta = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    alpha = alpha_vec[n_data-1];
    beta = beta_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void f_cdf_values ( int &n_data, int &a, int &b, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = FRatioDistribution [ dfn, dfd ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &A, int &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static int a_vec[N_MAX] = {
    1,
    1,
    5,
    1,
    2,
    4,
    1,
    6,
    8,
    1,
    3,
    6,
    1,
    1,
    1,
    1,
    2,
    3,
    4,
    5 };

  static int b_vec[N_MAX] = {
     1,
     5,
     1,
     5,
    10,
    20,
     5,
     6,
    16,
     5,
    10,
    12,
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     5 };

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.4999714850534485E+00,
     0.4996034370170990E+00,
     0.7496993658293228E+00,
     0.7504656462757382E+00,
     0.7514156325324275E+00,
     0.8999867031372156E+00,
     0.8997127554259699E+00,
     0.9002845660853669E+00,
     0.9500248817817622E+00,
     0.9500574946122442E+00,
     0.9501926400000000E+00,
     0.9750133887312993E+00,
     0.9900022327445249E+00,
     0.9949977837872073E+00,
     0.9989999621122122E+00,
     0.5687988496283079E+00,
     0.5351452100063650E+00,
     0.5143428032407864E+00,
     0.5000000000000000E+00 };

  static double x_vec[N_MAX] = {
      1.00E+00,
      0.528E+00,
      1.89E+00,
      1.69E+00,
      1.60E+00,
      1.47E+00,
      4.06E+00,
      3.05E+00,
      2.09E+00,
      6.61E+00,
      3.71E+00,
      3.00E+00,
     10.01E+00,
     16.26E+00,
     22.78E+00,
     47.18E+00,
      1.00E+00,
      1.00E+00,
      1.00E+00,
      1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0;
    b = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void f_noncentral_cdf_values ( int &n_data, int &n1, int &n2, double &lambda,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralFRatioDistribution [ n1, n2, lambda ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N1, int &N2, the numerator and denominator
//    degrees of freedom.
//
//    Output, double &LAMBDA, the noncentrality parameter.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.6367825323508774E+00,
     0.5840916116305482E+00,
     0.3234431872392788E+00,
     0.4501187879813550E+00,
     0.6078881441188312E+00,
     0.7059275551414605E+00,
     0.7721782003263727E+00,
     0.8191049017635072E+00,
     0.3170348430749965E+00,
     0.4327218008454471E+00,
     0.4502696915707327E+00,
     0.4261881186594096E+00,
     0.6753687206341544E+00,
     0.4229108778879005E+00,
     0.6927667261228938E+00,
     0.3632174676491226E+00,
     0.4210054012695865E+00,
     0.4266672258818927E+00,
     0.4464016600524644E+00,
     0.8445888579504827E+00,
     0.4339300273343604E+00 };

  static double lambda_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.25E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     1.00E+00,
     1.00E+00,
     0.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00 };

  static int n1_vec[N_MAX] = {
     1,  1,  1,  1,
     1,  1,  1,  1,
     1,  1,  2,  2,
     3,  3,  4,  4,
     5,  5,  6,  6,
     8, 16 };

  static int n2_vec[N_MAX] = {
     1,  5,  5,  5,
     5,  5,  5,  5,
     5,  5,  5, 10,
     5,  5,  5,  5,
     1,  5,  6, 12,
    16,  8 };

  static double x_vec[N_MAX] = {
     1.00E+00,
     1.00E+00,
     1.00E+00,
     0.50E+00,
     1.00E+00,
     2.00E+00,
     3.00E+00,
     4.00E+00,
     5.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     2.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n1 = 0;
    n2 = 0;
    lambda = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n1 = n1_vec[n_data-1];
    n2 = n2_vec[n_data-1];
    lambda = lambda_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void factorial_values ( int &n_data, int &n, int &fn )

//****************************************************************************80
//
//  Purpose:
//
//    FACTORIAL_VALUES returns values of the factorial function.
//
//  Discussion:
//
//    0! = 1
//    I! = Product ( 1 <= J <= I ) I
//
//    In Mathematica, the function can be evaluated by:
//
//      n!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the function.
//
//    Output, int &FN, the value of the function.
//
{
# define N_MAX 13

  static int fn_vec[N_MAX] = {
            1,
            1,
            2,
            6,
           24,
          120,
          720,
         5040,
        40320,
       362880,
      3628800,
     39916800,
    479001600 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void factorial2_values ( int &n_data, int &n, int &fn )

//****************************************************************************80
//
//  Purpose:
//
//    FACTORIAL2_VALUES returns values of the double factorial function.
//
//  Formula:
//
//    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
//                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
//
//    In Mathematica, the function can be evaluated by:
//
//      n!!
//
//  Example:
//
//     N    N!!
//
//     0     1
//     1     1
//     2     2
//     3     3
//     4     8
//     5    15
//     6    48
//     7   105
//     8   384
//     9   945
//    10  3840
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, page 16.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the function.
//
//    Output, int &FN, the value of the function.
//
{
# define N_MAX 16

  static int fn_vec[N_MAX] = {
          1,
          1,
          2,
          3,
          8,
         15,
         48,
        105,
        384,
        945,
       3840,
      10395,
      46080,
     135135,
     645120,
    2027025 };

  static int n_vec[N_MAX] = {
     0,
     1,  2,  3,  4,  5,
     6,  7,  8,  9, 10,
    11, 12, 13, 14, 15 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void factorial_rising_values ( int &n_data, int &m, int &n, int &fmn )

//****************************************************************************80
//
//  Purpose:
//
//    FACTORIAL_RISING_VALUES returns values of the int *Pochhammer function.
//
//  Discussion:
//
//    The int *Pochhammer function is sometimes symbolized by (m)_n.
//
//    The definition of the int *Pochhammer function is
//
//      (m)_n = (m-1+n)! / (m-1)!
//            = ( m ) * ( m + 1 ) * ( m + 2 ) ... * ( m - 1 + n )
//            = Gamma ( m + n ) / Gamma ( m )
//
//    We assume 0 <= N <= M.
//
//    In Mathematica, the function can be evaluated by:
//
//      Pochhammer[m,n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &M, &N, the arguments of the function.
//
//    Output, int &FMN, the value of the function.
//
{
# define N_MAX 8

  static int fmn_vec[N_MAX] = {
     1, 10, 4000, 110, 6840, 840, 970200, 5040 };

  static int m_vec[N_MAX] = {
    50, 10, 4000, 10, 18, 4, 98, 1 };

  static int n_vec[N_MAX] = {
    0,  1,   1,   2,  3, 4,  3, 7 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    m = 0;
    n = 0;
    fmn = 0;
  }
  else
  {
    m = m_vec[n_data-1];
    n = n_vec[n_data-1];
    fmn = fmn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void fresnel_cos_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    FRESNEL_COS_VALUES returns values of the Fresnel cosine integral function.
//
//  Discussion:
//
//    The Fresnel cosine integral is defined by:
//
//      C(X) = integral ( 0 <= T <= X ) cos ( PI * T^2 / 2 ) / T dT
//
//    In Mathematica, the function can be evaluated by:
//
//      FresnelC[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.1999210575944531E+00,
     0.3974807591723594E+00,
     0.5810954469916523E+00,
     0.7228441718963561E+00,
     0.7798934003768228E+00,
     0.7154377229230734E+00,
     0.5430957835462564E+00,
     0.3654616834404877E+00,
     0.3336329272215571E+00,
     0.4882534060753408E+00,
     0.6362860449033195E+00,
     0.5549614058564281E+00,
     0.3889374961919690E+00,
     0.4674916516989059E+00,
     0.6057207892976856E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void fresnel_sin_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    FRESNEL_SIN_VALUES returns some values of the Fresnel sine integral function.
//
//  Discussion:
//
//    The Fresnel sine integral is defined by
//
//      S(X) = integral ( 0 <= T <= X ) sin ( pi * T^2 / 2 ) / T dT
//
//    In Mathematica, the function can be evaluated by:
//
//      FresnelS[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 April 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.4187609161656762E-02,
     0.3335943266061318E-01,
     0.1105402073593870E+00,
     0.2493413930539178E+00,
     0.4382591473903548E+00,
     0.6234009185462497E+00,
     0.7135250773634121E+00,
     0.6388876835093809E+00,
     0.4509387692675831E+00,
     0.3434156783636982E+00,
     0.4557046121246569E+00,
     0.6196899649456836E+00,
     0.5499893231527195E+00,
     0.3915284435431718E+00,
     0.4963129989673750E+00 };

  static double x_vec[N_MAX] = {
     0.0E+00,
     0.2E+00,
     0.4E+00,
     0.6E+00,
     0.8E+00,
     1.0E+00,
     1.2E+00,
     1.4E+00,
     1.6E+00,
     1.8E+00,
     2.0E+00,
     2.2E+00,
     2.4E+00,
     2.6E+00,
     2.8E+00,
     3.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void frobenius_number_data_values ( int &n_data, int order, int c[], int &f )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_DATA_VALUES returns data for the Frobenius problem.
//
//  Discussion:
//
//    The user should first call FROBENIUS_NUMBER_ORDER_VALUES to get the
//    order or size of the "C" vector that will be returned by this routine.
//
//    The Frobenius number of order N and data C is the solution of the
//    Frobenius coin sum problem for N coin denominations C(1) through C(N).
//
//    The Frobenius coin sum problem assumes the existence of
//    N coin denominations, and asks for the largest value that cannot
//    be formed by any combination of coins of these denominations.
//
//    The coin denominations are assumed to be distinct positive integers.
//
//    For general N, this problem is fairly difficult to handle.
//
//    For N = 2, it is known that:
//
//    * if C1 and C2 are not relatively prime, then
//      there are infinitely large values that cannot be formed.
//
//    * otherwise, the largest value that cannot be formed is
//      C1 * C2 - C1 - C2, and that exactly half the values between
//      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
//
//    As a simple example, if C1 = 2 and C2 = 7, then the largest
//    unrepresentable value is 5, and there are (5+1)/2 = 3
//    unrepresentable values, namely 1, 3, and 5.
//
//    For a general N, and a set of coin denominations C1, C2, ..., CN,
//    the Frobenius number F(N, C(1:N) ) is defined as the largest value
//    B for which the equation
//
//      C1*X1 + C2*X2 + ... + CN*XN = B
//
//    has no nonnegative integer solution X(1:N).
//
//    In Mathematica, the Frobenius number can be determined by
//
//      FrobeniusNumber[ {C1,...,CN} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
//    Decomposition of Integer Programs and of Generating Sets,
//    in Algorithms - ESA 97,
//    Lecture Notes in Computer Science 1284,
//    edited by R Burkard, G Woeginger,
//    Springer, 1997, pages 92-103.
//
//    Robert Owens,
//    An Algorithm to Solve the Frobenius Problem,
//    Mathematics Magazine,
//    Volume 76, Number 4, October 2003, 264-275.
//
//    James Sylvester,
//    Question 7382,
//    Mathematical Questions with their Solutions,
//    Educational Times,
//    Volume 41, page 21, 1884.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input, int &N_DATA.  Unlike most other routines in this
//    library, this routine assumes that N_DATA has already been processed by a call
//    to FROBENIUS_NUMBER_ORDER_VALUE.  Therefore, this routine will return the
//    next set of data as long as N_DATA is in the expected range.
//
//    Input, int ORDER, the order of the problem.
//
//    Output, int C[ORDER], the denominations of the problem.
//
//    Output, int &F, the value of the function.
//
{
# define CVEC_MAX 77
# define N_MAX 19

  static int c_vec[CVEC_MAX] = {
        2,     5,
        3,    17,
        4,    19,
        5,    13,
       12,    11,
       99,   100,
        6,     9,    20,
        5,    17,    23,
      137,   251,   256,
       31,    41,    47,    61,
      271,   277,   281,   283,
       10,    18,    26,    33,    35,
       34,    37,    38,    40,    43,
    12223, 12224, 36674, 61119, 85569,
     1000,  1476,  3764,  4864,  4871,  7773,
    12228, 36679, 36682, 46908, 61139, 73365,
    12137, 36405, 24269, 36407, 84545, 60683,
    13211, 13212, 39638, 66060, 52864, 79268, 92482,
    13429, 26850, 26855, 40280, 40281, 53711, 53714, 67141
  };
  static int f_vec[N_MAX] = {
          3,
         31,
         53,
         47,
        109,
       9701,
         43,
         41,
       4948,
        240,
      13022,
         67,
        175,
   89643481,
      47350,
   89716838,
   58925134,
  104723595,
   45094583 };
  int i;
  static int v_data = 0;

  if ( n_data < 1 || N_MAX < n_data )
  {
    n_data = 0;
    v_data = 0;
    for ( i = 0; i < order; i++ )
    {
      c[i] = 0;
    }
    f = 0;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      c[i] = c_vec[v_data+i];
    }
    v_data = v_data + order;
    if ( n_data == N_MAX )
    {
      v_data = 0;
    }
    f = f_vec[n_data-1];
  }

  return;
# undef CVEC_MAX
# undef N_MAX
}
//****************************************************************************80

void frobenius_number_order_values ( int &n_data, int &order )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER_VALUES returns orders of the Frobenius problem.
//
//  Discussion:
//
//    This routine returns the order or size of a Frobenius problem.
//    To get the corresponding data, call FROBENIUS_NUMBER_DATA_VALUES.
//
//    The Frobenius number of order N and data C is the solution of a Frobenius
//    coin sum problem for N coin denominations C(1) through C(N).
//
//    The Frobenius coin sum problem assumes the existence of
//    N coin denominations, and asks for the largest value that cannot
//    be formed by any combination of coins of these denominations.
//
//    The coin denominations are assumed to be distinct positive integers.
//
//    For general order N, this problem is fairly difficult to handle.
//
//    For order N = 2, it is known that:
//
//    * if C1 and C2 are not relatively prime, then
//      there are infinitely large values that cannot be formed.
//
//    * otherwise, the largest value that cannot be formed is
//      C1 * C2 - C1 - C2, and that exactly half the values between
//      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
//
//    As a simple example, if C1 = 2 and C2 = 7, then the largest
//    unrepresentable value is 5, and there are (5+1)/2 = 3
//    unrepresentable values, namely 1, 3, and 5.
//
//    For a general N, and a set of coin denominations C1, C2, ..., CN,
//    the Frobenius number F(N, C(1:N) ) is defined as the largest value
//    B for which the equation
//
//      C1*X1 + C2*X2 + ... + CN*XN = B
//
//    has no nonnegative integer solution X(1:N).
//
//    In Mathematica, the Frobenius number can be determined by
//
//      FrobeniusNumber[ {C1,...,CN} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Gerard Cornuejols, Regina Urbaniak, Robert Weismantel, Laurence Wolsey,
//    Decomposition of Integer Programs and of Generating Sets,
//    in Algorithms - ESA 97,
//    Lecture Notes in Computer Science 1284,
//    edited by R Burkard, G Woeginger,
//    Springer, 1997, pages 92-103,
//    LC: QA76.9.A43.E83.
//
//    Robert Owens,
//    An Algorithm to Solve the Frobenius Problem,
//    Mathematics Magazine,
//    Volume 76, Number 4, October 2003, 264-275.
//
//    James Sylvester,
//    Question 7382,
//    Mathematical Questions with their Solutions,
//    Educational Times,
//    Volume 41, page 21, 1884.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &ORDER, the order of a Frobenius problem.
//
{
# define N_MAX 19

  static int order_vec[N_MAX] = {
    2,
    2,
    2,
    2,
    2,
    2,
    3,
    3,
    3,
    4,
    4,
    5,
    5,
    5,
    6,
    6,
    6,
    7,
    8
  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    order = 0;
  }
  else
  {
    order = order_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void frobenius_number_order2_values ( int &n_data, int &c1, int &c2, int &f )

//****************************************************************************80
//
//  Purpose:
//
//    FROBENIUS_NUMBER_ORDER2_VALUES returns values of the order 2 Frobenius number.
//
//  Discussion:
//
//    The Frobenius number of order N is the solution of the Frobenius
//    coin sum problem for N coin denominations.
//
//    The Frobenius coin sum problem assumes the existence of
//    N coin denominations, and asks for the largest value that cannot
//    be formed by any combination of coins of these denominations.
//
//    The coin denominations are assumed to be distinct positive integers.
//
//    For general N, this problem is fairly difficult to handle.
//
//    For N = 2, it is known that:
//
//    * if C1 and C2 are not relatively prime, then
//      there are infinitely large values that cannot be formed.
//
//    * otherwise, the largest value that cannot be formed is
//      C1 * C2 - C1 - C2, and that exactly half the values between
//      1 and C1 * C2 - C1 - C2 + 1 cannot be represented.
//
//    As a simple example, if C1 = 2 and C2 = 7, then the largest
//    unrepresentable value is 5, and there are (5+1)/2 = 3
//    unrepresentable values, namely 1, 3, and 5.
//
//    For a general N, and a set of coin denominations C1, C2, ..., CN,
//    the Frobenius number F(N, C(1:N) ) is defined as the largest value
//    B for which the equation
//
//      C1*X1 + C2*X2 + ... + CN*XN = B
//
//    has no nonnegative integer solution X(1:N).
//
//    In Mathematica, the Frobenius number can be determined by
//
//      FrobeniusNumber[ {C1,...,CN} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    James Sylvester,
//    Question 7382,
//    Mathematical Questions with their Solutions,
//    Educational Times,
//    Volume 41, page 21, 1884.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &C1, &C2, the parameters of the function.
//
//    Output, int &F, the value of the function.
//
{
# define N_MAX 6

  static int c1_vec[N_MAX] = {
     2,
     3,
     4,
     5,
    12,
    99 };
  static int c2_vec[N_MAX] = {
      5,
     17,
     19,
     13,
     11,
    100 };
  static int f_vec[N_MAX] = {
    3,
    31,
    53,
    47,
    109,
    9701 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    c1 = 0;
    c2 = 0;
    f = 0;
  }
  else
  {
    c1 = c1_vec[n_data-1];
    c2 = c2_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_VALUES returns some values of the Gamma function.
//
//  Discussion:
//
//    The Gamma function is defined as:
//
//      Gamma(Z) = Integral ( 0 <= T < +oo ) T^(Z-1) exp(-T) dT
//
//    It satisfies the recursion:
//
//      Gamma(X+1) = X * Gamma(X)
//
//    Gamma is undefined for nonpositive integral X.
//    Gamma(0.5) = sqrt(PI)
//    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
//
//    In Mathematica, the function can be evaluated by:
//
//      Gamma[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 25

  static double fx_vec[N_MAX] = {
     -0.3544907701811032E+01,
     -0.1005871979644108E+03,
      0.9943258511915060E+02,
      0.9513507698668732E+01,
      0.4590843711998803E+01,
      0.2218159543757688E+01,
      0.1772453850905516E+01,
      0.1489192248812817E+01,
      0.1164229713725303E+01,
      0.1000000000000000E+01,
      0.9513507698668732E+00,
      0.9181687423997606E+00,
      0.8974706963062772E+00,
      0.8872638175030753E+00,
      0.8862269254527580E+00,
      0.8935153492876903E+00,
      0.9086387328532904E+00,
      0.9313837709802427E+00,
      0.9617658319073874E+00,
      0.1000000000000000E+01,
      0.2000000000000000E+01,
      0.6000000000000000E+01,
      0.3628800000000000E+06,
      0.1216451004088320E+18,
      0.8841761993739702E+31 };

  static double x_vec[N_MAX] = {
     -0.50E+00,
     -0.01E+00,
      0.01E+00,
      0.10E+00,
      0.20E+00,
      0.40E+00,
      0.50E+00,
      0.60E+00,
      0.80E+00,
      1.00E+00,
      1.10E+00,
      1.20E+00,
      1.30E+00,
      1.40E+00,
      1.50E+00,
      1.60E+00,
      1.70E+00,
      1.80E+00,
      1.90E+00,
      2.00E+00,
      3.00E+00,
      4.00E+00,
     10.00E+00,
     20.00E+00,
     30.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_cdf_values ( int &n_data, double &mu, double &sigma, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_CDF_VALUES returns some values of the Gamma CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = GammaDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &SIGMA, the variance of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.8646647167633873E+00,
     0.9816843611112658E+00,
     0.9975212478233336E+00,
     0.9996645373720975E+00,
     0.6321205588285577E+00,
     0.4865828809674080E+00,
     0.3934693402873666E+00,
     0.3296799539643607E+00,
     0.4421745996289254E+00,
     0.1911531694619419E+00,
     0.6564245437845009E-01,
     0.1857593622214067E-01 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_inc_p_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_P_VALUES: values of the normalized incomplete Gamma function P(A,X).
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function is defined as:
//
//      P(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= P(A,X) <= 1
//
//    and
//
//      P(A,oo) = 1.0
//
//    In Mathematica, the function can be evaluated by:
//
//      1 - GammaRegularized[A,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.10E+01,
     0.10E+01,
     0.10E+01,
     0.11E+01,
     0.11E+01,
     0.11E+01,
     0.20E+01,
     0.20E+01,
     0.20E+01,
     0.60E+01,
     0.60E+01,
     0.11E+02,
     0.26E+02,
     0.41E+02  };

  static double fx_vec[N_MAX] = {
     0.7382350532339351E+00,
     0.9083579897300343E+00,
     0.9886559833621947E+00,
     0.3014646416966613E+00,
     0.7793286380801532E+00,
     0.9918490284064973E+00,
     0.9516258196404043E-01,
     0.6321205588285577E+00,
     0.9932620530009145E+00,
     0.7205974576054322E-01,
     0.5891809618706485E+00,
     0.9915368159845525E+00,
     0.1018582711118352E-01,
     0.4421745996289254E+00,
     0.9927049442755639E+00,
     0.4202103819530612E-01,
     0.9796589705830716E+00,
     0.9226039842296429E+00,
     0.4470785799755852E+00,
     0.7444549220718699E+00 };

  static double x_vec[N_MAX] = {
     0.30E-01,
     0.30E+00,
     0.15E+01,
     0.75E-01,
     0.75E+00,
     0.35E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.15E+00,
     0.15E+01,
     0.70E+01,
     0.25E+01,
     0.12E+02,
     0.16E+02,
     0.25E+02,
     0.45E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_inc_q_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_Q_VALUES: values of the normalized incomplete Gamma function Q(A,X).
//
//  Discussion:
//
//    The (normalized) incomplete Gamma function is defined as:
//
//      Q(A,X) = 1/Gamma(A) * Integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
//
//    With this definition, for all A and X,
//
//      0 <= Q(A,X) <= 1
//
//    and
//
//      Q(A,oo) = 0.0
//
//    In Mathematica, the function can be evaluated by:
//
//      GammaRegularized[A,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.10E+01,
     0.10E+01,
     0.10E+01,
     0.11E+01,
     0.11E+01,
     0.11E+01,
     0.20E+01,
     0.20E+01,
     0.20E+01,
     0.60E+01,
     0.60E+01,
     0.11E+02,
     0.26E+02,
     0.41E+02  };

  static double fx_vec[N_MAX] = {
    0.2617649467660649E+00,
    0.09164201026996572E+00,
    0.01134401663780527E+00,
    0.6985353583033387E+00,
    0.2206713619198468E+00,
    0.008150971593502700E+00,
    0.9048374180359596E+00,
    0.3678794411714423E+00,
    0.006737946999085467E+00,
    0.9279402542394568E+00,
    0.4108190381293515E+00,
    0.008463184015447498E+00,
    0.9898141728888165E+00,
    0.5578254003710746E+00,
    0.007295055724436130E+00,
    0.9579789618046939E+00,
    0.02034102941692837E+00,
    0.07739601577035708E+00,
    0.5529214200244148E+00,
    0.2555450779281301E+00  };

  static double x_vec[N_MAX] = {
     0.30E-01,
     0.30E+00,
     0.15E+01,
     0.75E-01,
     0.75E+00,
     0.35E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.15E+00,
     0.15E+01,
     0.70E+01,
     0.25E+01,
     0.12E+02,
     0.16E+02,
     0.25E+02,
     0.45E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_inc_tricomi_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_TRICOMI_VALUES: values of Tricomi's incomplete Gamma function.
//
//  Discussion:
//
//    Tricomi's incomplete Gamma function is defined as:
//
//      1/Gamma(A) * 1/X^A * Integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.10E+01,
     0.10E+01,
     0.10E+01,
     0.11E+01,
     0.11E+01,
     0.11E+01,
     0.20E+01,
     0.20E+01,
     0.20E+01,
     0.60E+01,
     0.60E+01,
     0.11E+02,
     0.26E+02,
     0.41E+02  };

  static double fx_vec[N_MAX] = {
    1.048292641463504E+00,
    1.024577737369574E+00,
    0.9493712443185374E+00,
    1.100793230316492E+00,
    0.8998911979655218E+00,
    0.5301656062431039E+00,
    0.9516258196404043E+00,
    0.6321205588285577E+00,
    0.1986524106001829E+00,
    0.9071784510537487E+00,
    0.5891809618706485E+00,
    0.1688269752193589E+00,
    0.4527034271637121E+00,
    0.1965220442795224E+00,
    0.02025928457705232E+00,
    0.0001721181724479739E+00,
    3.280858070850586E-07,
    5.244396471821590E-14,
    2.013462926183376E-37,
    1.230623887499875E-68 };

  static double x_vec[N_MAX] = {
     0.30E-01,
     0.30E+00,
     0.15E+01,
     0.75E-01,
     0.75E+00,
     0.35E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.15E+00,
     0.15E+01,
     0.70E+01,
     0.25E+01,
     0.12E+02,
     0.16E+02,
     0.25E+02,
     0.45E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_inc_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
//
//  Discussion:
//
//    The incomplete Gamma function is defined as:
//
//      Integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
//
//    In Mathematica, the function can be evaluated by:
//
//      Gamma[A,X]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.10E+01,
     0.10E+01,
     0.10E+01,
     0.11E+01,
     0.11E+01,
     0.11E+01,
     0.20E+01,
     0.20E+01,
     0.20E+01,
     0.60E+01,
     0.60E+01,
     0.11E+02,
     0.26E+02,
     0.41E+02  };

  static double fx_vec[N_MAX] = {
    2.490302836300570E+00,
    0.8718369702247978E+00,
    0.1079213896175866E+00,
    1.238121685818417E+00,
    0.3911298052193973E+00,
    0.01444722098952533E+00,
    0.9048374180359596E+00,
    0.3678794411714423E+00,
    0.006737946999085467E+00,
    0.8827966752611692E+00,
    0.3908330082003269E+00,
    0.008051456628620993E+00,
    0.9898141728888165E+00,
    0.5578254003710746E+00,
    0.007295055724436130E+00,
    114.9574754165633E+00,
    2.440923530031405E+00,
    280854.6620274718E+00,
    8.576480283455533E+24,
    2.085031346403364E+47 };

  static double x_vec[N_MAX] = {
     0.30E-01,
     0.30E+00,
     0.15E+01,
     0.75E-01,
     0.75E+00,
     0.35E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.10E+00,
     0.10E+01,
     0.50E+01,
     0.15E+00,
     0.15E+01,
     0.70E+01,
     0.25E+01,
     0.12E+02,
     0.16E+02,
     0.25E+02,
     0.45E+02 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gamma_log_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GAMMA_LOG_VALUES returns some values of the Log Gamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[Gamma[x]]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.1524063822430784E+01,
      0.7966778177017837E+00,
      0.3982338580692348E+00,
      0.1520596783998375E+00,
      0.0000000000000000E+00,
     -0.4987244125983972E-01,
     -0.8537409000331584E-01,
     -0.1081748095078604E+00,
     -0.1196129141723712E+00,
     -0.1207822376352452E+00,
     -0.1125917656967557E+00,
     -0.9580769740706586E-01,
     -0.7108387291437216E-01,
     -0.3898427592308333E-01,
     0.00000000000000000E+00,
     0.69314718055994530E+00,
     0.17917594692280550E+01,
     0.12801827480081469E+02,
     0.39339884187199494E+02,
     0.71257038967168009E+02 };

  static double x_vec[N_MAX] = {
      0.20E+00,
      0.40E+00,
      0.60E+00,
      0.80E+00,
      1.00E+00,
      1.10E+00,
      1.20E+00,
      1.30E+00,
      1.40E+00,
      1.50E+00,
      1.60E+00,
      1.70E+00,
      1.80E+00,
      1.90E+00,
      2.00E+00,
      3.00E+00,
      4.00E+00,
     10.00E+00,
     20.00E+00,
     30.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gegenbauer_poly_values ( int &n_data, int &n, double &a, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GEGENBAUER_POLY_VALUES returns some values of the Gegenbauer polynomials.
//
//  Discussion:
//
//    The Gegenbauer polynomials are also known as the "spherical
//    polynomials" or "ultraspherical polynomials".
//
//    In Mathematica, the function can be evaluated by:
//
//      GegenbauerC[n,m,x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order parameter of the function.
//
//    Output, double &A, the real parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 38

  static double a_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.0E+00,
      1.0E+00,
      2.0E+00,
      3.0E+00,
      4.0E+00,
      5.0E+00,
      6.0E+00,
      7.0E+00,
      8.0E+00,
      9.0E+00,
     10.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00,
      3.0E+00 };

  static double fx_vec[N_MAX] = {
       1.0000000000E+00,
       0.2000000000E+00,
      -0.4400000000E+00,
      -0.2800000000E+00,
       0.2320000000E+00,
       0.3075200000E+00,
      -0.0805760000E+00,
      -0.2935168000E+00,
      -0.0395648000E+00,
       0.2459712000E+00,
       0.1290720256E+00,
       0.0000000000E+00,
      -0.3600000000E+00,
      -0.0800000000E+00,
       0.8400000000E+00,
       2.4000000000E+00,
       4.6000000000E+00,
       7.4400000000E+00,
      10.9200000000E+00,
      15.0400000000E+00,
      19.8000000000E+00,
      25.2000000000E+00,
      -9.0000000000E+00,
      -0.1612800000E+00,
      -6.6729600000E+00,
      -8.3750400000E+00,
      -5.5267200000E+00,
       0.0000000000E+00,
       5.5267200000E+00,
       8.3750400000E+00,
       6.6729600000E+00,
       0.1612800000E+00,
      -9.0000000000E+00,
     -15.4252800000E+00,
      -9.6969600000E+00,
      22.4409600000E+00,
     100.8892800000E+00,
     252.0000000000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  2,
     2,  2,  2,
     2,  2,  2,
     2,  2,  2,
     2,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5,  5,
     5,  5 };

  static double x_vec[N_MAX] = {
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.20E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
      0.40E+00,
     -0.50E+00,
     -0.40E+00,
     -0.30E+00,
     -0.20E+00,
     -0.10E+00,
      0.00E+00,
      0.10E+00,
      0.20E+00,
      0.30E+00,
      0.40E+00,
      0.50E+00,
      0.60E+00,
      0.70E+00,
      0.80E+00,
      0.90E+00,
      1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void geometric_cdf_values ( int &n_data, int &x, double &p, double &cdf )

//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRIC_CDF_VALUES returns values of the geometric CDF.
//
//  Discussion:
//
//    The geometric or Pascal probability density function gives the
//    probability that the first success will happen on the X-th Bernoulli
//    trial, given that the probability of a success on a single trial is P.
//
//    The value of CDF ( X, P ) is the probability that the first success
//    will happen on or before the X-th trial.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = GeometricDistribution [ p ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger and Stephen Kokoska,
//    CRC Standard Probability and Statistics Tables and Formulae,
//    Chapman and Hall / CRC Press, 2000.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &X, the number of trials.
//
//    Output, double &P, the probability of success
//    on one trial.
//
//    Output, double &CDF, the cumulative density function.
//
{
# define N_MAX 14

  static double cdf_vec[N_MAX] = {
     0.1900000000000000E+00,
     0.2710000000000000E+00,
     0.3439000000000000E+00,
     0.6861894039100000E+00,
     0.3600000000000000E+00,
     0.4880000000000000E+00,
     0.5904000000000000E+00,
     0.9141006540800000E+00,
     0.7599000000000000E+00,
     0.8704000000000000E+00,
     0.9375000000000000E+00,
     0.9843750000000000E+00,
     0.9995117187500000E+00,
     0.9999000000000000E+00 };

  static double p_vec[N_MAX] = {
     0.1E+00,
     0.1E+00,
     0.1E+00,
     0.1E+00,
     0.2E+00,
     0.2E+00,
     0.2E+00,
     0.2E+00,
     0.3E+00,
     0.4E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.9E+00 };

  static int x_vec[N_MAX] = {
    1,  2,  3, 10, 1,
    2,  3, 10,  3, 3,
    3,  5, 10,  3 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0;
    p = 0.0;
    cdf = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    p = p_vec[n_data-1];
    cdf = cdf_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void goodwin_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GOODWIN_VALUES returns some values of the Goodwin and Staton function.
//
//  Discussion:
//
//    The function is defined by:
//
//      GOODWIN(x) = Integral ( 0 <= t < +oo ) exp ( -t^2 ) / ( t + x ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.59531540040441651584E+01,
      0.45769601268624494109E+01,
      0.32288921331902217638E+01,
      0.19746110873568719362E+01,
      0.96356046208697728563E+00,
      0.60513365250334458174E+00,
      0.51305506459532198016E+00,
      0.44598602820946133091E+00,
      0.37344458206879749357E+00,
      0.35433592884953063055E+00,
      0.33712156518881920994E+00,
      0.29436170729362979176E+00,
      0.25193499644897222840E+00,
      0.22028778222123939276E+00,
      0.19575258237698917033E+00,
      0.17616303166670699424E+00,
      0.16015469479664778673E+00,
      0.14096116876193391066E+00,
      0.13554987191049066274E+00,
      0.11751605060085098084E+00 };

  static double x_vec[N_MAX] = {
      0.0019531250E+00,
      0.0078125000E+00,
      0.0312500000E+00,
      0.1250000000E+00,
      0.5000000000E+00,
      1.0000000000E+00,
      1.2500000000E+00,
      1.5000000000E+00,
      1.8750000000E+00,
      2.0000000000E+00,
      2.1250000000E+00,
      2.5000000000E+00,
      3.0000000000E+00,
      3.5000000000E+00,
      4.0000000000E+00,
      4.5000000000E+00,
      5.0000000000E+00,
      5.7500000000E+00,
      6.0000000000E+00,
      7.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void gud_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    GUD_VALUES returns some values of the Gudermannian function.
//
//  Discussion:
//
//    The Gudermannian function relates the hyperbolic and trigonomentric
//    functions.  For any argument X, there is a corresponding value
//    GD so that
//
//      SINH(X) = TAN(GD).
//
//    This value GD is called the Gudermannian of X and symbolized
//    GD(X).  The inverse Gudermannian function is given as input a value
//    GD and computes the corresponding value X.
//
//    GD(X) = 2 * arctan ( exp ( X ) ) - PI / 2
//
//    In Mathematica, the function can be evaluated by:
//
//      2 * Atan[Exp[x]] - Pi/2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     -0.1301760336046015E+01,
     -0.8657694832396586E+00,
      0.0000000000000000E+00,
      0.9983374879348662E-01,
      0.1986798470079397E+00,
      0.4803810791337294E+00,
      0.8657694832396586E+00,
      0.1131728345250509E+01,
      0.1301760336046015E+01,
      0.1406993568936154E+01,
      0.1471304341117193E+01,
      0.1510419907545700E+01,
      0.1534169144334733E+01 };

  static double x_vec[N_MAX] = {
     -2.0E+00,
     -1.0E+00,
      0.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      1.5E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hermite_function_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_FUNCTION_VALUES returns some values of the Hermite function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Hf(n,x) = HermiteH[n,x] 
//        * Exp [ -1/2 * x^2] / Sqrt [ 2^n * n! * Sqrt[Pi] ]
//
//    The Hermite functions are orthonormal:
//
//      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 23

  static double fx_vec[N_MAX] = {
    0.7511255444649425E+00,  0.0000000000000000E+00, -0.5311259660135985E+00, 
    0.0000000000000000E+00,  0.4599685791773266E+00,  0.0000000000000000E+00, 
    0.4555806720113325E+00,  0.6442883651134752E+00,  0.3221441825567376E+00, 
   -0.2630296236233334E+00, -0.4649750762925110E+00, -0.5881521185179581E-01, 
    0.3905052515434106E+00,  0.2631861423064045E+00, -0.2336911435996523E+00, 
   -0.3582973361472840E+00,  0.6146344487883041E-01,  0.3678312067984882E+00, 
    0.9131969309166278E-01,  0.4385750950032321E+00, -0.2624689527931006E-01, 
    0.5138426125477819E+00,  0.9355563118061758E-01 };

  static int n_vec[N_MAX] = {
    0,  1,  2,  
    3,  4,  5,  
    0,  1,  2,  
    3,  4,  5,  
    6,  7,  8,  
    9, 10, 11,  
   12,  5,  5,  
    5,  5 };

  static double x_vec[N_MAX] = {
    0.0E+00, 0.0E+00, 0.0E+00, 
    0.0E+00, 0.0E+00, 0.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 1.0E+00, 1.0E+00, 
    1.0E+00, 0.5E+00, 2.0E+00, 
    3.0E+00, 4.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hermite_poly_phys_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PHYS_VALUES returns some values of the physicist's Hermite polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      HermiteH[n,x]
//
//  Differential equation:
//
//    Y'' - 2 X Y' + 2 N Y = 0;
//
//  First terms:
//
//      1
//      2 X
//      4 X^2     -  2
//      8 X^3     - 12 X
//     16 X^4     - 48 X^2     + 12
//     32 X^5    - 160 X^3    + 120 X
//     64 X^6    - 480 X^4    + 720 X^2    - 120
//    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
//    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
//    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
//   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
//
//  Recursion:
//
//    H(0,X) = 1,
//    H(1,X) = 2*X,
//    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
//    = sqrt ( PI ) * 2^N * N!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.1000000000000000E+02,
      0.9800000000000000E+02,
      0.9400000000000000E+03,
      0.8812000000000000E+04,
      0.8060000000000000E+05,
      0.7178800000000000E+06,
      0.6211600000000000E+07,
      0.5206568000000000E+08,
      0.4212712000000000E+09,
      0.3275529760000000E+10,
      0.2432987360000000E+11,
      0.1712370812800000E+12,
      0.0000000000000000E+00,
      0.4100000000000000E+02,
     -0.8000000000000000E+01,
      0.3816000000000000E+04,
      0.3041200000000000E+07 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5,  5 };

  static double x_vec[N_MAX] = {
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     0.0E+00,
     0.5E+00,
     1.0E+00,
     3.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hermite_poly_prob_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HERMITE_POLY_PROB_VALUES: values of the probabilist's Hermite polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      He(n,x) = HermiteH[n,x/Sqrt[2]] / Sqrt [ 2^n ] 
//
//  First terms:
//
//   1
//   X
//   X^2  -  1
//   X^3  -  3 X
//   X^4  -  6 X^2 +   3
//   X^5  - 10 X^3 +  15 X
//   X^6  - 15 X^4 +  45 X^2 -   15
//   X^7  - 21 X^5 + 105 X^3 -  105 X
//   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
//   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
//   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
//
//  Recursion:
//
//    He(0,X) = 1,
//    He(1,X) = X,
//    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
//
//  Norm:
//
//    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
//    = sqrt ( 2 * pi ) * N! * delta ( M, N )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
    1.000000000000000E+00, 
    5.000000000000000E+00, 
    24.00000000000000E+00, 
    110.0000000000000E+00,
    478.0000000000000E+00,
    1950.000000000000E+00, 
    7360.000000000000E+00, 
    25100.00000000000E+00, 
    73980.00000000000E+00, 
    169100.0000000000E+00, 
    179680.0000000000E+00, 
   -792600.0000000000E+00, 
   -5939480.000000000E+00, 
    0.000000000000000E+00, 
    6.281250000000000E+00, 
    6.000000000000000E+00, 
    18.00000000000000E+00, 
    90150.00000000000E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5,  5 };

  static double x_vec[N_MAX] = {
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     5.0E+00,
     0.0E+00,
     0.5E+00,
     1.0E+00,
     3.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hyper_1f1_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPER_1F1_VALUES returns some values of the hypergeometric function 1F1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = Hypergeometric1F1 [ a, b, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, &X, the parameters of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 24

  static double a_vec[N_MAX] = {
    -2.500,
    -0.500,
     0.500,
     2.500,
    -2.500,
    -0.500,
     0.500,
     2.500,
    -2.500,
    -0.500,
     0.500,
     2.500,
     0.825,
     1.100,
     1.650,
     3.300,
     0.825,
     1.100,
     1.650,
     3.300,
     0.825,
     1.100,
     1.650,
     3.300 };
  static double b_vec[N_MAX] = {
     3.3,
     1.1,
     1.1,
     3.3,
     3.3,
     1.1,
     1.1,
     3.3,
     3.3,
     1.1,
     1.1,
     3.3,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7 };
  static double fx_vec[N_MAX] = {
     0.81879926689265186854,
     0.88283984828032972070,
     1.1245023764952626690,
     1.2101049301639599598,
     0.12723045536781567174,
     0.12326016871544045107,
     2.3297954665128293051,
     3.3890020264468009733,
    -0.18819510282516768874,
    -1.0764203806547022727,
     5.7521824680907968433,
     9.9998567403304086593,
     1.0317208964319891384,
     1.0424867029249952040,
     1.0643112000949092012,
     1.1321844369742336326,
     1.2328402688568452181,
     1.3200654482027340732,
     1.5104811522310825217,
     2.2307520785940524365,
     1.5197286298183137741,
     1.7364938170250847619,
     2.2492330307668135926,
     4.6377737119178965298 };
  static double x_vec[N_MAX] = {
     0.25,
     0.25,
     0.25,
     0.25,
     1.55,
     1.55,
     1.55,
     1.55,
     2.85,
     2.85,
     2.85,
     2.85,
     0.25,
     0.25,
     0.25,
     0.25,
     1.55,
     1.55,
     1.55,
     1.55,
     2.85,
     2.85,
     2.85,
     2.85 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hyper_2f1_values ( int &n_data, double &a, double &b, double &c,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPER_2F1_VALUES returns some values of the hypergeometric function 2F1.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = Hypergeometric2F1 [ a, b, c, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 September 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Shanjie Zhang, Jianming Jin,
//    Computation of Special Functions,
//    Wiley, 1996,
//    ISBN: 0-471-11963-6,
//    LC: QA351.C45
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, &C, &X, the parameters of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 24

  static double a_vec[N_MAX] = {
   -2.5,
   -0.5,
    0.5,
    2.5,
   -2.5,
   -0.5,
    0.5,
    2.5,
   -2.5,
   -0.5,
    0.5,
    2.5,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3 };
  static double b_vec[N_MAX] = {
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    3.3,
    1.1,
    1.1,
    3.3,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7 };
  static double c_vec[N_MAX] = {
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
    6.7,
   -5.5,
   -0.5,
    0.5,
    4.5,
   -5.5,
   -0.5,
    0.5,
    4.5,
   -5.5,
   -0.5,
    0.5,
    4.5 };
  static double fx_vec[N_MAX] = {
    0.72356129348997784913,
    0.97911109345277961340,
    1.0216578140088564160,
    1.4051563200112126405,
    0.46961431639821611095,
    0.95296194977446325454,
    1.0512814213947987916,
    2.3999062904777858999,
    0.29106095928414718320,
    0.92536967910373175753,
    1.0865504094806997287,
    5.7381565526189046578,
    15090.669748704606754,
   -104.31170067364349677,
    21.175050707768812938,
    4.1946915819031922850,
    1.0170777974048815592E+10,
   -24708.635322489155868,
    1372.2304548384989560,
    58.092728706394652211,
    5.8682087615124176162E+18,
   -4.4635010147295996680E+08,
    5.3835057561295731310E+06,
    20396.913776019659426 };
  static double x_vec[N_MAX] = {
    0.25,
    0.25,
    0.25,
    0.25,
    0.55,
    0.55,
    0.55,
    0.55,
    0.85,
    0.85,
    0.85,
    0.85,
    0.25,
    0.25,
    0.25,
    0.25,
    0.55,
    0.55,
    0.55,
    0.55,
    0.85,
    0.85,
    0.85,
    0.85 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    c = c_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hypergeometric_cdf_values ( int &n_data, int &sam, int &suc, int &pop,
  int &n, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_CDF_VALUES returns some values of the hypergeometric CDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of at most X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      CDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &SAM, int &SUC, int &POP, the sample size,
//    success size, and population parameters of the function.
//
//    Output, int &N, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.6001858177500578E-01,
     0.2615284665839845E+00,
     0.6695237889132748E+00,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.5332595856827856E+00,
     0.1819495964117640E+00,
     0.4448047017527730E-01,
     0.9999991751316731E+00,
     0.9926860896560750E+00,
     0.8410799901444538E+00,
     0.3459800113391901E+00,
     0.0000000000000000E+00,
     0.2088888139634505E-02,
     0.3876752992448843E+00,
     0.9135215248834896E+00 };

  static int n_vec[N_MAX] = {
     7,  8,  9, 10,
     6,  6,  6,  6,
     6,  6,  6,  6,
     0,  0,  0,  0 };

  static int pop_vec[N_MAX] = {
    100, 100, 100, 100,
    100, 100, 100, 100,
    100, 100, 100, 100,
    90,  200, 1000, 10000 };

  static int sam_vec[N_MAX] = {
    10, 10, 10, 10,
     6,  7,  8,  9,
    10, 10, 10, 10,
    10, 10, 10, 10 };

  static int suc_vec[N_MAX] = {
    90, 90, 90, 90,
    90, 90, 90, 90,
    10, 30, 50, 70,
    90, 90, 90, 90 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    sam = 0;
    suc = 0;
    pop = 0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    sam = sam_vec[n_data-1];
    suc = suc_vec[n_data-1];
    pop = pop_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hypergeometric_pdf_values ( int &n_data, int &sam, int &suc, int &pop,
  int &n, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_PDF_VALUES returns some values of the hypergeometric PDF.
//
//  Discussion:
//
//    CDF(X)(A,B) is the probability of X successes in A trials,
//    given that the probability of success on a single trial is B.
//
//    In Mathematica, the function can be evaluated by:
//
//      dist = HypergeometricDistribution [ sam, suc, pop ]
//      PDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 January 2008
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &SAM, int &SUC, int &POP, the sample size,
//    success size, and population parameters of the function.
//
//    Output, int &N, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
    0.05179370533242827E+00,
    0.2015098848089788E+00,
    0.4079953223292903E+00,
    0.3304762110867252E+00,
    0.5223047493549780E+00,
    0.3889503452643453E+00,
    0.1505614239732950E+00,
    0.03927689321042477E+00,
    0.00003099828465518108E+00,
    0.03145116093938197E+00,
    0.2114132170316862E+00,
    0.2075776621999210E+00,
    0.0000000000000000E+00,
    0.002088888139634505E+00,
    0.3876752992448843E+00,
    0.9135215248834896E+00 };

  static int n_vec[N_MAX] = {
     7,  8,  9, 10,
     6,  6,  6,  6,
     6,  6,  6,  6,
     0,  0,  0,  0 };

  static int pop_vec[N_MAX] = {
    100, 100, 100, 100,
    100, 100, 100, 100,
    100, 100, 100, 100,
    90,  200, 1000, 10000 };

  static int sam_vec[N_MAX] = {
    10, 10, 10, 10,
     6,  7,  8,  9,
    10, 10, 10, 10,
    10, 10, 10, 10 };

  static int suc_vec[N_MAX] = {
    90, 90, 90, 90,
    90, 90, 90, 90,
    10, 30, 50, 70,
    90, 90, 90, 90 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    sam = 0;
    suc = 0;
    pop = 0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    sam = sam_vec[n_data-1];
    suc = suc_vec[n_data-1];
    pop = pop_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void hypergeometric_u_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    HYPERGEOMETRIC_U_VALUES: some values of the hypergeometric function U(a,b,x).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = HypergeometricU [ a, b, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 651-652.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, &X, the parameters of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 24

  static double a_vec[N_MAX] = {
    -2.500,
    -0.500,
     0.500,
     2.500,
    -2.500,
    -0.500,
     0.500,
     2.500,
    -2.500,
    -0.500,
     0.500,
     2.500,
     0.825,
     1.100,
     1.650,
     3.300,
     0.825,
     1.100,
     1.650,
     3.300,
     0.825,
     1.100,
     1.650,
     3.300 };
  static double b_vec[N_MAX] = {
     3.3,
     1.1,
     1.1,
     3.3,
     3.3,
     1.1,
     1.1,
     3.3,
     3.3,
     1.1,
     1.1,
     3.3,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7,
     6.7 };
  static double fx_vec[N_MAX] = {
         -68.693628728078601389,
          -0.0029710551374761070801,
           1.5008631742177797301,
          20.614688244200596134,
           7.4563815469305551938,
           1.0155793767749293733,
           0.73446538936622668912,
           0.28046404941879399225,
           3.4508153741446547607,
           1.5156637368753063495,
           0.56042118587934993510,
           0.064897147735134223341,
      223432.02356977463356,
      263079.25980740811495,
      269802.90319351274132,
       82809.311335606553425,
          26.465684783131844524,
          28.093506172516056560,
          23.889164624518872504,
           4.5338847857070388229,
           3.0224469362694842535,
           2.8040650913713359934,
           1.9262578111480172682,
           0.23020518115860909098 };
  static double x_vec[N_MAX] = {
     0.25,
     0.25,
     0.25,
     0.25,
     1.55,
     1.55,
     1.55,
     1.55,
     2.85,
     2.85,
     2.85,
     2.85,
     0.25,
     0.25,
     0.25,
     0.25,
     1.55,
     1.55,
     1.55,
     1.55,
     2.85,
     2.85,
     2.85,
     2.85 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void i0ml0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    I0ML0_VALUES returns some values of the I0ML0 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      I0ML0(x) = I0(x) - L0(x)
//
//    I0(x) is the modified Bessel function of the first kind of order 0,
//    L0(x) is the modified Struve function of order 0.
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.99875755515461749793E+00,
     0.99011358230706643807E+00,
     0.92419435310023947018E+00,
     0.73624267134714273902E+00,
     0.55582269181411744686E+00,
     0.34215154434462160628E+00,
     0.17087174888774706539E+00,
     0.81081008709219208918E-01,
     0.53449421441089580702E-01,
     0.39950321008923244846E-01,
     0.39330637437584921392E-01,
     0.37582274342808670750E-01,
     0.31912486554480390343E-01,
     0.25506146883504738403E-01,
     0.21244480317825292412E-01,
     0.15925498348551684335E-01,
     0.12737506927242585015E-01,
     0.84897750814784916847E-02,
     0.63668349178454469153E-02,
     0.50932843163122551114E-02 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0156250000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
       4.0000000000E+00,
       8.0000000000E+00,
      12.0000000000E+00,
      16.0000000000E+00,
      16.2500000000E+00,
      17.0000000000E+00,
      20.0000000000E+00,
      25.0000000000E+00,
      30.0000000000E+00,
      40.0000000000E+00,
      50.0000000000E+00,
      75.0000000000E+00,
     100.0000000000E+00,
     125.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void i1ml1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    I1ML1_VALUES returns some values of the I1ML1 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      I1ML1(x) = I1(x) - L1(x)
//
//    I1(x) is the modified Bessel function of the first kind of order 1,
//    L1(x) is the modified Struve function of order 1.
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.97575346155386267134E-03,
     0.77609293280609272733E-02,
     0.59302966404545373770E-01,
     0.20395212276737365307E+00,
     0.33839472293667639038E+00,
     0.48787706726961324579E+00,
     0.59018734196576517506E+00,
     0.62604539530312149476E+00,
     0.63209315274909764698E+00,
     0.63410179313235359215E+00,
     0.63417966797578128188E+00,
     0.63439268632392089434E+00,
     0.63501579073257770690E+00,
     0.63559616677359459337E+00,
     0.63591001826697110312E+00,
     0.63622113181751073643E+00,
     0.63636481702133606597E+00,
     0.63650653499619902120E+00,
     0.63655609126300261851E+00,
     0.63657902087183929223E+00 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0156250000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
       4.0000000000E+00,
       8.0000000000E+00,
      12.0000000000E+00,
      16.0000000000E+00,
      16.2500000000E+00,
      17.0000000000E+00,
      20.0000000000E+00,
      25.0000000000E+00,
      30.0000000000E+00,
      40.0000000000E+00,
      50.0000000000E+00,
      75.0000000000E+00,
     100.0000000000E+00,
     125.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void int_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    INT_VALUES returns some values of the "integer part" function.
//
//  Discussion:
//
//    INT(X) = returns the integer part of a real number.
//
//    The result is returned as a real number.
//
//    The result is computed by rounding the absolute value of the
//   input towards 0, and then restoring the sign.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output double FX, the value of the function.
//
{
# define N_MAX 25

  static double fx_vec[N_MAX] = {
     -2.00E+00, 
     -1.00E+00, 
     -1.00E+00, 
     -1.00E+00,     
     -1.00E+00,     
     -1.00E+00,       
      0.00E+00, 
      0.00E+00,     
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      0.00E+00, 
      1.00E+00, 
      1.00E+00, 
      1.00E+00, 
      1.00E+00, 
      1.00E+00, 
      2.00E+00 };

  static double x_vec[N_MAX] = {
     -2.01E+00, 
     -1.99E+00, 
     -1.50E+00, 
     -1.10E+00,     
     -1.01E+00,     
     -1.00E+00,       
     -0.99E+00, 
     -0.90E+00,     
     -0.51E+00, 
     -0.50E+00, 
     -0.49E+00, 
     -0.01E+00, 
      0.00E+00, 
      0.01E+00, 
      0.49E+00, 
      0.50E+00, 
      0.51E+00, 
      0.90E+00, 
      0.99E+00, 
      1.00E+00, 
      1.01E+00, 
      1.10E+00, 
      1.50E+00, 
      1.99E+00, 
      2.01E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jacobi_cn_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_CN_VALUES returns some values of the Jacobi elliptic function CN(A,X).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiCN[ x, a ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00 };

  static double fx_vec[N_MAX] = {
      0.9950041652780258E+00,
      0.9800665778412416E+00,
      0.8775825618903727E+00,
      0.5403023058681397E+00,
     -0.4161468365471424E+00,
      0.9950124626090582E+00,
      0.9801976276784098E+00,
      0.8822663948904403E+00,
      0.5959765676721407E+00,
     -0.1031836155277618E+00,
      0.9950207489532265E+00,
      0.9803279976447253E+00,
      0.8868188839700739E+00,
      0.6480542736638854E+00,
      0.2658022288340797E+00,
      0.3661899347368653E-01,
      0.9803279976447253E+00,
      0.8868188839700739E+00,
      0.6480542736638854E+00,
      0.2658022288340797E+00 };

  static double x_vec[N_MAX] = {
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      4.0E+00,
     -0.2E+00,
     -0.5E+00,
     -1.0E+00,
     -2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jacobi_dn_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_DN_VALUES returns some values of the Jacobi elliptic function DN(A,X).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiDN[ x, a ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00 };

  static double fx_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.9975093485144243E+00,
     0.9901483195224800E+00,
     0.9429724257773857E+00,
     0.8231610016315963E+00,
     0.7108610477840873E+00,
     0.9950207489532265E+00,
     0.9803279976447253E+00,
     0.8868188839700739E+00,
     0.6480542736638854E+00,
     0.2658022288340797E+00,
     0.3661899347368653E-01,
     0.9803279976447253E+00,
     0.8868188839700739E+00,
     0.6480542736638854E+00,
     0.2658022288340797E+00  };

  static double x_vec[N_MAX] = {
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      4.0E+00,
     -0.2E+00,
     -0.5E+00,
     -1.0E+00,
     -2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jacobi_poly_values ( int &n_data, int &n, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_POLY_VALUES returns some values of the Jacobi polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiP[ n, a, b, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the degree of the polynomial.
//
//    Output, double &A, &B, parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 26

  static int a_vec[N_MAX] = {
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 1.0, 2.0,
     3.0, 4.0, 5.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0, 0.0, 0.0,
     0.0, 0.0 };

  static int b_vec[N_MAX] = {
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 2.0,
    3.0, 4.0, 5.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0, 1.0, 1.0,
    1.0, 1.0 };

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.3750000000000000E+00,
     -0.4843750000000000E+00,
     -0.1328125000000000E+00,
      0.2753906250000000E+00,
     -0.1640625000000000E+00,
     -0.1174804687500000E+01,
     -0.2361328125000000E+01,
     -0.2616210937500000E+01,
      0.1171875000000000E+00,
      0.4218750000000000E+00,
      0.5048828125000000E+00,
      0.5097656250000000E+00,
      0.4306640625000000E+00,
     -0.6000000000000000E+01,
      0.3862000000000000E-01,
      0.8118400000000000E+00,
      0.3666000000000000E-01,
     -0.4851200000000000E+00,
     -0.3125000000000000E+00,
      0.1891200000000000E+00,
      0.4023400000000000E+00,
      0.1216000000000000E-01,
     -0.4396200000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0, 1, 2, 3,
     4, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5, 5, 5,
     5, 5 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
      0.5E+00,
     -1.0E+00,
     -0.8E+00,
     -0.6E+00,
     -0.4E+00,
     -0.2E+00,
      0.0E+00,
      0.2E+00,
      0.4E+00,
      0.6E+00,
      0.8E+00,
      1.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jacobi_sn_values ( int &n_data, double &a, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    JACOBI_SN_VALUES returns some values of the Jacobi elliptic function SN(A,X).
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      JacobiSN[ x, a ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 May 2001
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00 };

  static double fx_vec[N_MAX] = {
      0.9983341664682815E-01,
      0.1986693307950612E+00,
      0.4794255386042030E+00,
      0.8414709848078965E+00,
      0.9092974268256817E+00,
      0.9975068547462484E-01,
      0.1980217429819704E+00,
      0.4707504736556573E+00,
      0.8030018248956439E+00,
      0.9946623253580177E+00,
      0.9966799462495582E-01,
      0.1973753202249040E+00,
      0.4621171572600098E+00,
      0.7615941559557649E+00,
      0.9640275800758169E+00,
      0.9993292997390670E+00,
     -0.1973753202249040E+00,
     -0.4621171572600098E+00,
     -0.7615941559557649E+00,
     -0.9640275800758169E+00  };

  static double x_vec[N_MAX] = {
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      0.1E+00,
      0.2E+00,
      0.5E+00,
      1.0E+00,
      2.0E+00,
      4.0E+00,
     -0.2E+00,
     -0.5E+00,
     -1.0E+00,
     -2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jed_ce_values ( int &n_data, double &jed, int &y, int &m, int &d,
  double &f )

//****************************************************************************80
//
//  Purpose:
//
//    JED_CE_VALUES returns the Common Era dates for Julian Ephemeris Dates.
//
//  Discussion:
//
//    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
//    starting from noon on 1 January 4713 BCE.
//
//    The CE or Common Era is the day, month and year under the
//    hybrid Julian/Gregorian Calendar, with a transition from Julian
//    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
//
//    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
//    years BC/BCE are indicated by a negative year value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 May 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold and Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &JED, the Julian Ephemeris Date.
//
//    Output, int &Y, &M, &D, the Common Era date.
//
//    Output, double &F, the fractional part of the day.
//
{
# define N_MAX 51

  static int d_vec[N_MAX] = {
    1,
    2,
    26,
    8,
    6,
    18,
    8,
    9,
    1,
    26,
    26,
    1,
    1,
    29,
    31,
    1,
    3,
    3,
    29,
    24,
    24,
    29,
    3,
    11,
    12,
    24,
    19,
    15,
    16,
    16,
    21,
    17,
    9,
    4,
    15,
    4,
    13,
    14,
    18,
    22,
    21,
    24,
    17,
    31,
    1,
    6,
    25,
    1,
    9,
    23,
    1 };
  static double f_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.00,
    0.00,
    0.25,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.50,
    0.50,
    0.00,
    0.50,
    0.50,
    0.00,
    0.00,
    0.00,
    0.00,
    0.00,
    0.81,
    0.00,
    0.00,
    0.00,
    0.00,
    0.33,
    0.00,
    0.50 };
  static double jed_vec[N_MAX] = {
           0.00,
           1.00,
      259261.00,
      347998.50,
      584282.50,
      588465.75,
      758325.50,
     1438178.50,
     1446389.50,
     1448637.50,
     1448637.50,
     1607708.50,
     1607738.50,
     1713262.50,
     1721422.50,
     1721423.50,
     1721425.50,
     1721425.50,
     1724220.50,
     1741959.50,
     1749994.50,
     1825029.50,
     1862836.50,
     1922867.50,
     1936747.50,
     1940351.50,
     1948320.50,
     1948438.50,
     1948439.50,
     1952062.50,
     1952067.50,
     2114872.50,
     2289425.50,
     2299160.00,
     2299161.00,
     2333269.50,
     2361221.00,
     2361222.00,
     2372547.50,
     2375839.50,
     2394646.50,
     2394710.50,
     2400000.50,
     2415020.31,
     2440587.50,
     2444244.50,
     2450138.50,
     2451544.50,
     2453073.83,
     2456284.50,
     2913943.00 };
  static int m_vec[N_MAX] = {
     1,
     1,
     10,
     10,
     9,
     2,
     3,
     7,
     1,
     2,
     2,
     9,
     10,
     8,
     12,
     1,
     1,
     1,
     8,
     3,
     3,
     8,
     3,
     7,
     7,
     5,
     3,
     7,
     7,
     6,
     6,
     3,
     2,
     10,
     10,
     3,
     9,
     9,
     9,
     9,
     3,
     5,
     11,
     12,
     1,
     1,
     2,
     1,
     3,
     12,
     1 };
  static int y_vec[N_MAX] = {
    -4713,
    -4713,
    -4004,
    -3761,
    -3114,
    -3102,
    -2637,
     -776,
     -753,
     -747,
     -747,
     -312,
     -312,
      -23,
       -1,
        1,
        1,
        1,
        8,
       57,
       79,
      284,
      388,
      552,
      590,
      600,
      622,
      622,
      622,
      632,
      632,
     1078,
     1556,
     1582,
     1582,
     1676,
     1752,
     1752,
     1783,
     1792,
     1844,
     1844,
     1858,
     1899,
     1970,
     1980,
     1996,
     2000,
     2004,
     2012,
     3266 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    jed = 0.0;
    y = 0;
    m = 0;
    d = 0;
    f = 0.0;
  }
  else
  {
    jed = jed_vec[n_data-1];
    y = y_vec[n_data-1];
    m = m_vec[n_data-1];
    d = d_vec[n_data-1];
    f = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jed_mjd_values ( int &n_data, double &jed, double &mjd )

//****************************************************************************80
//
//  Purpose:
//
//    JED_MJD_VALUES returns the MJD for Julian Ephemeris Dates.
//
//  Discussion:
//
//    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
//    starting from noon on 1 January 4713 BCE.
//
//    The MJD (Modified Julian Day) counts days starting from midnight,
//    17 November 1858.  This essentially subtracts 2400000.5 days from the JED.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold and Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &JED, the Julian Ephemeris Date.
//
//    Output, double &MJD, the Modified Julian Ephemeris Date.
//
{
# define N_MAX 33

  static double jed_vec[N_MAX] = {
     1507231.5E+00,
     1660037.5E+00,
     1746893.5E+00,
     1770641.5E+00,
     1892731.5E+00,
     1931579.5E+00,
     1974851.5E+00,
     2091164.5E+00,
     2121509.5E+00,
     2155779.5E+00,
     2174029.5E+00,
     2191584.5E+00,
     2195261.5E+00,
     2229274.5E+00,
     2245580.5E+00,
     2266100.5E+00,
     2288542.5E+00,
     2290901.5E+00,
     2323140.5E+00,
     2334848.5E+00,
     2348020.5E+00,
     2366978.5E+00,
     2385648.5E+00,
     2392825.5E+00,
     2416223.5E+00,
     2425848.5E+00,
     2430266.5E+00,
     2430833.5E+00,
     2431004.5E+00,
     2448698.5E+00,
     2450138.5E+00,
     2465737.5E+00,
     2486076.5E+00 };

  static double mjd_vec[N_MAX] = {
     -892769.0E+00,
     -739963.0E+00,
     -653107.0E+00,
     -629359.0E+00,
     -507269.0E+00,
     -468421.0E+00,
     -425149.0E+00,
     -308836.0E+00,
     -278491.0E+00,
     -244221.0E+00,
     -225971.0E+00,
     -208416.0E+00,
     -204739.0E+00,
     -170726.0E+00,
     -154420.0E+00,
     -133900.0E+00,
     -111458.0E+00,
     -109099.0E+00,
      -76860.0E+00,
      -65152.0E+00,
      -51980.0E+00,
      -33022.0E+00,
      -14352.0E+00,
       -7175.0E+00,
       16223.0E+00,
       25848.0E+00,
       30266.0E+00,
       30833.0E+00,
       31004.0E+00,
       48698.0E+00,
       50138.0E+00,
       65737.0E+00,
       86076.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    jed = 0.0;
    mjd = 0.0;
  }
  else
  {
    jed = jed_vec[n_data-1];
    mjd = mjd_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jed_rd_values ( int &n_data, double &jed, double &rd )

//****************************************************************************80
//
//  Purpose:
//
//    JED_RD_VALUES returns the RD for Julian Ephemeris Dates.
//
//  Discussion:
//
//    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
//    starting from noon on 1 January 4713 BCE.
//
//    The RD is the Reingold Dershowitz Date, which counts days from
//    midnight, 1 January year 1 in the Gregorian calendar.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold and Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &JED, the Julian Ephemeris Date.
//
//    Output, double &RD, the Reingold Dershowitz Date.
//
{
# define N_MAX 33

  static double jed_vec[N_MAX] = {
    1507231.5E+00,
    1660037.5E+00,
    1746893.5E+00,
    1770641.5E+00,
    1892731.5E+00,
    1931579.5E+00,
    1974851.5E+00,
    2091164.5E+00,
    2121509.5E+00,
    2155779.5E+00,
    2174029.5E+00,
    2191584.5E+00,
    2195261.5E+00,
    2229274.5E+00,
    2245580.5E+00,
    2266100.5E+00,
    2288542.5E+00,
    2290901.5E+00,
    2323140.5E+00,
    2334848.5E+00,
    2348020.5E+00,
    2366978.5E+00,
    2385648.5E+00,
    2392825.5E+00,
    2416223.5E+00,
    2425848.5E+00,
    2430266.5E+00,
    2430833.5E+00,
    2431004.5E+00,
    2448698.5E+00,
    2450138.5E+00,
    2465737.5E+00,
    2486076.5E+00 };

  static double rd_vec[N_MAX] = {
    -214193.0E+00,
     -61387.0E+00,
      25469.0E+00,
      49217.0E+00,
     171307.0E+00,
     210155.0E+00,
     253427.0E+00,
     369740.0E+00,
     400085.0E+00,
     434355.0E+00,
     452605.0E+00,
     470160.0E+00,
     473837.0E+00,
     507850.0E+00,
     524156.0E+00,
     544676.0E+00,
     567118.0E+00,
     569477.0E+00,
     601716.0E+00,
     613424.0E+00,
     626596.0E+00,
     645554.0E+00,
     664224.0E+00,
     671401.0E+00,
     694799.0E+00,
     704424.0E+00,
     708842.0E+00,
     709409.0E+00,
     709580.0E+00,
     727274.0E+00,
     728714.0E+00,
     744313.0E+00,
     764652.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    jed = 0.0;
    rd = 0.0;
  }
  else
  {
    jed = jed_vec[n_data-1];
    rd = rd_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void jed_weekday_values ( int &n_data, double &jed, int &weekday )

//****************************************************************************80
//
//  Purpose:
//
//    JED_WEEKDAY_VALUES returns the day of the week for Julian Ephemeris Dates.
//
//  Discussion:
//
//    The JED (Julian Ephemeris Date) is a calendrical system which counts days,
//    starting from noon on 1 January 4713 BCE.
//
//    Weekdays are numbered as follows:
//
//    1  Sunday
//    2  Monday
//    3  Tuesday
//    4  Wednesday
//    5  Thursday
//    6  Friday
//    7  Saturday
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold and Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &JED, the Julian Ephemeris Date.
//
//    Output, int &WEEKDAY, the day of the week.
//
{
# define N_MAX 33

  static double jed_vec[N_MAX] = {
    1507231.5E+00,
    1660037.5E+00,
    1746893.5E+00,
    1770641.5E+00,
    1892731.5E+00,
    1931579.5E+00,
    1974851.5E+00,
    2091164.5E+00,
    2121509.5E+00,
    2155779.5E+00,
    2174029.5E+00,
    2191584.5E+00,
    2195261.5E+00,
    2229274.5E+00,
    2245580.5E+00,
    2266100.5E+00,
    2288542.5E+00,
    2290901.5E+00,
    2323140.5E+00,
    2334848.5E+00,
    2348020.5E+00,
    2366978.5E+00,
    2385648.5E+00,
    2392825.5E+00,
    2416223.5E+00,
    2425848.5E+00,
    2430266.5E+00,
    2430833.5E+00,
    2431004.5E+00,
    2448698.5E+00,
    2450138.5E+00,
    2465737.5E+00,
    2486076.5E+00 };

  static int weekday_vec[N_MAX] = {
    1, 4, 4, 1, 4,
    2, 7, 1, 1, 6,
    7, 6, 1, 1, 4,
    7, 7, 7, 4, 1,
    6, 1, 2, 4, 1,
    1, 2, 2, 5, 3,
    1, 4, 1 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    jed = 0.0;
    weekday = 0;
  }
  else
  {
    jed = jed_vec[n_data-1];
    weekday = weekday_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void kei0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    KEI0_VALUES returns some values of the Kelvin KEI function of order NU = 0.
//
//  Discussion:
//
//    The function is defined by:
//
//      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
//
//    where K(NU,X) is the K Bessel function.
//
//    In Mathematica, KEI(NU,X) can be defined by:
//
//      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
    0.0000000000000000,
    -0.6715816950943676,
    -0.4949946365187199,
    -0.3313955623385585,
    -0.2024000677647043,
    -0.1106960991556749,
   -0.05112188404598678,
   -0.01600256851827124,
   0.002198399294972520,
   0.009720918540151990,
   0.01118758650986964 };
  static double x_vec[N_MAX] = {
    0.0,
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void kei1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    KEI1_VALUES returns some values of the Kelvin KEI function of order NU = 1.
//
//  Discussion:
//
//    The function is defined by:
//
//      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
//
//    where K(NU,X) is the K Bessel function.
//
//    In Mathematica, KEI(NU,X) can be defined by:
//
//      Im [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 10

  static double fx_vec[N_MAX] = {
    -1.051182085412523,
    -0.2419959664297382,
    0.001008680985009855,
    0.08004939780706674,
    0.09331378813535750,
   0.08027022252392219,
   0.05937625647622691,
   0.03916601076917133,
   0.02300216024690250,
   0.01157775439325247 };
  static double x_vec[N_MAX] = {
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ker0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    KER0_VALUES returns some values of the Kelvin KER function of order NU = 0.
//
//  Discussion:
//
//    The function is defined by:
//
//      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
//
//    where K(NU,X) is the K Bessel function.
//
//    In Mathematica, KER(NU,X) can be defined by:
//
//      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 10

  static double fx_vec[N_MAX] = {
    0.8559058721186342,
    0.2867062087283160,
    0.05293491548771044,
   -0.04166451399150953,
   -0.06968797258904534,
   -0.06702923330379870,
   -0.05263927724224119,
   -0.03617884789954761,
   -0.02199987504667382,
   -0.01151172719949066 };
  static double x_vec[N_MAX] = {
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void ker1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    KER1_VALUES returns some values of the Kelvin KER function of order NU = 1.
//
//  Discussion:
//
//    The function is defined by:
//
//      KER(NU,X) + i * KEI(NU,X) = exp(-nu*Pi*I/2) * K(NU,X*exp(PI*I/4))
//
//    where K(NU,X) is the K Bessel function.
//
//    In Mathematica, KER(NU,X) can be defined by:
//
//      Re [ Exp [ -NU * Pi * I / 2 ] * BesselK [ NU, X * Exp[ Pi * I / 4 ] ] ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 10

  static double fx_vec[N_MAX] = {
   -1.522403406532090,
   -0.7403222768419827,
   -0.4170442851662574,
   -0.2308059295181230,
   -0.1172561358598705,
   -0.04989830778751491,
   -0.01272324936181659,
    0.005351296460277448,
    0.01209090413515866,
    0.01273739048421857 };
  static double x_vec[N_MAX] = {
    0.5,
    1.0,
    1.5,
    2.0,
    2.5,
    3.0,
    3.5,
    4.0,
    4.5,
    5.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void laguerre_associated_values ( int &n_data, int &n, int &m, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_ASSOCIATED_VALUES returns some values of the associated Laguerre polynomials.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,m,x]
//
//    The associated Laguerre polynomials may be generalized so that the
//    parameter M is allowed to take on arbitrary nonint *values.
//    The resulting function is known as the generalized Laguerre function.
//
//    The polynomials satisfy the differential equation:
//
//      X * Y'' + (M+1-X) * Y' + (N-M) * Y = 0;
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, int &M, the parameter.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1500000000000000E+01,
     0.1625000000000000E+01,
     0.1479166666666667E+01,
     0.1148437500000000E+01,
     0.4586666666666667E+00,
     0.2878666666666667E+01,
     0.8098666666666667E+01,
     0.1711866666666667E+02,
     0.1045328776041667E+02,
     0.1329019368489583E+02,
     0.5622453647189670E+02,
     0.7484729341779436E+02,
     0.3238912982762806E+03,
     0.4426100000097533E+03,
     0.1936876572288250E+04 };

  static int m_vec[N_MAX] = {
    0, 0, 0, 0,
    0, 1, 1, 1,
    1, 0, 1, 2,
    3, 2, 2, 3,
    3, 4, 4, 5 };

  static int n_vec[N_MAX] = {
    1,  2,  3,  4,
    5,  1,  2,  3,
    4,  3,  3,  3,
    3,  4,  5,  6,
    7,  8,  9, 10 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void laguerre_general_values ( int &n_data, int &n, double &a, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_GENERAL_VALUES returns some values of the generalized Laguerre function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,a,x]
//
//    The functions satisfy the following differential equation:
//
//      X * Y'' + (ALPHA+1-X) * Y' + N * Y = 0;
//
//    Function values can be generated by the recursion:
//
//      L(0,ALPHA)(X) = 1
//      L(1,ALPHA)(X) = 1+ALPHA-X
//
//      L(N,ALPHA)(X) = ( (2*N-1+ALPHA-X) * L(N-1,ALPHA)(X)
//                     - (N-1+ALPHA) * L(N-2,ALPHA)(X) ) / N
//
//    The parameter ALPHA is required to be greater than -1.
//
//    For ALPHA = 0, the generalized Laguerre function L(N,ALPHA)(X)
//    is equal to the Laguerre polynomial L(N)(X).
//
//    For ALPHA integral, the generalized Laguerre function
//    L(N,ALPHA)(X) equals the associated Laguerre polynomial L(N,ALPHA)(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &A, the parameter.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double a_vec[N_MAX] = {
     0.00E+00,
     0.25E+00,
     0.50E+00,
     0.75E+00,
     1.50E+00,
     2.50E+00,
     5.00E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     1.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00,
     5.20E+00 };

  static double fx_vec[N_MAX] = {
      0.3726399739583333E-01,
      0.3494791666666667E+00,
      0.8710042317708333E+00,
      0.1672395833333333E+01,
      0.6657625325520833E+01,
      0.2395726725260417E+02,
      0.2031344319661458E+03,
      0.1284193996800000E+02,
      0.5359924801587302E+01,
      0.9204589064126984E+00,
     -0.1341585114857143E+01,
     -0.2119726307555556E+01,
     -0.1959193658349206E+01,
      0.1000000000000000E+01,
      0.5450000000000000E+01,
      0.1720125000000000E+02,
      0.4110393750000000E+02,
      0.8239745859375000E+02,
      0.1460179186171875E+03,
      0.2359204608298828E+03 };

  static int n_vec[N_MAX] = {
     5,
     5,
     5,
     5,
     5,
     5,
     5,
     8,
     8,
     8,
     8,
     8,
     8,
     0,
     1,
     2,
     3,
     4,
     5,
     6 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.20E+00,
     0.40E+00,
     0.60E+00,
     0.80E+00,
     1.00E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00,
     0.75E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    a = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void laguerre_polynomial_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAGUERRE_POLYNOMIAL_VALUES returns some values of the Laguerre polynomial.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LaguerreL[n,x]
//
//  Differential equation:
//
//    X * Y'' + (1-X) * Y' + N * Y = 0;
//
//  First terms:
//
//      1
//     -X    +  1
//   (  X^2 -  4 X     +  2 ) / 2
//   ( -X^3 +  9 X^2 -  18 X    +    6 ) / 6
//   (  X^4 - 16 X^3 +  72 X^2 -   96 X +      24 ) / 24
//   ( -X^5 + 25 X^4 - 200 X^3 +  600 X^2 -  600 x    +  120 ) / 120
//   (  X^6 - 36 X^5 + 450 X^4 - 2400 X^3 + 5400 X^2 - 4320 X + 720 ) / 720
//   ( -X^7 + 49 X^6 - 882 X^5 + 7350 X^4 - 29400 X^3
//      + 52920 X^2 - 35280 X + 5040 ) / 5040
//
//  Recursion:
//
//    L(0)(X) = 1,
//    L(1)(X) = 1-X,
//    N * L(N)(X) = (2*N-1-X) * L(N-1)(X) - (N-1) * L(N-2)(X)
//
//  Orthogonality:
//
//    Integral ( 0 <= X < +oo ) exp ( - X ) * L(N)(X) * L(M)(X) dX
//    = 0 if N /= M
//    = 1 if N == M
//
//  Special values:
//
//    L(N)(0) = 1.
//
//  Relations:
//
//    L(N)(X) = (-1)**N / N! * exp ( x ) * (d/dx)**n ( exp ( - x ) * x**n )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the polynomial.
//
//    Output, double &X, the point where the polynomial is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 17

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.0000000000000000E+00,
     -0.5000000000000000E+00,
     -0.6666666666666667E+00,
     -0.6250000000000000E+00,
     -0.4666666666666667E+00,
     -0.2569444444444444E+00,
     -0.4047619047619048E-01,
      0.1539930555555556E+00,
      0.3097442680776014E+00,
      0.4189459325396825E+00,
      0.4801341790925124E+00,
      0.4962122235082305E+00,
     -0.4455729166666667E+00,
      0.8500000000000000E+00,
     -0.3166666666666667E+01,
      0.3433333333333333E+02  };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10, 11,
    12,  5,  5,
     5,  5 };

  static double x_vec[N_MAX] = {
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     0.5E+00,
     3.0E+00,
     5.0E+00,
     1.0E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void lambert_w_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAMBERT_W_VALUES returns some values of the Lambert W function.
//
//  Discussion:
//
//    The function W(X) is defined implicitly by:
//
//      W(X) * e^W(X) = X
//
//    The function is also known as the "Omega" function.
//
//    In Mathematica, the function can be evaluated by:
//
//      W = ProductLog [ X ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Brian Hayes,
//    "Why W?",
//    The American Scientist,
//    Volume 93, March-April 2005, pages 104-108.
//
//    Eric Weisstein,
//    "Lambert's W-Function",
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 1998.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
    0.0000000000000000E+00,
    0.3517337112491958E+00,
    0.5671432904097839E+00,
    0.7258613577662263E+00,
    0.8526055020137255E+00,
    0.9585863567287029E+00,
    0.1000000000000000E+01,
    0.1049908894964040E+01,
    0.1130289326974136E+01,
    0.1202167873197043E+01,
    0.1267237814307435E+01,
    0.1326724665242200E+01,
    0.1381545379445041E+01,
    0.1432404775898300E+01,
    0.1479856830173851E+01,
    0.1524345204984144E+01,
    0.1566230953782388E+01,
    0.1605811996320178E+01,
    0.1745528002740699E+01,
    0.3385630140290050E+01,
    0.5249602852401596E+01,
    0.1138335808614005E+02 };
  static double x_vec[N_MAX] = {
    0.0000000000000000E+00,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.1500000000000000E+01,
    0.2000000000000000E+01,
    0.2500000000000000E+01,
    0.2718281828459045E+01,
    0.3000000000000000E+01,
    0.3500000000000000E+01,
    0.4000000000000000E+01,
    0.4500000000000000E+01,
    0.5000000000000000E+01,
    0.5500000000000000E+01,
    0.6000000000000000E+01,
    0.6500000000000000E+01,
    0.7000000000000000E+01,
    0.7500000000000000E+01,
    0.8000000000000000E+01,
    0.1000000000000000E+02,
    0.1000000000000000E+03,
    0.1000000000000000E+04,
    0.1000000000000000E+07 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void laplace_cdf_values ( int &n_data, double &mu, double &beta, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LAPLACE_CDF_VALUES returns some values of the Laplace CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LaplaceDistribution [ mu, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &BETA, the shape parameter.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double beta_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.8160602794142788E+00,
     0.9323323583816937E+00,
     0.9751064658160680E+00,
     0.6967346701436833E+00,
     0.6417343447131054E+00,
     0.6105996084642976E+00,
     0.5906346234610091E+00,
     0.5000000000000000E+00,
     0.3032653298563167E+00,
     0.1839397205857212E+00,
     0.1115650800742149E+00 };

  static double mu_vec[N_MAX] = {
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.0000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.0000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    beta = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    beta = beta_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated_values ( int &n_data, int &n, int &m, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_VALUES returns values of associated Legendre functions.
//
//  Discussion:
//
//    The function considered is the associated Legendre polynomial P^M_N(X).
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//  Differential equation:
//
//    (1-X*X) * Y'' - 2 * X * Y + ( N (N+1) - (M*M/(1-X*X)) * Y = 0;
//
//  First terms:
//
//    M = 0  ( = Legendre polynomials of first kind P(N)(X) )
//
//    P00 =    1
//    P10 =    1 X
//    P20 = (  3 X^2 -   1)/2
//    P30 = (  5 X^3 -   3 X)/2
//    P40 = ( 35 X^4 -  30 X^2 +   3)/8
//    P50 = ( 63 X^5 -  70 X^3 +  15 X)/8
//    P60 = (231 X^6 - 315 X^4 + 105 X^2 -  5)/16
//    P70 = (429 X^7 - 693 X^5 + 315 X^3 - 35 X)/16
//
//    M = 1
//
//    P01 =   0
//    P11 =   1 * SQRT(1-X*X)
//    P21 =   3 * SQRT(1-X*X) * X
//    P31 = 1.5 * SQRT(1-X*X) * (5*X*X-1)
//    P41 = 2.5 * SQRT(1-X*X) * (7*X*X*X-3*X)
//
//    M = 2
//
//    P02 =   0
//    P12 =   0
//    P22 =   3 * (1-X*X)
//    P32 =  15 * (1-X*X) * X
//    P42 = 7.5 * (1-X*X) * (7*X*X-1)
//
//    M = 3
//
//    P03 =   0
//    P13 =   0
//    P23 =   0
//    P33 =  15 * (1-X*X)^1.5
//    P43 = 105 * (1-X*X)^1.5 * X
//
//    M = 4
//
//    P04 =   0
//    P14 =   0
//    P24 =   0
//    P34 =   0
//    P44 = 105 * (1-X*X)^2
//
//  Recursion:
//
//    if N < M:
//      P(N,M) = 0;
//    if N = M:
//      P(N,M) = (2*M-1)!! * (1-X*X)**(M/2) where N!! means the product of
//      all the odd integers less than or equal to N.
//    if N = M+1:
//      P(N,M) = X*(2*M+1)*P(M,M)
//    if M+1 < N:
//      P(N,M) = ( X*(2*N-1)*P(N-1,M) - (N+M-1)*P(N-2,M) )/(N-M)
//
//  Restrictions:
//
//    -1 <= X <= 1
//     0 <= M <= N
//
//  Special values:
//
//    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
//    polynomial of the first kind equals the Legendre polynomial of the
//    first kind.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.0000000000000000E+00,
     -0.5000000000000000E+00,
      0.0000000000000000E+00,
      0.3750000000000000E+00,
      0.0000000000000000E+00,
     -0.8660254037844386E+00,
     -0.1299038105676658E+01,
     -0.3247595264191645E+00,
      0.1353164693413185E+01,
     -0.2800000000000000E+00,
      0.1175755076535925E+01,
      0.2880000000000000E+01,
     -0.1410906091843111E+02,
     -0.3955078125000000E+01,
     -0.9997558593750000E+01,
      0.8265311444100484E+02,
      0.2024442836815152E+02,
     -0.4237997531890869E+03,
      0.1638320624828339E+04,
     -0.2025687389227225E+05  };

  static int m_vec[N_MAX] = {
    0, 0, 0, 0,
    0, 1, 1, 1,
    1, 0, 1, 2,
    3, 2, 2, 3,
    3, 4, 4, 5 };

  static int n_vec[N_MAX] = {
    1,  2,  3,  4,
    5,  1,  2,  3,
    4,  3,  3,  3,
    3,  4,  5,  6,
    7,  8,  9, 10 };

  static double x_vec[N_MAX] = {
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.00E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.20E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated_normalized_sphere_values ( int &n_data, int &n, int &m,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_SPHERE_VALUES: normalized associated Legendre.
//
//  Discussion:
//
//    The function considered is the associated Legendre polynomial P^M_N(X).
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//    The function is normalized for the sphere by dividing by
//
//      sqrt ( 4 * pi * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
     0.2820947917738781,
     0.2443012559514600,
    -0.2992067103010745,
    -0.07884789131313000,
    -0.3345232717786446,
     0.2897056515173922,
    -0.3265292910163510,
    -0.06997056236064664,
     0.3832445536624809,
    -0.2709948227475519,
    -0.2446290772414100,
     0.2560660384200185,
     0.1881693403754876,
    -0.4064922341213279,
     0.2489246395003027,
     0.08405804426339821,
     0.3293793022891428,
    -0.1588847984307093,
    -0.2808712959945307,
     0.4127948151484925,
    -0.2260970318780046 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_associated_normalized_values ( int &n_data, int &n, int &m,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_ASSOCIATED_NORMALIZED_VALUES: normalized associated Legendre.
//
//  Discussion:
//
//    The function considered is the associated Legendre polynomial P^M_N(X).
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreP [ n, m, x ]
//
//    The function is normalized by dividing by
//
//      sqrt ( 2 * ( n + m )! / ( 2 * n + 1 ) / ( n - m )! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 March 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, int &M, double &X,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
    0.7071067811865475E+00, 
    0.6123724356957945E+00, 
   -0.7500000000000000E+00, 
   -0.1976423537605237E+00, 
   -0.8385254915624211E+00, 
    0.7261843774138907E+00, 
   -0.8184875533567997E+00, 
   -0.1753901900050285E+00, 
    0.9606516343087123E+00, 
   -0.6792832849776299E+00, 
   -0.6131941618102092E+00, 
    0.6418623720763665E+00, 
    0.4716705890038619E+00, 
   -0.1018924927466445E+01, 
    0.6239615396237876E+00, 
    0.2107022704608181E+00, 
    0.8256314721961969E+00, 
   -0.3982651281554632E+00, 
   -0.7040399320721435E+00, 
    0.1034723155272289E+01, 
   -0.5667412129155530E+00 };

  static int m_vec[N_MAX] = {
    0, 0, 1, 0,
    1, 2, 0, 1,
    2, 3, 0, 1,
    2, 3, 4, 0,
    1, 2, 3, 4,
    5 };

  static int n_vec[N_MAX] = {
    0,  1,  1,  2,
    2,  2,  3,  3,
    3,  3,  4,  4,
    4,  4,  4,  5,
    5,  5,  5,  5,
    5 };

  static double x_vec[N_MAX] = {
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50,
    0.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_function_q_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_FUNCTION_Q_VALUES returns values of the Legendre Q function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreQ[n,x]
//
//  Differential equation:
//
//    (1-X*X) Y'' - 2 X Y' + N (N+1) = 0;
//
//  First terms:
//
//    Q(0)(X) = 0.5 * log((1+X)/(1-X))
//    Q(1)(X) = Q(0)(X)*X - 1
//    Q(2)(X) = Q(0)(X)*(3*X*X-1)/4 - 1.5*X
//    Q(3)(X) = Q(0)(X)*(5*X*X*X-3*X)/4 - 2.5*X^2 + 2/3
//    Q(4)(X) = Q(0)(X)*(35*X^4-30*X^2+3)/16 - 35/8 * X^3 + 55/24 * X
//    Q(5)(X) = Q(0)(X)*(63*X^5-70*X^3+15*X)/16 - 63/8*X^4 + 49/8*X^2 - 8/15
//
//  Recursion:
//
//    Q(0) = 0.5 * log ( (1+X) / (1-X) )
//    Q(1) = 0.5 * X * log ( (1+X) / (1-X) ) - 1.0
//
//    Q(N) = ( (2*N-1) * X * Q(N-1) - (N-1) * Q(N-2) ) / N
//
//  Restrictions:
//
//    -1 < X < 1
//
//  Special values:
//
//    Note that the Legendre function Q(N)(X) is equal to the
//    associated Legendre function of the second kind,
//    Q(N,M)(X) with M = 0.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double fx_vec[N_MAX] = {
      0.2554128118829953E+00,
     -0.9361467970292512E+00,
     -0.4787614548274669E+00,
      0.4246139251747229E+00,
      0.5448396833845414E+00,
     -0.9451328261673470E-01,
     -0.4973516573531213E+00,
     -0.1499018843853194E+00,
      0.3649161918783626E+00,
      0.3055676545072885E+00,
     -0.1832799367995643E+00,
      0.6666666666666667E+00,
      0.6268672028763330E+00,
      0.5099015515315237E+00,
      0.3232754180589764E+00,
      0.8026113738148187E-01,
     -0.1986547714794823E+00,
     -0.4828663183349136E+00,
     -0.7252886849144386E+00,
     -0.8454443502398846E+00,
     -0.6627096245052618E+00 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void legendre_poly_values ( int &n_data, int &n, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_POLY_VALUES returns values of the Legendre polynomials.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      LegendreP [ n, x ]
//
//  Differential equation:
//
//    (1-X*X) * P(N)(X)'' - 2 * X * P(N)(X)' + N * (N+1) = 0;
//
//  First terms:
//
//    P( 0)(X) =       1
//    P( 1)(X) =       1 X
//    P( 2)(X) =  (    3 X^2 -       1)/2
//    P( 3)(X) =  (    5 X^3 -     3 X)/2
//    P( 4)(X) =  (   35 X^4 -    30 X^2 +     3)/8
//    P( 5)(X) =  (   63 X^5 -    70 X^3 +    15 X)/8
//    P( 6)(X) =  (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
//    P( 7)(X) =  (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
//    P( 8)(X) =  ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
//    P( 9)(X) =  (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
//    P(10)(X) =  (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2
//                 -63 ) /256
//
//  Recursion:
//
//    P(0)(X) = 1
//    P(1)(X) = X
//    P(N)(X) = ( (2*N-1)*X*P(N-1)(X)-(N-1)*P(N-2)(X) ) / N
//
//    P'(0)(X) = 0;
//    P'(1)(X) = 1
//    P'(N)(X) = ( (2*N-1)*(P(N-1)(X)+X*P'(N-1)(X)-(N-1)*P'(N-2)(X) ) / N
//
//  Formula:
//
//    P(N)(X) = (1/2**N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
//
//  Orthogonality:
//
//    Integral ( -1 <= X <= 1 ) P(I)(X) * P(J)(X) dX
//      = 0 if I =/= J
//      = 2 / ( 2*I+1 ) if I = J.
//
//  Approximation:
//
//    A function F(X) defined on [-1,1] may be approximated by the series
//
//      C0*P(0)(X) + C1*P(1)(X) + ... + CN*P(N)(X)
//
//    where
//
//      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I)(X) dx.
//
//  Special values:
//
//    P(N)(1) = 1.
//    P(N)(-1) = (-1)**N.
//    | P(N)(X) | <= 1 in [-1,1].
//
//    P(N,0)(X) = P(N)(X), that is, for M=0, the associated Legendre
//    function of the first kind and order N equals the Legendre polynomial
//    of the first kind and order N.
//
//    The N zeroes of P(N)(X) are the abscissas used for Gauss-Legendre
//    quadrature of the integral of a function F(X) with weight function 1
//    over the interval [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the order of the function.
//
//    Output, double &X, the point where the function is evaluated.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
      0.1000000000000000E+01,
      0.2500000000000000E+00,
     -0.4062500000000000E+00,
     -0.3359375000000000E+00,
      0.1577148437500000E+00,
      0.3397216796875000E+00,
      0.2427673339843750E-01,
     -0.2799186706542969E+00,
     -0.1524540185928345E+00,
      0.1768244206905365E+00,
      0.2212002165615559E+00,
      0.0000000000000000E+00,
     -0.1475000000000000E+00,
     -0.2800000000000000E+00,
     -0.3825000000000000E+00,
     -0.4400000000000000E+00,
     -0.4375000000000000E+00,
     -0.3600000000000000E+00,
     -0.1925000000000000E+00,
      0.8000000000000000E-01,
      0.4725000000000000E+00,
      0.1000000000000000E+01 };

  static int n_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     6,  7,  8,
     9, 10,  3,
     3,  3,  3,
     3,  3,  3,
     3,  3,  3,
     3 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.00E+00,
     0.10E+00,
     0.20E+00,
     0.30E+00,
     0.40E+00,
     0.50E+00,
     0.60E+00,
     0.70E+00,
     0.80E+00,
     0.90E+00,
     1.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void lerch_values ( int &n_data, double &z, int &s, double &a, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LERCH_VALUES returns some values of the Lerch transcendent function.
//
//  Discussion:
//
//    The Lerch function is defined as
//
//      Phi(z,s,a) = Sum ( 0 <= k < +oo ) z^k / ( a + k )^s
//
//    omitting any terms with ( a + k ) = 0.
//
//    In Mathematica, the function can be evaluated by:
//
//      LerchPhi[z,s,a]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &Z, the parameters of the function.
//
//    Output, int &S, the parameters of the function.
//
//    Output, double &A, the parameters of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double a_vec[N_MAX] = {
    0.0E+00,
    0.0E+00,
    0.0E+00,
    1.0E+00,
    1.0E+00,
    1.0E+00,
    2.0E+00,
    2.0E+00,
    2.0E+00,
    3.0E+00,
    3.0E+00,
    3.0E+00 };

  static double fx_vec[N_MAX] = {
    0.1644934066848226E+01,
    0.1202056903159594E+01,
    0.1000994575127818E+01,
    0.1164481052930025E+01,
    0.1074426387216080E+01,
    0.1000492641212014E+01,
    0.2959190697935714E+00,
    0.1394507503935608E+00,
    0.9823175058446061E-03,
    0.1177910993911311E+00,
    0.3868447922298962E-01,
    0.1703149614186634E-04 };

  static int s_vec[N_MAX] = {
     2, 3, 10,
     2, 3, 10,
     2, 3, 10,
     2, 3, 10 };

  static double z_vec[N_MAX] = {
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.3333333333333333E+00,
    0.3333333333333333E+00,
    0.3333333333333333E+00,
    0.1000000000000000E+00,
    0.1000000000000000E+00,
    0.1000000000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    z = 0.0;
    s = 0;
    a = 0.0;
    fx = 0.0;
  }
  else
  {
    z = z_vec[n_data-1];
    s = s_vec[n_data-1];
    a = a_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void linear_system_values ( int &n_data, int &nrow, int &ncol, int &nsys,
  double *a, double *x, double *b )

//****************************************************************************80
//
//  Purpose:
//
//    LINEAR_SYSTEM_VALUES returns some linear systems.
//
//  Discussion:
//
//    Each call to this routine returns scalars NROW, NCOL and NSYS,
//    which give the dimensions of the linear system
//
//      A(NROW,NCOL) * X(NCOL,NSYS) = B(NROW,NSYS)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 May 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &NROW, int NCOL, the number of rows and columns of A.
//
//    Output, int &NSYS, the number of systems.
//
//    Output, double **A[NROW*NCOL], the matrix.
//
//    Output, double **X[NCOL*NSYS], the solutions of the linear system.
//
//    Output, double **B[NROW*NSYS], the right hand sides.
//
{
# define N_MAX 4

  if ( a )
  {
    delete [] a;
  }
  if ( b )
  {
    delete [] b;
  }
  if ( x )
  {
    delete [] x;
  }

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    nrow = 0;
    ncol = 0;
    nsys = 0;
    a = NULL;
    b = NULL;
    x = NULL;
  }
  else if ( n_data == 1 )
  {
    nrow = 3;
    ncol = 3;
    nsys = 2;

    a = new double[(nrow)*(ncol)];
    x = new double[(ncol)*(nsys)];
    b = new double[(nrow)*(nsys)];

    a[0+0*(nrow)] = 1.0;
    a[1+0*(nrow)] = 0.0;
    a[2+0*(nrow)] = 0.0;
    a[0+1*(nrow)] = 0.0;
    a[1+1*(nrow)] = 2.0;
    a[2+1*(nrow)] = 0.0;
    a[0+2*(nrow)] = 0.0;
    a[1+2*(nrow)] = 0.0;
    a[2+2*(nrow)] = 3.0;

    x[0+0*(ncol)] = 1.0;
    x[1+0*(ncol)] = 0.0;
    x[2+0*(ncol)] = 0.0;
    x[0+1*(ncol)] = 1.0;
    x[1+1*(ncol)] = 1.0;
    x[2+1*(ncol)] = 1.0;

    b[0+0*(nrow)] = 1.0;
    b[1+0*(nrow)] = 0.0;
    b[2+0*(nrow)] = 0.0;
    b[0+1*(nrow)] = 1.0;
    b[1+1*(nrow)] = 2.0;
    b[2+1*(nrow)] = 3.0;
  }
  else if ( n_data == 2 )
  {
    nrow = 3;
    ncol = 3;
    nsys = 2;

    a = new double[(nrow)*(ncol)];
    x = new double[(ncol)*(nsys)];
    b = new double[(nrow)*(nsys)];

    a[0+0*(nrow)] = 1.0;
    a[1+0*(nrow)] = 2.0;
    a[2+0*(nrow)] = 3.0;
    a[0+1*(nrow)] = 2.0;
    a[1+1*(nrow)] = 2.0;
    a[2+1*(nrow)] = 3.0;
    a[0+2*(nrow)] = 3.0;
    a[1+2*(nrow)] = 3.0;
    a[2+2*(nrow)] = 3.0;

    x[0+0*(ncol)] = 1.0;
    x[1+0*(ncol)] = 1.0;
    x[2+0*(ncol)] = 1.0;
    x[0+1*(ncol)] = 1.0;
    x[1+1*(ncol)] = 2.0;
    x[2+1*(ncol)] = 3.0;

    b[0+0*(nrow)] =  6.0;
    b[1+0*(nrow)] =  7.0;
    b[2+0*(nrow)] =  9.0;
    b[0+1*(nrow)] = 14.0;
    b[1+1*(nrow)] = 15.0;
    b[2+1*(nrow)] = 18.0;
  }
  else if ( n_data == 3 )
  {
    nrow = 5;
    ncol = 5;
    nsys = 2;

    a = new double[(nrow)*(ncol)];
    x = new double[(ncol)*(nsys)];
    b = new double[(nrow)*(nsys)];

    a[0+0*(nrow)] = 1.0;
    a[1+0*(nrow)] = 2.0;
    a[2+0*(nrow)] = 3.0;
    a[3+0*(nrow)] = 4.0;
    a[4+0*(nrow)] = 5.0;
    a[0+1*(nrow)] = 2.0;
    a[1+1*(nrow)] = 3.0;
    a[2+1*(nrow)] = 4.0;
    a[3+1*(nrow)] = 5.0;
    a[4+1*(nrow)] = 1.0;
    a[0+2*(nrow)] = 3.0;
    a[1+2*(nrow)] = 4.0;
    a[2+2*(nrow)] = 5.0;
    a[3+2*(nrow)] = 1.0;
    a[4+2*(nrow)] = 2.0;
    a[0+3*(nrow)] = 4.0;
    a[1+3*(nrow)] = 5.0;
    a[2+3*(nrow)] = 1.0;
    a[3+3*(nrow)] = 2.0;
    a[4+3*(nrow)] = 3.0;
    a[0+4*(nrow)] = 5.0;
    a[1+4*(nrow)] = 1.0;
    a[2+4*(nrow)] = 2.0;
    a[3+4*(nrow)] = 3.0;
    a[4+4*(nrow)] = 4.0;

    x[0+0*(ncol)] = 0.066667;
    x[1+0*(ncol)] = 0.066667;
    x[2+0*(ncol)] = 0.066667;
    x[3+0*(ncol)] = 0.066667;
    x[4+0*(ncol)] = 0.066667;
    x[0+1*(ncol)] = 1.0;
    x[1+1*(ncol)] = 0.0;
    x[2+1*(ncol)] = 0.0;
    x[3+1*(ncol)] = 0.0;
    x[4+1*(ncol)] = 0.0;

    b[0+0*(nrow)] = 1.0;
    b[1+0*(nrow)] = 1.0;
    b[2+0*(nrow)] = 1.0;
    b[3+0*(nrow)] = 1.0;
    b[4+0*(nrow)] = 1.0;
    b[0+1*(nrow)] = 1.0;
    b[1+1*(nrow)] = 2.0;
    b[2+1*(nrow)] = 3.0;
    b[3+1*(nrow)] = 4.0;
    b[4+1*(nrow)] = 5.0;
  }
  else if ( n_data == 4 )
  {
    nrow = 5;
    ncol = 5;
    nsys = 2;

    a = new double[(nrow)*(ncol)];
    x = new double[(ncol)*(nsys)];
    b = new double[(nrow)*(nsys)];

    a[0+0*(nrow)] = 1.4;
    a[1+0*(nrow)] = 1.6;
    a[2+0*(nrow)] = 3.8;
    a[3+0*(nrow)] = 4.6;
    a[4+0*(nrow)] = 2.6;
    a[0+1*(nrow)] = 2.1;
    a[1+1*(nrow)] = 1.5;
    a[2+1*(nrow)] = 8.0;
    a[3+1*(nrow)] = 8.2;
    a[4+1*(nrow)] = 2.9;
    a[0+2*(nrow)] = 2.1;
    a[1+2*(nrow)] = 1.1;
    a[2+2*(nrow)] = 9.6;
    a[3+2*(nrow)] = 8.4;
    a[4+2*(nrow)] = 0.1;
    a[0+3*(nrow)] = 7.4;
    a[1+3*(nrow)] = 0.7;
    a[2+3*(nrow)] = 5.4;
    a[3+3*(nrow)] = 0.4;
    a[4+3*(nrow)] = 9.6;
    a[0+4*(nrow)] = 9.6;
    a[1+4*(nrow)] = 5.0;
    a[2+4*(nrow)] = 8.8;
    a[3+4*(nrow)] = 8.0;
    a[4+4*(nrow)] = 7.7;

    x[0+0*(ncol)] =  -5.313077;
    x[1+0*(ncol)] =   5.735670;
    x[2+0*(ncol)] =  -2.507606;
    x[3+0*(ncol)] =  -1.058741;
    x[4+0*(ncol)] =   0.999381;
    x[0+1*(ncol)] =  31.601006;
    x[1+1*(ncol)] = -28.594793;
    x[2+1*(ncol)] =  13.389395;
    x[3+1*(ncol)] =   2.780322;
    x[4+1*(ncol)] =  -3.008797;

    b[0+0*(nrow)] = 1.1;
    b[1+0*(nrow)] = 1.6;
    b[2+0*(nrow)] = 4.7;
    b[3+0*(nrow)] = 9.1;
    b[4+0*(nrow)] = 0.1;
    b[0+1*(nrow)] = 4.0;
    b[1+1*(nrow)] = 9.3;
    b[2+1*(nrow)] = 8.4;
    b[3+1*(nrow)] = 0.4;
    b[4+1*(nrow)] = 4.1;
  }
  return;
# undef N_MAX
}
//****************************************************************************80

void lobachevsky_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOBACHEVSKY_VALUES returns some values of the Lobachevsky function.
//
//  Discussion:
//
//    The function is defined by:
//
//      LOBACHEVSKY(x) = Integral ( 0 <= t <= x ) -ln ( abs ( cos ( t ) ) dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.12417639065161393857E-08,
     0.79473344770001088225E-07,
     0.50867598186208834198E-05,
     0.32603097901207200319E-03,
     0.21380536815408214419E-01,
     0.18753816902083824050E+00,
     0.83051199971883645115E+00,
     0.18854362426679034904E+01,
     0.21315988986516411053E+01,
     0.21771120185613427221E+01,
     0.22921027921896650849E+01,
     0.39137195028784495586E+01,
     0.43513563983836427904E+01,
     0.44200644968478185898E+01,
     0.65656013133623829156E+01,
     0.10825504661504599479E+02,
     0.13365512855474227325E+02,
     0.21131002685639959927E+02,
     0.34838236589449117389E+02,
     0.69657062437837394278E+02 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0078125000E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       5.0000000000E+00,
       6.0000000000E+00,
       7.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00,
     100.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void log_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_VALUES returns some values of the natural logarithm function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
    -11.512925464970228420E+00,
     -4.6051701859880913680E+00,
     -2.3025850929940456840E+00,
     -1.6094379124341003746E+00,
     -1.2039728043259359926E+00,
     -0.91629073187415506518E+00,
     -0.69314718055994530942E+00,
     -0.51082562376599068321E+00,
     -0.35667494393873237891E+00,
     -0.22314355131420975577E+00,
     -0.10536051565782630123E+00,
      0.00000000000000000000E+00,
      0.69314718055994530942E+00,
      1.0986122886681096914E+00,
      1.1447298858494001741E+00,
      1.6094379124341003746E+00,
      2.3025850929940456840E+00,
      2.9957322735539909934E+00,
      4.6051701859880913680E+00,
      18.631401766168018033E+00 };

  static double x_vec[N_MAX] = {
    1.0E-05,
    1.0E-02,
    0.1E+00,
    0.2E+00,
    0.3E+00,
    0.4E+00,
    0.5E+00,
    0.6E+00,
    0.7E+00,
    0.8E+00,
    0.9E+00,
    1.0E+00,
    2.0E+00,
    3.0E+00,
    3.1415926535897932385E+00,
    5.0E+00,
    10.0E+00,
    20.0E+00,
    100.0E+00,
    123456789.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void log_normal_cdf_values ( int &n_data, double &mu, double &sigma,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_NORMAL_CDF_VALUES returns some values of the Log Normal CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LogNormalDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &SIGMA, the shape parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.2275013194817921E-01,
     0.2697049307349095E+00,
     0.5781741008028732E+00,
     0.7801170895122241E+00,
     0.4390310097476894E+00,
     0.4592655190218048E+00,
     0.4694258497695908E+00,
     0.4755320473858733E+00,
     0.3261051056816658E+00,
     0.1708799040927608E+00,
     0.7343256357952060E-01,
     0.2554673736161761E-01 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void log_series_cdf_values ( int &n_data, double &t, int &n, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG_SERIES_CDF_VALUES returns some values of the log series CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = LogSeriesDistribution [ t ]
//      CDF [ dist, n ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    27 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &T, the parameter of the function.
//
//    Output, int &N, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 29

  static double fx_vec[N_MAX] = {
     0.9491221581029903E+00,
     0.9433541128559735E+00,
     0.9361094611773272E+00,
     0.9267370278044118E+00,
     0.9141358246245129E+00,
     0.8962840235449100E+00,
     0.8690148741955517E+00,
     0.8221011541254772E+00,
     0.7213475204444817E+00,
     0.6068261510845583E+00,
     0.5410106403333613E+00,
     0.4970679476476894E+00,
     0.4650921887927060E+00,
     0.4404842934597863E+00,
     0.4207860535926143E+00,
     0.4045507673897055E+00,
     0.3908650337129266E+00,
     0.2149757685421097E+00,
     0.0000000000000000E+00,
     0.2149757685421097E+00,
     0.3213887739704539E+00,
     0.3916213575531612E+00,
     0.4437690508633213E+00,
     0.4850700239649681E+00,
     0.5191433267738267E+00,
     0.5480569580144867E+00,
     0.5731033910767085E+00,
     0.5951442521714636E+00,
     0.6147826594068904E+00 };

  static int n_vec[N_MAX] = {
     1, 1, 1, 1, 1,
     1, 1, 1, 1, 1,
     1, 1, 1, 1, 1,
     1, 1, 1, 0, 1,
     2, 3, 4, 5, 6,
     7, 8, 9, 10 };

  static double t_vec[N_MAX] = {
     0.1000000000000000E+00,
     0.1111111111111111E+00,
     0.1250000000000000E+00,
     0.1428571428571429E+00,
     0.1666666666666667E+00,
     0.2000000000000000E+00,
     0.2500000000000000E+00,
     0.3333333333333333E+00,
     0.5000000000000000E+00,
     0.6666666666666667E+00,
     0.7500000000000000E+00,
     0.8000000000000000E+00,
     0.8333333333333333E+00,
     0.8571485714857149E+00,
     0.8750000000000000E+00,
     0.8888888888888889E+00,
     0.9000000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00,
     0.9900000000000000E+00  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    t = 0.0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    t = t_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void log10_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOG10_VALUES returns some values of the logarithm base 10 function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
   -5.0000000000000000000,
   -2.0000000000000000000,
   -1.0000000000000000000,
   -0.69897000433601880479,
   -0.52287874528033756270,
   -0.39794000867203760957,
   -0.30102999566398119521,
   -0.22184874961635636749,
   -0.15490195998574316929,
   -0.096910013008056414359,
   -0.045757490560675125410,
    0.000000000000000000000,
    0.30102999566398119521,
    0.47712125471966243730,
    0.49714987269413385435,
    0.69897000433601880479,
    1.0000000000000000000,
    1.3010299956639811952,
    2.0000000000000000000,
    8.0915149771692704475 };

  static double x_vec[N_MAX] = {
    1.0E-05,
    1.0E-02,
    0.1E+00,
    0.2E+00,
    0.3E+00,
    0.4E+00,
    0.5E+00,
    0.6E+00,
    0.7E+00,
    0.8E+00,
    0.9E+00,
    1.0E+00,
    2.0E+00,
    3.0E+00,
    3.1415926535897932385E+00,
    5.0E+00,
    10.0E+00,
    20.0E+00,
    100.0E+00,
    123456789.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void logarithmic_integral_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOGARITHMIC_INTEGRAL_VALUES returns values of the logarithmic integral LI(X).
//
//  Discussion:
//
//    The logarithmic integral is defined as:
//
//      LI(X) = integral ( 0 <= T <= Z ) dT / log ( T )
//
//    The principal value of the integral is taken.  There is a
//    branch cut discontinuity in the complex plane from -oo to +1.
//
//    Abramowitz and Stegun assume 1 < X.
//
//    In Mathematica, the function can be evaluated by:
//
//      LogIntegral[x]
//
//    There is a simple relationship with the exponential integral EI:
//
//      LI(X) = EI(LN(X))
//
//    The function LI(X) provides a good approximation to PI(X),
//    the number of primes less than or equal to X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 28

  static double fx_vec[N_MAX] = {
      0.0000000000000000E+00,
     -0.3238978959329102E-01,
     -0.8512648672879405E-01,
     -0.1574149028946895E+00,
     -0.2529494192126213E+00,
     -0.3786710430610880E+00,
     -0.5468514142104170E+00,
     -0.7809468775455607E+00,
     -0.1134011957382327E+01,
     -0.1775800683423525E+01,
     -0.2443622553873225E+01,
     -0.3124190050507211E+01,
     -0.2872935510329120E+01,
     -0.2164282524138207E+01,
     -0.1440351296279408E+01,
     -0.6864884538258716E+00,
      0.1250649863152964E+00,
      0.1045163780117493E+01,
      0.2967585095039051E+01,
      0.5253718299558931E+01,
      0.8519716463711059E+01,
      0.1360509217709172E+02,
      0.2193466832805100E+02,
      0.3604254831722944E+02,
      0.6051306533791733E+02,
      0.1037211171690373E+03,
      0.1810780396816945E+03,
      0.3211144156746837E+03 };

  static double x_vec[N_MAX] = {
     0.000000E+00,
     0.100000E+00,
     0.200000E+00,
     0.300000E+00,
     0.400000E+00,
     0.500000E+00,
     0.600000E+00,
     0.700000E+00,
     0.800000E+00,
     0.900000E+00,
     0.950000E+00,
     0.975000E+00,
     0.103125E+01,
     0.106250E+01,
     0.112500E+01,
     0.125000E+01,
     0.150000E+01,
     0.200000E+01,
     0.400000E+01,
     0.800000E+01,
     0.160000E+02,
     0.320000E+02,
     0.640000E+02,
     0.128000E+03,
     0.256000E+03,
     0.512000E+03,
     0.102400E+04,
     0.204800E+04 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void logistic_cdf_values ( int &n_data, double &mu, double &beta, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    LOGISTIC_CDF_VALUES returns some values of the Logistic CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = LogisticDistribution [ mu, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &BETA, the shape parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double beta_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.8807970779778824E+00,
     0.9820137900379084E+00,
     0.9975273768433652E+00,
     0.6224593312018546E+00,
     0.5825702064623147E+00,
     0.5621765008857981E+00,
     0.5498339973124779E+00,
     0.6224593312018546E+00,
     0.5000000000000000E+00,
     0.3775406687981454E+00,
     0.2689414213699951E+00 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    beta = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    beta = beta_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void mertens_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    MERTENS_VALUES returns some values of the Mertens function.
//
//  Discussion:
//
//    The Mertens function M(N) is the sum from 1 to N of the Moebius
//    function MU.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 October 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    M Deleglise, J Rivat,
//    Computing the Summation of the Moebius Function,
//    Experimental Mathematics,
//    Volume 5, 1996, pages 291-295.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 2002,
//    Second edition,
//    ISBN: 1584883472,
//    LC: QA5.W45.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and N_DATA
//    is set to 1.  On each subsequent call, the input value of N_DATA is
//    incremented and that test data item is returned, if available.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the argument of the Mertens function.
//
//    Output, int &C, the value of the Mertens function.
//
{
# define N_MAX 15

  static int c_vec[N_MAX] = {
      1,   0,  -1,   -1,  -2,  -1,  -2,  -2,   -2,  -1,
     -2,  -2,   1,    2, -23 };
  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     11,  12,  100, 1000, 10000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  if ( N_MAX <= n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data];
    c = c_vec[n_data];
    n_data = n_data + 1;
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void moebius_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    MOEBIUS_VALUES returns some values of the Moebius function.
//
//  Discussion:
//
//    MU(N) is defined as follows:
//
//      MU(N) = 1 if N = 1;
//              0 if N is divisible by the square of a prime;
//              (-1)**K, if N is the product of K distinct primes.
//
//    In Mathematica, the function can be evaluated by:
//
//      MoebiusMu[n]
//
//  First values:
//
//     N  MU(N)
//
//     1    1
//     2   -1
//     3   -1
//     4    0
//     5   -1
//     6    1
//     7   -1
//     8    0
//     9    0
//    10    1
//    11   -1
//    12    0
//    13   -1
//    14    1
//    15    1
//    16    0
//    17   -1
//    18    0
//    19   -1
//    20    0
//
//  Note:
//
//    As special cases, MU(N) is -1 if N is a prime, and MU(N) is 0
//    if N is a square, cube, etc.
//
//  Formula:
//
//    The Moebius function is related to Euler's totient function:
//
//      PHI(N) = Sum ( D divides N ) MU(D) * ( N / D ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the Moebius function.
//
//    Output, int &C, the value of the Moebius function.
//
{
# define N_MAX 20

  static int c_vec[N_MAX] = {
      1,  -1,  -1,   0,  -1,   1,  -1,   0,   0,   1,
     -1,   0,  -1,   1,   1,   0,  -1,   0,  -1,   0 };

  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     11,  12,  13,  14,  15,  16,  17,  18,  19,  20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void negative_binomial_cdf_values ( int &n_data, int &f, int &s, double &p,
  double &cdf )

//****************************************************************************80
//
//  Purpose:
//
//    NEGATIVE_BINOMIAL_CDF_VALUES returns values of the negative binomial CDF.
//
//  Discussion:
//
//    Assume that a coin has a probability P of coming up heads on
//    any one trial.  Suppose that we plan to flip the coin until we
//    achieve a total of S heads.  If we let F represent the number of
//    tails that occur in this process, then the value of F satisfies
//    a negative binomial PDF:
//
//      PDF(F,S,P) = Choose ( F from F+S-1 ) * P**S * (1-P)**F
//
//    The negative binomial CDF is the probability that there are F or
//    fewer failures upon the attainment of the S-th success.  Thus,
//
//      CDF(F,S,P) = sum ( 0 <= G <= F ) PDF(G,S,P)
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = NegativeBinomialDistribution [ s, p ]
//      CDF [ dist, f ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    F C Powell,
//    Statistical Tables for Sociology, Biology and Physical Sciences,
//    Cambridge University Press, 1982.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &F, the maximum number of failures.
//
//    Output, int &S, the number of successes.
//
//    Output, double &P, the probability of a success on one trial.
//
//    Output, double &CDF, the probability of at most F failures
//    before the S-th success.
//
{
# define N_MAX 27

  static double cdf_vec[N_MAX] = {
     0.6367187500000000E+00,
     0.3632812500000000E+00,
     0.1445312500000000E+00,
     0.5000000000000000E+00,
     0.2265625000000000E+00,
     0.6250000000000000E-01,
     0.3437500000000000E+00,
     0.1093750000000000E+00,
     0.1562500000000000E-01,
     0.1792000000000000E+00,
     0.4096000000000000E-01,
     0.4096000000000000E-02,
     0.7047000000000000E-01,
     0.1093500000000000E-01,
     0.7290000000000000E-03,
     0.9861587127990000E+00,
     0.9149749500510000E+00,
     0.7471846521450000E+00,
     0.8499053647030009E+00,
     0.5497160941090026E+00,
     0.2662040052146710E+00,
     0.6513215599000000E+00,
     0.2639010709000000E+00,
     0.7019082640000000E-01,
     0.1000000000000000E+01,
     0.1990000000000000E-01,
     0.1000000000000000E-03 };

  static int f_vec[N_MAX] = {
     4,  3,  2,
     3,  2,  1,
     2,  1,  0,
     2,  1,  0,
     2,  1,  0,
    11, 10,  9,
    17, 16, 15,
     9,  8,  7,
     2,  1,  0 };

  static double p_vec[N_MAX] = {
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     0.40E+00,
     0.40E+00,
     0.40E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.30E+00,
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.10E+00,
     0.10E-01,
     0.10E-01,
     0.10E-01 };

  static int s_vec[N_MAX] = {
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    4, 5, 6,
    1, 2, 3,
    1, 2, 3,
    1, 2, 3,
    0, 1, 2 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    f = 0;
    s = 0;
    p = 0.0;
    cdf = 0.0;
  }
  else
  {
    f = f_vec[n_data-1];
    s = s_vec[n_data-1];
    p = p_vec[n_data-1];
    cdf = cdf_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void nine_j_values ( int &n_data, double &j1, double &j2, double &j3,
  double &j4, double &j5, double &j6, double &j7, double &j8, double &j9,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    NINE_J_VALUES returns some values of the Wigner 9J function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &J1, &J2, &J3, &J4, &J5, &J6, &J7, &J8, &J9,
//    the arguments of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 9

  static double fx_vec[N_MAX] = {
     0.0004270039294528318,
    -0.001228915451058514,
    -0.0001944260688400887,
     0.003338419923885592,
    -0.0007958936865080434,
    -0.004338208690251972,
     0.05379143536399187,
     0.006211299937499411,
     0.03042903097250921 };
  static double j1_vec[N_MAX] = {
    1.0,
    1.5,
    2.0,
    1.0,
    1.5,
    2.0,
    0.5,
    1.0,
    1.5 };
  static double j2_vec[N_MAX] = {
    8.0,
    8.0,
    8.0,
    3.0,
    3.0,
    3.0,
    0.5,
    0.5,
    0.5 };
  static double j3_vec[N_MAX] = {
    7.0,
    7.0,
    7.0,
    2.0,
    2.0,
    2.0,
    1.0,
    1.0,
    1.0 };
  static double j4_vec[N_MAX] = {
    6.5,
    6.5,
    6.5,
    4.0,
    4.0,
    4.0,
    2.0,
    2.0,
    2.0 };
  static double j5_vec[N_MAX] = {
    7.5,
    7.5,
    7.5,
    1.5,
    1.5,
    1.5,
    1.0,
    1.0,
    1.0 };
  static double j6_vec[N_MAX] = {
    7.5,
    7.5,
    7.5,
    3.0,
    3.0,
    3.0,
    1.5,
    1.5,
    1.5 };
  static double j7_vec[N_MAX] = {
    6.0,
    6.0,
    6.0,
    3.5,
    3.5,
    3.5,
    1.5,
    1.5,
    1.5 };
  static double j8_vec[N_MAX] = {
    10.0,
    10.0,
    10.0,
     2.0,
     2.0,
     2.0,
     0.5,
     0.5,
     0.5 };
  static double j9_vec[N_MAX] = {
    6.0,
    6.0,
    6.0,
    2.0,
    2.0,
    2.0,
    1.5,
    1.5,
    1.5 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    j1 = 0.0;
    j2 = 0.0;
    j3 = 0.0;
    j4 = 0.0;
    j5 = 0.0;
    j6 = 0.0;
    j7 = 0.0;
    j8 = 0.0;
    j9 = 0.0;
    fx = 0.0;
  }
  else
  {
    j1 = j1_vec[n_data-1];
    j2 = j2_vec[n_data-1];
    j3 = j3_vec[n_data-1];
    j4 = j4_vec[n_data-1];
    j5 = j5_vec[n_data-1];
    j6 = j6_vec[n_data-1];
    j7 = j7_vec[n_data-1];
    j8 = j8_vec[n_data-1];
    j9 = j9_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void normal_cdf_values ( int &n_data, double &mu, double &sigma, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_CDF_VALUES returns some values of the Normal CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NormalDistribution [ mu, sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &MU, the mean of the distribution.
//
//    Output, double &SIGMA, the variance of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.9772498680518208E+00,
     0.9999683287581669E+00,
     0.9999999990134124E+00,
     0.6914624612740131E+00,
     0.6305586598182364E+00,
     0.5987063256829237E+00,
     0.5792597094391030E+00,
     0.6914624612740131E+00,
     0.5000000000000000E+00,
     0.3085375387259869E+00,
     0.1586552539314571E+00 };

  static double mu_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    mu = 0.0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    mu = mu_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void normal_01_cdf_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NormalDistribution [ 0, 1 ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 17

  static double fx_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5398278372770290E+00,
     0.5792597094391030E+00,
     0.6179114221889526E+00,
     0.6554217416103242E+00,
     0.6914624612740131E+00,
     0.7257468822499270E+00,
     0.7580363477769270E+00,
     0.7881446014166033E+00,
     0.8159398746532405E+00,
     0.8413447460685429E+00,
     0.9331927987311419E+00,
     0.9772498680518208E+00,
     0.9937903346742239E+00,
     0.9986501019683699E+00,
     0.9997673709209645E+00,
     0.9999683287581669E+00 };

  static double x_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.1000000000000000E+00,
     0.2000000000000000E+00,
     0.3000000000000000E+00,
     0.4000000000000000E+00,
     0.5000000000000000E+00,
     0.6000000000000000E+00,
     0.7000000000000000E+00,
     0.8000000000000000E+00,
     0.9000000000000000E+00,
     0.1000000000000000E+01,
     0.1500000000000000E+01,
     0.2000000000000000E+01,
     0.2500000000000000E+01,
     0.3000000000000000E+01,
     0.3500000000000000E+01,
     0.4000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void omega_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    OMEGA_VALUES returns some values of the OMEGA function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by
//
//      Length [ FactorInteger [ n ] ]
//
//  First values:
//
//     N   OMEGA(N)
//
//     1    1
//     2    1
//     3    1
//     4    1
//     5    1
//     6    2
//     7    1
//     8    1
//     9    1
//    10    2
//    11    1
//    12    2
//    13    1
//    14    2
//    15    2
//    16    1
//    17    1
//    18    2
//    19    1
//    20    2
//
//  Formula:
//
//    If N = 1, then
//
//      OMEGA(N) = 1
//
//    else if the prime factorization of N is
//
//      N = P1**E1 * P2**E2 * ... * PM**EM,
//
//    then
//
//      OMEGA(N) = M
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the OMEGA function.
//
//    Output, int &C, the value of the OMEGA function.
//
{
# define N_MAX 23

  static int c_vec[N_MAX] = {
      1,   1,   1,   1,   1,
      2,   1,   1,   1,   2,
      3,   1,   4,   4,   3,
      1,   5,   2,   2,   1,
      6,   7,   8 };

  static int n_vec[N_MAX] = {
           1,
           2,
           3,
           4,
           5,
           6,
           7,
           8,
           9,
          10,
          30,
         101,
         210,
        1320,
        1764,
        2003,
        2310,
        2827,
        8717,
       12553,
       30030,
      510510,
     9699690 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void owen_values ( int &n_data, double &h, double &a, double &t )

//****************************************************************************80
//
//  Purpose:
//
//    OWEN_VALUES returns some values of Owen's T function.
//
//  Discussion:
//
//    Owen's T function is useful for computation of the bivariate normal
//    distribution and the distribution of a skewed normal distribution.
//
//    Although it was originally formulated in terms of the bivariate
//    normal function, the function can be defined more directly as
//
//      T(H,A) = 1 / ( 2 * pi ) *
//        Integral ( 0 <= X <= A ) e^(H^2*(1+X^2)/2) / (1+X^2) dX
//
//    In Mathematica, the function can be evaluated by:
//
//      fx = 1/(2*Pi) * Integrate [ E^(-h^2*(1+x^2)/2)/(1+x^2), {x,0,a} ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 May 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Mike Patefield, David Tandy,
//    Fast and Accurate Calculation of Owen's T Function,
//    Journal of Statistical Software,
//    Volume 5, Number 5, 2000, pages 1-25.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &H, a parameter.
//
//    Output, double &A, the upper limit of the integral.
//
//    Output, double &T, the value of the function.
//
{
# define N_MAX 28

  static double a_vec[N_MAX] = {
    0.2500000000000000E+00,
    0.4375000000000000E+00,
    0.9687500000000000E+00,
    0.0625000000000000E+00,
    0.5000000000000000E+00,
    0.9999975000000000E+00,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.5000000000000000E+00,
    0.1000000000000000E+01,
    0.2000000000000000E+01,
    0.3000000000000000E+01,
    0.1000000000000000E+02,
    0.1000000000000000E+03 };

  static double h_vec[N_MAX] = {
    0.0625000000000000E+00,
    6.5000000000000000E+00,
    7.0000000000000000E+00,
    4.7812500000000000E+00,
    2.0000000000000000E+00,
    1.0000000000000000E+00,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.2500000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.1250000000000000E+00,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02,
    0.7812500000000000E-02 };

  static double t_vec[N_MAX] = {
    3.8911930234701366E-02,
    2.0005773048508315E-11,
    6.3990627193898685E-13,
    1.0632974804687463E-07,
    8.6250779855215071E-03,
    6.6741808978228592E-02,
    0.4306469112078537E-01,
    0.6674188216570097E-01,
    0.7846818699308410E-01,
    0.7929950474887259E-01,
    0.6448860284750376E-01,
    0.1066710629614485E+00,
    0.1415806036539784E+00,
    0.1510840430760184E+00,
    0.7134663382271778E-01,
    0.1201285306350883E+00,
    0.1666128410939293E+00,
    0.1847501847929859E+00,
    0.7317273327500385E-01,
    0.1237630544953746E+00,
    0.1737438887583106E+00,
    0.1951190307092811E+00,
    0.7378938035365546E-01,
    0.1249951430754052E+00,
    0.1761984774738108E+00,
    0.1987772386442824E+00,
    0.2340886964802671E+00,
    0.2479460829231492E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    h = 0.0;
    a = 0.0;
    t = 0.0;
  }
  else
  {
    h = h_vec[n_data-1];
    a = a_vec[n_data-1];
    t = t_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void partition_count_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_COUNT_VALUES returns some values of the int *partition count.
//
//  Discussion:
//
//    A partition of an int *N is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  The number of partitions of N is symbolized
//    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
//    following partitions:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//      = 3 + 1 + 1
//      = 2 + 2 + 1
//      = 2 + 1 + 1 + 1
//      = 1 + 1 + 1 + 1 + 1
//
//    In Mathematica, the function can be evaluated by
//
//      PartitionsP[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the integer.
//
//    Output, int &C, the number of partitions of the integer.
//
{
# define N_MAX 21

  static int c_vec[N_MAX] = {
      1,
      1,   2,   3,   5,   7,  11,  15,  22,  30,  42,
     56,  77, 101, 135, 176, 231, 297, 385, 490, 627 };

  static int n_vec[N_MAX] = {
     0,
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void partition_distinct_count_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    PARTITION_DISTINCT_COUNT_VALUES returns some values of Q(N).
//
//  Discussion:
//
//    A partition of an int *N is a representation of the integer
//    as the sum of nonzero positive integers.  The order of the summands
//    does not matter.  The number of partitions of N is symbolized
//    by P(N).  Thus, the number 5 has P(N) = 7, because it has the
//    following partitions:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//      = 3 + 1 + 1
//      = 2 + 2 + 1
//      = 2 + 1 + 1 + 1
//      = 1 + 1 + 1 + 1 + 1
//
//    However, if we require that each member of the partition
//    be distinct, so that no nonzero summand occurs more than once,
//    we are computing something symbolized by Q(N).
//    The number 5 has Q(N) = 3, because it has the following partitions
//    into distinct parts:
//
//    5 = 5
//      = 4 + 1
//      = 3 + 2
//
//    In Mathematica, the function can be evaluated by
//
//      PartitionsQ[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the integer.
//
//    Output, int &C, the number of partitions of the integer
//    into distinct parts.
//
{
# define N_MAX 21

  static int c_vec[N_MAX] = {
      1,
      1,   1,   2,   2,   3,   4,   5,   6,   8,  10,
     12,  15,  18,  22,  27,  32,  38,  46,  54,  64 };

  static int n_vec[N_MAX] = {
     0,
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void phi_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    PHI_VALUES returns some values of the PHI function.
//
//  Discussion:
//
//    PHI(N) is the number of integers between 1 and N which are
//    relatively prime to N.  I and J are relatively prime if they
//    have no common factors.  The function PHI(N) is known as
//    "Euler's totient function".
//
//    By convention, 1 and N are relatively prime.
//
//    In Mathematica, the function can be evaluated by:
//
//      EulerPhi[n]
//
//  First values:
//
//     N  PHI(N)
//
//     1    1
//     2    1
//     3    2
//     4    2
//     5    4
//     6    2
//     7    6
//     8    4
//     9    6
//    10    4
//    11   10
//    12    4
//    13   12
//    14    6
//    15    8
//    16    8
//    17   16
//    18    6
//    19   18
//    20    8
//
//  Formula:
//
//    PHI(U*V) = PHI(U) * PHI(V) if U and V are relatively prime.
//
//    PHI(P**K) = P**(K-1) * ( P - 1 ) if P is prime.
//
//    PHI(N) = N * Product ( P divides N ) ( 1 - 1 / P )
//
//    N = Sum ( D divides N ) PHI(D).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the PHI function.
//
//    Output, int &C, the value of the PHI function.
//
{
# define N_MAX 20

  static int c_vec[N_MAX] = {
      1,   1,   2,   2,   4,   2,   6,   4,   6,   4,
      8,   8,  16,  20,  16,  40, 148, 200, 200, 648 };

  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     20,  30,  40,  50,  60, 100, 149, 500, 750, 999 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void pi_values ( int &n_data, int &n, int &p )

//****************************************************************************80
//
//  Purpose:
//
//    PI_VALUES returns values of the Pi function.
//
//  Discussion:
//
//    Pi[n] is the number of primes less than or equal to n.
//
//    In Mathematica, the function can be evaluated by:
//
//      PrimePi[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument.
//
//    Output, int &P, the value of the function.
//
{
# define N_MAX 17

  static int n_vec[N_MAX] = {
            10,
            20,
            30,
            40,
            50,
            60,
            70,
            80,
            90,
           100,
          1000,
         10000,
        100000,
       1000000,
      10000000,
     100000000,
    1000000000 };

  static int p_vec[N_MAX] = {
             4,
             8,
            10,
            12,
            15,
            17,
            19,
            22,
            24,
            25,
           168,
          1229,
          9592,
         78498,
        664579,
       5761455,
      50847534 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    p = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    p = p_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void pochhammer_values ( int &n_data, double &x, double &y, double &fxy )

//****************************************************************************80
//
//  Purpose:
//
//    POCHHAMMER_VALUES returns some values of the Pochhammer function.
//
//  Discussion:
//
//    Pochhammer(X,Y) = Gamma(X+Y) / Gamma(X)
//
//    For integer arguments, Pochhammer(M,N) = ( M + N - 1 )! / ( N - 1 )!
//
//    In Mathematica, the function can be evaluated by:
//
//      Pochhammer[X,Y]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 April 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, &Y, the arguments of the function.
//
//    Output, double &FXY, the value of the function.
//
{
# define N_MAX 19

  static double f_vec[N_MAX] = {
     720.0000000000000,
       1.875000000000000,
       1.000000000000000,
       4.500000000000000,
      24.75000000000000,
     110.0000000000000,
     377.1036305819165,
       4.362197352456253,
       1.000000000000000,
       1.467150493866654,
       2.180949074356397,
       3.282686710888467,
       5.000000000000000,
       7.702540092799931,
      11.99521990896018,
      18.87544858760869,
      30.00000000000000,
      48.14087557999957,
      77.96892940824118 };

  static double x_vec[N_MAX] = {
     1.00,
     0.50,
     4.50,
     4.50,
     4.50,
    10.00,
    10.00,
     7.25,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00,
     5.00 };

  static double y_vec[N_MAX] = {
    6.00,
    3.00,
    0.00,
    1.00,
    2.00,
    2.00,
    2.50,
    0.75,
    0.00,
    0.25,
    0.50,
    0.75,
    1.00,
    1.25,
    1.50,
    1.75,
    2.00,
    2.25,
    2.50 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    y = 0.0;
    fxy = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    y = y_vec[n_data-1];
    fxy = f_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void poisson_cdf_values ( int &n_data, double &a, int &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    POISSON_CDF_VALUES returns some values of the Poisson CDF.
//
//  Discussion:
//
//    CDF(X)(A) is the probability of at most X successes in unit time,
//    given that the expected mean number of successes is A.
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`DiscreteDistributions`]
//      dist = PoissonDistribution [ a ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition, CRC Press, 1996, pages 653-658.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, the parameter of the function.
//
//    Output, int *X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 21

  static double a_vec[N_MAX] = {
     0.02E+00,
     0.10E+00,
     0.10E+00,
     0.50E+00,
     0.50E+00,
     0.50E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     1.00E+00,
     2.00E+00,
     2.00E+00,
     2.00E+00,
     2.00E+00,
     5.00E+00,
     5.00E+00,
     5.00E+00,
     5.00E+00,
     5.00E+00,
     5.00E+00,
     5.00E+00 };

  static double fx_vec[N_MAX] = {
     0.9801986733067553E+00,
     0.9048374180359596E+00,
     0.9953211598395555E+00,
     0.6065306597126334E+00,
     0.9097959895689501E+00,
     0.9856123220330293E+00,
     0.3678794411714423E+00,
     0.7357588823428846E+00,
     0.9196986029286058E+00,
     0.9810118431238462E+00,
     0.1353352832366127E+00,
     0.4060058497098381E+00,
     0.6766764161830635E+00,
     0.8571234604985470E+00,
     0.6737946999085467E-02,
     0.4042768199451280E-01,
     0.1246520194830811E+00,
     0.2650259152973617E+00,
     0.4404932850652124E+00,
     0.6159606548330631E+00,
     0.7621834629729387E+00 };

  static int x_vec[N_MAX] = {
     0, 0, 1, 0,
     1, 2, 0, 1,
     2, 3, 0, 1,
     2, 3, 0, 1,
     2, 3, 4, 5,
     6 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    x = 0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void polylogarithm_values ( int &n_data, int &n, double &z, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    POLYLOGARITHM_VALUES returns some values of the polylogarithm.
//
//  Discussion:
//
//    The polylogarithm of n and z is defined as
//
//      f[n,z] = Sum ( 1 <= k < +oo ) z^k / k^n
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyLog[n,z]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    14 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the exponent of the denominator.
//
//    Output, double &Z, the base of the numerator.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double fx_vec[N_MAX] = {
    0.1644934066848226E+01,
    0.1202056903159594E+01,
    0.1000994575127818E+01,
    0.5822405264650125E+00,
    0.5372131936080402E+00,
    0.5002463206060068E+00,
    0.3662132299770635E+00,
    0.3488278611548401E+00,
    0.3334424797228716E+00,
    0.1026177910993911E+00,
    0.1012886844792230E+00,
    0.1000097826564961E+00 };

  static int n_vec[N_MAX] = {
     2, 3, 10, 2, 3, 10, 2, 3, 10, 2, 3, 10 };

  static double z_vec[N_MAX] = {
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.1000000000000000E+01,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.5000000000000000E+00,
    0.3333333333333333E+00,
    0.3333333333333333E+00,
    0.3333333333333333E+00,
    0.1000000000000000E+00,
    0.1000000000000000E+00,
    0.1000000000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    z = 0.0;
    fx = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    z = z_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void prandtl_values ( int &n_data, double &tc, double &p, double &pr )

//****************************************************************************80
//
//  Purpose:
//
//    PRANDTL_VALUES returns some values of the Prandtl number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, page 265.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &PR, the Prandtl number, dimensionless.
//
{
# define N_MAX 35

  static double pr_vec[N_MAX] = {
     13.50E+00,
     13.48E+00,
     13.46E+00,
     13.39E+00,
     13.27E+00,
     13.15E+00,
     13.04E+00,
     12.93E+00,
     12.83E+00,
     12.73E+00,
     12.63E+00,
     12.53E+00,
     12.43E+00,
     12.34E+00,
     12.25E+00,
     12.08E+00,
     11.92E+00,
     11.77E+00,
     11.62E+00,
     11.48E+00,
     11.36E+00,
     11.23E+00,
     11.12E+00,
     10.91E+00,
     10.72E+00,
     10.55E+00,
      6.137E+00,
      3.555E+00,
      2.378E+00,
      1.000E+00,
      0.974E+00,
      0.960E+00,
      0.924E+00,
      0.899E+00,
      0.882E+00 };

  static double p_vec[N_MAX] = {
        1.0E+00,
        5.0E+00,
       10.0E+00,
       25.0E+00,
       50.0E+00,
       75.0E+00,
      100.0E+00,
      125.0E+00,
      150.0E+00,
      175.0E+00,
      200.0E+00,
      225.0E+00,
      250.0E+00,
      275.0E+00,
      300.0E+00,
      350.0E+00,
      400.0E+00,
      450.0E+00,
      500.0E+00,
      550.0E+00,
      600.0E+00,
      650.0E+00,
      700.0E+00,
      800.0E+00,
      900.0E+00,
     1000.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00 };

  static double tc_vec[N_MAX] = {
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
      25.0E+00,
      50.0E+00,
      75.0E+00,
     100.0E+00,
     150.0E+00,
     200.0E+00,
     400.0E+00,
     600.0E+00,
     800.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    pr = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    pr = pr_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void prime_values ( int &n_data, int &n, int &p )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME_VALUES returns values of the prime function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Prime[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the index of the prime.
//
//    Output, int &P, the value of the prime.
//
{
# define N_MAX 24

  static int n_vec[N_MAX] = {
          1,
          2,
          4,
          8,
         16,
         32,
         64,
        128,
        256,
        512,
       1000,
       2000,
       4000,
       8000,
      16000,
      32000,
      64000,
     128000,
     256000,
     512000,
    1024000,
    2048000,
    4096000,
    8129000 };

  static int p_vec[N_MAX] = {
            2,
            3,
            7,
           19,
           53,
          131,
          311,
          719,
         1619,
         3671,
         7919,
        17389,
        37813,
        81799,
       176081,
       376127,
       800573,
      1698077,
      3588941,
      7559173,
     15881419,
     33283031,
     69600977,
    145253029 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    p = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    p = p_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void psat_values ( int &n_data, double &tc, double &p )

//****************************************************************************80
//
//  Purpose:
//
//    PSAT_VALUES returns some values of the saturation pressure.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, pages 9-15.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the saturation pressure, in bar.
//
{
# define N_MAX 12

  static double p_vec[N_MAX] = {
     0.0061173E+00,
     0.0065716E+00,
     0.0087260E+00,
     0.12344E+00,
     1.0132E+00,
     2.3201E+00,
     4.7572E+00,
     15.537E+00,
     39.737E+00,
     85.838E+00,
     165.21E+00,
     220.55E+00 };

  static double tc_vec[N_MAX] = {
     0.100000E-01,
     0.100000E+01,
     0.500000E+01,
     0.500000E+02,
     0.100000E+03,
     0.125000E+03,
     0.150000E+03,
     0.200000E+03,
     0.250000E+03,
     0.300000E+03,
     0.350000E+03,
     0.373976E+03  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void psi_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    PSI_VALUES returns some values of the Psi or Digamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyGamma[x]
//
//    or
//
//      Polygamma[0,x]
//
//    PSI(X) = d ln ( Gamma ( X ) ) / d X = Gamma'(X) / Gamma(X)
//
//    PSI(1) = -Euler's constant.
//
//    PSI(X+1) = PSI(X) + 1 / X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
     -0.5772156649015329E+00,
     -0.4237549404110768E+00,
     -0.2890398965921883E+00,
     -0.1691908888667997E+00,
     -0.6138454458511615E-01,
      0.3648997397857652E-01,
      0.1260474527734763E+00,
      0.2085478748734940E+00,
      0.2849914332938615E+00,
      0.3561841611640597E+00,
      0.4227843350984671E+00 };

  static double x_vec[N_MAX] = {
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8_factorial_values ( int &n_data, int &n, double &fn )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_VALUES returns values of the real factorial function.
//
//  Discussion:
//
//    0! = 1
//    I! = Product ( 1 <= J <= I ) J
//
//    Although the factorial is an int *valued function, it quickly
//    becomes too large for an int *to hold.  This routine still accepts
//    an int *as the input argument, but returns the function value
//    as a real number.
//
//    In Mathematica, the function can be evaluated by:
//
//      n!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the function.
//
//    Output, double &FN, the value of the function.
//
{
# define N_MAX 25

  static double fn_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.6000000000000000E+01,
     0.2400000000000000E+02,
     0.1200000000000000E+03,
     0.7200000000000000E+03,
     0.5040000000000000E+04,
     0.4032000000000000E+05,
     0.3628800000000000E+06,
     0.3628800000000000E+07,
     0.3991680000000000E+08,
     0.4790016000000000E+09,
     0.6227020800000000E+10,
     0.8717829120000000E+11,
     0.1307674368000000E+13,
     0.2092278988800000E+14,
     0.3556874280960000E+15,
     0.6402373705728000E+16,
     0.1216451004088320E+18,
     0.2432902008176640E+19,
     0.1551121004333099E+26,
     0.3041409320171338E+65,
     0.9332621544394415E+158,
     0.5713383956445855E+263 };

  static int n_vec[N_MAX] = {
       0,
       1,
       2,
       3,
       4,
       5,
       6,
       7,
       8,
       9,
      10,
      11,
      12,
      13,
      14,
      15,
      16,
      17,
      18,
      19,
      20,
      25,
      50,
     100,
     150 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void r8_factorial_log_values ( int &n_data, int &n, double &fn )

//****************************************************************************80
//
//  Purpose:
//
//    R8_FACTORIAL_LOG_VALUES returns values of log(n!).
//
//  Discussion:
//
//    The function log(n!) can be written as
//
//     log(n!) = sum ( 1 <= i <= n ) log ( i )
//
//    In Mathematica, the function can be evaluated by:
//
//      Log[n!]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//    Daniel Zwillinger, editor,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the function.
//
//    Output, double &FN, the value of the function.
//
{
# define N_MAX 27

  static double fn_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.0000000000000000E+00,
     0.6931471805599453E+00,
     0.1791759469228055E+01,
     0.3178053830347946E+01,
     0.4787491742782046E+01,
     0.6579251212010101E+01,
     0.8525161361065414E+01,
     0.1060460290274525E+02,
     0.1280182748008147E+02,
     0.1510441257307552E+02,
     0.1750230784587389E+02,
     0.1998721449566189E+02,
     0.2255216385312342E+02,
     0.2519122118273868E+02,
     0.2789927138384089E+02,
     0.3067186010608067E+02,
     0.3350507345013689E+02,
     0.3639544520803305E+02,
     0.3933988418719949E+02,
     0.4233561646075349E+02,
     0.5800360522298052E+02,
     0.1484777669517730E+03,
     0.3637393755555635E+03,
     0.6050201058494237E+03,
     0.2611330458460156E+04,
     0.5912128178488163E+04 };

  static int n_vec[N_MAX] = {
       0,
       1,
       2,
       3,
       4,
       5,
       6,
       7,
       8,
       9,
      10,
      11,
      12,
      13,
      14,
      15,
      16,
      17,
      18,
      19,
      20,
      25,
      50,
     100,
     150,
     500,
    1000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void rayleigh_cdf_values ( int &n_data, double &sigma, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    RAYLEIGH_CDF_VALUES returns some values of the Rayleigh CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = RayleighDistribution [ sigma ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &SIGMA, the shape parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 9

  static double fx_vec[N_MAX] = {
     0.8646647167633873E+00,
     0.9996645373720975E+00,
     0.9999999847700203E+00,
     0.999999999999987E+00,
     0.8646647167633873E+00,
     0.3934693402873666E+00,
     0.1992625970831920E+00,
     0.1175030974154046E+00,
     0.7688365361336422E-01 };

  static double sigma_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    sigma = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    sigma = sigma_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void secvir_values ( int &n_data, double &tc, double &vir )

//****************************************************************************80
//
//  Purpose:
//
//    SECVIR_VALUES returns some values of the second virial coefficient.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, pages 24-25.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &VIR, the second virial coefficient, in
//    m^3/kg.
//
{
# define N_MAX 19

  static double tc_vec[N_MAX] = {
        0.0E+00,
        5.0E+00,
       10.0E+00,
       20.0E+00,
       30.0E+00,
       40.0E+00,
       60.0E+00,
       90.0E+00,
      120.0E+00,
      150.0E+00,
      180.0E+00,
      210.0E+00,
      240.0E+00,
      300.0E+00,
      400.0E+00,
      500.0E+00,
      700.0E+00,
     1000.0E+00,
     2000.0E+00 };

  static double vir_vec[N_MAX] = {
     -98.96E+00,
     -90.08E+00,
     -82.29E+00,
     -69.36E+00,
     -59.19E+00,
     -51.07E+00,
     -39.13E+00,
     -27.81E+00,
     -20.83E+00,
     -16.21E+00,
     -12.98E+00,
     -10.63E+00,
      -8.85E+00,
      -6.39E+00,
      -4.03E+00,
      -2.71E+00,
      -1.32E+00,
      -0.39E+00,
       0.53E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    vir = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    vir = vir_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void shi_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SHI_VALUES returns some values of the hyperbolic sine integral function.
//
//  Discussion:
//
//    SHI(X) = integral ( 0 <= T <= X ) sinh ( T ) / T dt
//
//    In Mathematica, the function can be evaluated by:
//
//      SinhIntegral[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
    0.5069967498196672,
    0.6121303965633808,
    0.7193380189288998,
    0.8289965633789345,
    0.9414978265114335,
    1.057250875375729,
    1.300250361022057,
    1.561713388361002,
    1.845814141358504,
    2.157290343425901,
    2.501567433354976,
    3.549340406224435,
    4.973440475859807,
    6.966162067504942,
    9.817326911233034,
    13.96788504934715 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void si_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SI_VALUES returns some values of the sine integral function.
//
//  Discussion:
//
//    SI(X) = integral ( 0 <= T <= X ) sin ( T ) / T dt
//
//    In Mathematica, the function can be evaluated by:
//
//      SinIntegral[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 16

  static double fx_vec[N_MAX] = {
     0.4931074180430667E+00,
     0.5881288096080801E+00,
     0.6812222391166113E+00,
     0.7720957854819966E+00,
     0.8604707107452929E+00,
     0.9460830703671830E+00,
     0.1108047199013719E+01,
     0.1256226732779218E+01,
     0.1389180485870438E+01,
     0.1505816780255579E+01,
     0.1605412976802695E+01,
     0.1778520173443827E+01,
     0.1848652527999468E+01,
     0.1833125398665997E+01,
     0.1758203138949053E+01,
     0.1654140414379244E+01 };

  static double x_vec[N_MAX] = {
      0.5E+00,
      0.6E+00,
      0.7E+00,
      0.8E+00,
      0.9E+00,
      1.0E+00,
      1.2E+00,
      1.4E+00,
      1.6E+00,
      1.8E+00,
      2.0E+00,
      2.5E+00,
      3.0E+00,
      3.5E+00,
      4.0E+00,
      4.5E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sigma_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA_VALUES returns some values of the Sigma function.
//
//  Discussion:
//
//    SIGMA(N) is the sum of the distinct divisors of N, including 1 and N.
//
//    In Mathematica, the function can be evaluated by:
//
//      DivisorSigma[1,n]
//
//  First values:
//
//     N  SIGMA(N)
//
//     1    1
//     2    3
//     3    4
//     4    7
//     5    6
//     6   12
//     7    8
//     8   15
//     9   13
//    10   18
//    11   12
//    12   28
//    13   14
//    14   24
//    15   24
//    16   31
//    17   18
//    18   39
//    19   20
//    20   42
//
//  Formula:
//
//    SIGMA(U*V) = SIGMA(U) * SIGMA(V) if U and V are relatively prime.
//
//    SIGMA(P**K) = ( P**(K+1) - 1 ) / ( P - 1 ) if P is prime.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the Sigma function.
//
//    Output, int &C, the value of the Sigma function.
//
{
# define N_MAX 20

  static int c_vec[N_MAX] = {
     1,    3,    4,    7,    6,   12,    8,   15,   13,   18,
    72,  128,  255,  176,  576, 1170,  618,  984, 2232, 2340 };

  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,   10,
     30, 127, 128, 129, 210, 360, 617, 815, 816, 1000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sin_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_VALUES returns some values of the sine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sin[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double fx_vec[N_MAX] = {
     0.00000000000000000000,
     0.25881904510252076235,
     0.47942553860420300027,
     0.50000000000000000000,
     0.70710678118654752440,
     0.84147098480789650665,
     0.86602540378443864676,
     1.00000000000000000000,
     0.90929742682568169540,
     0.14112000805986722210,
     0.00000000000000000000,
    -0.75680249530792825137,
    -0.95892427466313846889 };

  static double x_vec[N_MAX] = {
    0.0000000000000000000,
    0.26179938779914943654,
    0.50000000000000000000,
    0.52359877559829887308,
    0.78539816339744830962,
    1.0000000000000000000,
    1.0471975511965977462,
    1.5707963267948966192,
    2.0000000000000000000,
    3.0000000000000000000,
    3.1415926535897932385,
    4.0000000000000000000,
    5.0000000000000000000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sin_degree_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_DEGREE_VALUES: the sine function with argument in degrees.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sin[x Degree]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 March 2010
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 22

  static double fx_vec[N_MAX] = {
    -0.087155742747658173558,
     0.000000000000000000000,
     0.017452406437283512819,
     0.034899496702500971646,
     0.052335956242943832722,
     0.069756473744125300776,
     0.087155742747658173558,
     0.17364817766693034885,
     0.25881904510252076235,
     0.50000000000000000000,
     0.70710678118654752440,
     0.86602540378443864676,
     0.96592582628906828675,
     0.99619469809174553230,
     0.99756405025982424761,
     0.99862953475457387378,
     0.99939082701909573001,
     0.99984769515639123916,
     1.0000000000000000000,
     0.99984769515639123916,
     0.96592582628906828675,
     0.00000000000000000000 };
  static double x_vec[N_MAX] = {
     -5.0,
      0.0,
      1.0,
      2.0,
      3.0,
      4.0,
      5.0,
     10.0,
     15.0,
     30.0,
     45.0,
     60.0,
     75.0,
     85.0,
     86.0,
     87.0,
     88.0,
     89.0,
     90.0,
     91.0,
    105.0,
    180.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sin_power_int_values ( int &n_data, double &a, double &b, int &n,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIN_POWER_INT_VALUES returns some values of the sine power integral.
//
//  Discussion:
//
//    The function has the form
//
//      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
//
//    In Mathematica, the function can be evaluated by:
//
//      Integrate [ ( Sin[x] )^n, { x, a, b } ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the limits of integration.
//
//    Output, int &N, the power.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 10

  static double a_vec[N_MAX] = {
      0.10E+02,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.00E+00,
      0.10E+01,
      0.00E+00,
      0.00E+00 };

  static double b_vec[N_MAX] = {
      0.20E+02,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.10E+01,
      0.20E+01,
      0.20E+01,
      0.10E+01,
      0.10E+01 };

  static double fx_vec[N_MAX] = {
     0.10000000000000000000E+02,
     0.45969769413186028260E+00,
     0.27267564329357957615E+00,
     0.17894056254885809051E+00,
     0.12402556531520681830E+00,
     0.88974396451575946519E-01,
     0.90393123848149944133E+00,
     0.81495684202992349481E+00,
     0.21887522421729849008E-01,
     0.17023439374069324596E-01 };

  static int n_vec[N_MAX] = {
     0,
     1,
     2,
     3,
     4,
     5,
     5,
     5,
    10,
    11 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    n = 0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    n = n_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sinh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SINH_VALUES returns some values of the hyperbolic sine function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Sinh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
      -74.203210577788758977,
       -1.1752011936438014569,
        0.00000000000000000000,
        0.10016675001984402582,
        0.20133600254109398763,
        0.30452029344714261896,
        0.41075232580281550854,
        0.52109530549374736162,
        0.63665358214824127112,
        0.75858370183953350346,
        0.88810598218762300657,
        1.0265167257081752760,
        1.1752011936438014569,
        3.6268604078470187677,
       10.017874927409901899,
       27.289917197127752449,
       74.203210577788758977,
    11013.232874703393377 };

  static double x_vec[N_MAX] = {
   -5.0,
   -1.0,
    0.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
   10.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void six_j_values ( int &n_data, double &j1, double &j2, double &j3,
  double &j4, double &j5, double &j6, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SIX_J_VALUES returns some values of the Wigner 6J function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      SixJSymbol[{j1,j2,j3},{j4,j5,j6}]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &J1, &J2, &J3, &J4, &J5, &J6, the arguments
//    of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double fx_vec[N_MAX] = {
     0.03490905138373300,
    -0.03743025039659792,
     0.01890866390959560,
     0.007342448254928643,
    -0.02358935185081794,
     0.01913476955215437,
     0.001288017397724172,
    -0.01930018366290527,
     0.01677305949382889,
     0.005501147274850949,
    -0.02135439790896831,
     0.003460364451435387,
     0.02520950054795585,
     0.01483990561221713,
     0.002708577680633186 };
  static double j1_vec[N_MAX] = {
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
   10.0,
   11.0,
   12.0,
   13.0,
   14.0,
   15.0 };
  static double j2_vec[N_MAX] = {
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0,
    8.0 };
  static double j3_vec[N_MAX] = {
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0,
    7.0 };
  static double j4_vec[N_MAX] = {
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5,
    6.5 };
  static double j5_vec[N_MAX] = {
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5 };
  static double j6_vec[N_MAX] = {
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5,
    7.5 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    j1 = 0.0;
    j2 = 0.0;
    j3 = 0.0;
    j4 = 0.0;
    j5 = 0.0;
    j6 = 0.0;
    fx = 0.0;
  }
  else
  {
    j1 = j1_vec[n_data-1];
    j2 = j2_vec[n_data-1];
    j3 = j3_vec[n_data-1];
    j4 = j4_vec[n_data-1];
    j5 = j5_vec[n_data-1];
    j6 = j6_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sound_values ( int &n_data, double &tc, double &p, double &c )

//****************************************************************************80
//
//  Purpose:
//
//    SOUND_VALUES returns some values of the speed of sound.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, page 238-246.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &C, the speed of sound, in m/s.
//
{
# define N_MAX 20

  static double c_vec[N_MAX] = {
     1401.0E+00,
      472.8E+00,
      533.7E+00,
      585.7E+00,
      609.5E+00,
      632.2E+00,
      674.6E+00,
      713.9E+00,
      802.0E+00,
      880.1E+00,
     1017.8E+00,
     1115.9E+00,
     1401.7E+00,
     1402.6E+00,
     1409.6E+00,
     1418.1E+00,
     1443.1E+00,
     1484.6E+00,
     1577.1E+00,
     1913.4E+00 };

  static double p_vec[N_MAX] = {
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        5.0E+00,
       10.0E+00,
       50.0E+00,
      100.0E+00,
      250.0E+00,
      500.0E+00,
     1000.0E+00,
     2500.0E+00 };

  static double tc_vec[N_MAX] = {
        0.0E+00,
      100.0E+00,
      200.0E+00,
      300.0E+00,
      350.0E+00,
      400.0E+00,
      500.0E+00,
      600.0E+00,
      850.0E+00,
     1100.0E+00,
     1600.0E+00,
     2000.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00,
        0.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    c = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sphere_unit_area_values ( int &n_data, int &n, double &area )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_AREA_VALUES returns some areas of the unit sphere in ND.
//
//  Discussion:
//
//    The formula for the surface area of the unit sphere in N dimensions is:
//
//      Sphere_Unit_Area ( N ) = 2 * PI^(N/2) / Gamma ( N / 2 )
//
//    Some values of the function include:
//
//       N   Area
//
//       2    2        * PI
//       3  ( 4 /    ) * PI
//       4  ( 2 /   1) * PI^2
//       5  ( 8 /   3) * PI^2
//       6  ( 1 /   1) * PI^3
//       7  (16 /  15) * PI^3
//       8  ( 1 /   3) * PI^4
//       9  (32 / 105) * PI^4
//      10  ( 1 /  12) * PI^5
//
//    For the unit sphere, Area(N) = N * Volume(N)
//
//    In Mathematica, the function can be evaluated by:
//
//      2 * Pi^(n/2) / Gamma[n/2]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the spatial dimension.
//
//    Output, double &AREA, the area of the unit sphere
//    in that dimension.
//
{
# define N_MAX 20

  static double area_vec[N_MAX] = {
     0.2000000000000000E+01,
     0.6283185307179586E+01,
     0.1256637061435917E+02,
     0.1973920880217872E+02,
     0.2631894506957162E+02,
     0.3100627668029982E+02,
     0.3307336179231981E+02,
     0.3246969701133415E+02,
     0.2968658012464836E+02,
     0.2550164039877345E+02,
     0.2072514267328890E+02,
     0.1602315322625507E+02,
     0.1183817381218268E+02,
     0.8389703410491089E+01,
     0.5721649212349567E+01,
     0.3765290085742291E+01,
     0.2396678817591364E+01,
     0.1478625959000308E+01,
     0.8858104195716824E+00,
     0.5161378278002812E+00 };

  static int n_vec[N_MAX] = {
     1,
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10,
    11,
    12,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    area = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    area = area_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sphere_unit_volume_values ( int &n_data, int &n, double &volume )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERE_UNIT_VOLUME_VALUES returns some volumes of the unit sphere in ND.
//
//  Discussion:
//
//    The formula for the volume of the unit sphere in N dimensions is
//
//      Volume(N) = 2 * PI^(N/2) / ( N * Gamma ( N / 2 ) )
//
//    This function satisfies the relationships:
//
//      Volume(N) = 2 * PI * Volume(N-2) / N
//      Volume(N) = Area(N) / N
//
//    Some values of the function include:
//
//       N  Volume
//
//       1    1
//       2    1        * PI
//       3  ( 4 /   3) * PI
//       4  ( 1 /   2) * PI^2
//       5  ( 8 /  15) * PI^2
//       6  ( 1 /   6) * PI^3
//       7  (16 / 105) * PI^3
//       8  ( 1 /  24) * PI^4
//       9  (32 / 945) * PI^4
//      10  ( 1 / 120) * PI^5
//
//    In Mathematica, the function can be evaluated by:
//
//      2 * Pi^(n/2) / ( n * Gamma[n/2] )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.
//    On input, if N_DATA is 0, the first test data is returned, and
//    N_DATA is set to the index of the test data.  On each subsequent
//    call, N_DATA is incremented and that test data is returned.  When
//    there is no more test data, N_DATA is set to 0.
//
//    Output, int &N, the spatial dimension.
//
//    Output, double &VOLUME, the volume of the unit
//    sphere in that dimension.
//
{
# define N_MAX 20

  static int n_vec[N_MAX] = {
     1,  2,
     3,  4,
     5,  6,
     7,  8,
     9, 10,
    11, 12,
    13, 14,
    15, 16,
    17, 18,
    19, 20 };

  static double volume_vec[N_MAX] = {
     0.2000000000000000E+01,
     0.3141592653589793E+01,
     0.4188790204786391E+01,
     0.4934802200544679E+01,
     0.5263789013914325E+01,
     0.5167712780049970E+01,
     0.4724765970331401E+01,
     0.4058712126416768E+01,
     0.3298508902738707E+01,
     0.2550164039877345E+01,
     0.1884103879389900E+01,
     0.1335262768854589E+01,
     0.9106287547832831E+00,
     0.5992645293207921E+00,
     0.3814432808233045E+00,
     0.2353306303588932E+00,
     0.1409811069171390E+00,
     0.8214588661112823E-01,
     0.4662160103008855E-01,
     0.2580689139001406E-01  };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    volume = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    volume = volume_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void spherical_harmonic_values ( int &n_data, int &l, int &m, double &theta,
  double &phi, double &yr, double &yi )

//****************************************************************************80
//
//  Purpose:
//
//    SPHERICAL_HARMONIC_VALUES returns values of spherical harmonic functions.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by
//
//      SphericalHarmonicY [ l, m, theta, phi ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 March 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Eric Weisstein,
//    CRC Concise Encyclopedia of Mathematics,
//    CRC Press, 1998.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &L, int &M, double &THETA, &PHI, the arguments
//    of the function.
//
//    Output, double &YR, &YI, the real and imaginary parts of
//    the function.
//
{
# define N_MAX 20

  static int l_vec[N_MAX] = {
     0,  1,  2,
     3,  4,  5,
     5,  5,  5,
     5,  4,  4,
     4,  4,  4,
     3,  3,  3,
     3,  3 };
  static int m_vec[N_MAX] = {
     0,  0,  1,
     2,  3,  5,
     4,  3,  2,
     1,  2,  2,
     2,  2,  2,
    -1, -1, -1,
    -1, -1 };
  static double phi_vec[N_MAX] = {
    0.1047197551196598E+01, 0.1047197551196598E+01, 0.1047197551196598E+01,
    0.1047197551196598E+01, 0.1047197551196598E+01, 0.6283185307179586E+00,
    0.6283185307179586E+00, 0.6283185307179586E+00, 0.6283185307179586E+00,
    0.6283185307179586E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
    0.7853981633974483E+00, 0.7853981633974483E+00, 0.7853981633974483E+00,
    0.4487989505128276E+00, 0.8975979010256552E+00, 0.1346396851538483E+01,
    0.1795195802051310E+01, 0.2243994752564138E+01 };
  static double theta_vec[N_MAX] = {
    0.5235987755982989E+00, 0.5235987755982989E+00, 0.5235987755982989E+00,
    0.5235987755982989E+00, 0.5235987755982989E+00, 0.2617993877991494E+00,
    0.2617993877991494E+00, 0.2617993877991494E+00, 0.2617993877991494E+00,
    0.2617993877991494E+00, 0.6283185307179586E+00, 0.1884955592153876E+01,
    0.3141592653589793E+01, 0.4398229715025711E+01, 0.5654866776461628E+01,
    0.3926990816987242E+00, 0.3926990816987242E+00, 0.3926990816987242E+00,
    0.3926990816987242E+00, 0.3926990816987242E+00 };
  static double yi_vec[N_MAX] = {
    0.0000000000000000E+00,  0.0000000000000000E+00, -0.2897056515173922E+00,
    0.1916222768312404E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
    0.3739289485283311E-02, -0.4219517552320796E-01,  0.1876264225575173E+00,
   -0.3029973424491321E+00,  0.4139385503112256E+00, -0.1003229830187463E+00,
    0.0000000000000000E+00, -0.1003229830187463E+00,  0.4139385503112256E+00,
   -0.1753512375142586E+00, -0.3159720118970196E+00, -0.3940106541811563E+00,
   -0.3940106541811563E+00, -0.3159720118970196E+00 };
  static double yr_vec[N_MAX] = {
   0.2820947917738781E+00,  0.4231421876608172E+00, -0.1672616358893223E+00,
  -0.1106331731112457E+00,  0.1354974113737760E+00,  0.5390423109043568E-03,
  -0.5146690442951909E-02,  0.1371004361349490E-01,  0.6096352022265540E-01,
  -0.4170400640977983E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
   0.0000000000000000E+00,  0.0000000000000000E+00,  0.0000000000000000E+00,
   0.3641205966137958E+00,  0.2519792711195075E+00,  0.8993036065704300E-01,
  -0.8993036065704300E-01, -0.2519792711195075E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    l = 0;
    m = 0;
    theta = 0.0;
    phi = 0.0;
    yr = 0.0;
    yi = 0.0;
  }
  else
  {
    l = l_vec[n_data-1];
    m = m_vec[n_data-1];
    theta = theta_vec[n_data-1];
    phi = phi_vec[n_data-1];
    yr = yr_vec[n_data-1];
    yi = yi_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void sqrt_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SQRT_VALUES returns some values of the square root function.
//
//  Discussion:
//
//    SQRT(X) = positive real number Y such that Y * Y = X.
//
//    In Mathematica, the function can be evaluated by:
//
//      Sqrt[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output double FX, the value of the function.
//
{
# define N_MAX 14

  static double fx_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.9000000040950000E-04,
     0.3000000000000000E+00,
     0.3162277660168379E+00,
     0.6324555320336759E+00,
     0.1000000000000000E+01,
     0.1414213562373095E+01,
     0.1732050807568877E+01,
     0.1772453850905516E+01,
     0.4358898943540674E+01,
     0.5385164807134504E+01,
     0.8426149773176359E+01,
     0.9848857801796105E+01,
     0.1111111106055556E+05 };

  static double x_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.8100000073710001E-08,
     0.9000000000000000E-01,
     0.1000000000000000E+00,
     0.4000000000000000E+00,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3141592653589793E+01,
     0.1900000000000000E+02,
     0.2900000000000000E+02,
     0.7100000000000000E+02,
     0.9700000000000000E+02,
     0.1234567890000000E+09 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void stirling1_values ( int &n_data, int &n, int &m, int &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING1_VALUES returns some values of the Stirling numbers, kind 1.
//
//  Discussion:
//
//    The absolute value of the Stirling number S1(N,M) gives the number
//    of permutations on N objects having exactly M cycles, while the
//    sign of the Stirling number records the sign (odd or even) of
//    the permutations.  For example, there are six permutations on 3 objects:
//
//      A B C   3 cycles (A) (B) (C)
//      A C B   2 cycles (A) (BC)
//      B A C   2 cycles (AB) (C)
//      B C A   1 cycle  (ABC)
//      C A B   1 cycle  (ABC)
//      C B A   2 cycles (AC) (B)
//
//    There are
//
//      2 permutations with 1 cycle, and S1(3,1) = 2
//      3 permutations with 2 cycles, and S1(3,2) = -3,
//      1 permutation with 3 cycles, and S1(3,3) = 1.
//
//    Since there are N! permutations of N objects, the sum of the absolute
//    values of the Stirling numbers in a given row,
//
//      sum ( 1 <= I <= N ) abs ( S1(N,I) ) = N!
//
//  First terms:
//
//    N/M:  1     2      3     4     5    6    7    8
//
//    1     1     0      0     0     0    0    0    0
//    2    -1     1      0     0     0    0    0    0
//    3     2    -3      1     0     0    0    0    0
//    4    -6    11     -6     1     0    0    0    0
//    5    24   -50     35   -10     1    0    0    0
//    6  -120   274   -225    85   -15    1    0    0
//    7   720 -1764   1624  -735   175  -21    1    0
//    8 -5040 13068 -13132  6769 -1960  322  -28    1
//
//    In Mathematica, the function can be evaluated by:
//
//      StirlingS1[n,m]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, &M, the arguments of the function.
//
//    Output, int &FX, the value of the function.
//
{
# define N_MAX 16

  static int fx_vec[N_MAX] = {
           0,
           1,
          -3,
          11,
         -50,
         274,
       -1764,
       13068,
     -109584,
     1026576,
      -13132,
        6769,
       -1960,
         322,
         -28,
           1 };

  static int m_vec[N_MAX] = {
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8 };

  static int n_vec[N_MAX] = {
     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    fx = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void stirling2_values ( int &n_data, int &n, int &m, int &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STIRLING2_VALUES returns some values of the Stirling numbers, kind 2.
//
//  Discussion:
//
//    S2(N,M) represents the number of distinct partitions of N elements
//    into M nonempty sets.  For a fixed N, the sum of the Stirling
//    numbers S2(N,M) is represented by B(N), called "Bell's number",
//    and represents the number of distinct partitions of N elements.
//
//    For example, with 4 objects, there are:
//
//    1 partition into 1 set:
//
//      (A,B,C,D)
//
//    7 partitions into 2 sets:
//
//      (A,B,C) (D)
//      (A,B,D) (C)
//      (A,C,D) (B)
//      (A) (B,C,D)
//      (A,B) (C,D)
//      (A,C) (B,D)
//      (A,D) (B,C)
//
//    6 partitions into 3 sets:
//
//      (A,B) (C) (D)
//      (A) (B,C) (D)
//      (A) (B) (C,D)
//      (A,C) (B) (D)
//      (A,D) (B) (C)
//      (A) (B,D) (C)
//
//    1 partition into 4 sets:
//
//      (A) (B) (C) (D)
//
//    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
//
//
//  First terms:
//
//    N/M: 1    2    3    4    5    6    7    8
//
//    1    1    0    0    0    0    0    0    0
//    2    1    1    0    0    0    0    0    0
//    3    1    3    1    0    0    0    0    0
//    4    1    7    6    1    0    0    0    0
//    5    1   15   25   10    1    0    0    0
//    6    1   31   90   65   15    1    0    0
//    7    1   63  301  350  140   21    1    0
//    8    1  127  966 1701 1050  266   28    1
//
//    In Mathematica, the function can be evaluated by:
//
//      StirlingS2[n,m]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, &M, the arguments of the function.
//
//    Output, int &FX, the value of the function.
//
{
# define N_MAX 16

  static int fx_vec[N_MAX] = {
           0,
           1,
           3,
           7,
          15,
          31,
          63,
         127,
         255,
         511,
         966,
        1701,
        1050,
         266,
          28,
           1 };

  static int m_vec[N_MAX] = {
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7, 8 };

  static int n_vec[N_MAX] = {
     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 8, 8, 8, 8, 8, 8 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    m = 0;
    fx = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    m = m_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void stromgen_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STROMGEN_VALUES returns some values of the Stromgen function.
//
//  Discussion:
//
//    The function is defined by:
//
//      STROMGEN(X) = integral ( 0 <= t <= X ) t^7 * exp(2*t) / (exp(t)-1)^3 dt
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.21901065985698662316E-15,
     0.22481399438625244761E-12,
     0.23245019579558857124E-09,
     0.24719561475975007037E-06,
     0.28992610989833245669E-03,
     0.10698146390809715091E-01,
     0.89707650964424730705E-01,
     0.40049605719592888440E+00,
     0.30504104398079096598E+01,
     0.11367704858439426431E+02,
     0.12960679405324786954E+02,
     0.18548713944748505675E+02,
     0.27866273821903121400E+02,
     0.51963334071699323351E+02,
     0.10861016747891228129E+03,
     0.15378903316556621624E+03,
     0.19302665532558721516E+03,
     0.19636850166006541482E+03,
     0.19651946766008214217E+03,
     0.19651956920868316152E+03 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0078125000E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.1250000000E+00,
       4.5000000000E+00,
       5.0000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void struve_h0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_H0_VALUES returns some values of the Struve H0 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      HO(x) = 2/pi * integral ( 0 <= t <= pi/2 ) sin ( x * cos ( t ) ) dt
//
//    In Mathematica, the function can be evaluated by:
//
//      StruveH[0,x]
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.12433974658847434366E-02,
     -0.49735582423748415045E-02,
      0.39771469054536941564E-01,
     -0.15805246001653314198E+00,
      0.56865662704828795099E+00,
      0.66598399314899916605E+00,
      0.79085884950809589255E+00,
     -0.13501457342248639716E+00,
      0.20086479668164503137E+00,
     -0.11142097800261991552E+00,
     -0.17026804865989885869E+00,
     -0.13544931808186467594E+00,
      0.94393698081323450897E-01,
     -0.10182482016001510271E+00,
      0.96098421554162110012E-01,
     -0.85337674826118998952E-01,
     -0.76882290637052720045E-01,
      0.47663833591418256339E-01,
     -0.70878751689647343204E-01,
      0.65752908073352785368E-01 };

  static double x_vec[N_MAX] = {
        0.0019531250E+00,
       -0.0078125000E+00,
        0.0625000000E+00,
       -0.2500000000E+00,
        1.0000000000E+00,
        1.2500000000E+00,
        2.0000000000E+00,
       -4.0000000000E+00,
        7.5000000000E+00,
       11.0000000000E+00,
       11.5000000000E+00,
      -16.0000000000E+00,
       20.0000000000E+00,
       25.0000000000E+00,
      -30.0000000000E+00,
       50.0000000000E+00,
       75.0000000000E+00,
      -80.0000000000E+00,
      100.0000000000E+00,
     -125.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void struve_h1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_H1_VALUES returns some values of the Struve H1 function.
//
//  Discussion:
//
//    The function is defined by:
//
//      H1(x) = 2*x/pi * integral ( 0 <= t <= pi/2 )
//        sin ( x * cos ( t ) )^2 * sin ( t ) dt
//
//    In Mathematica, the function can be evaluated by:
//
//      StruveH[1,x]
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.80950369576367526071E-06,
     0.12952009724113229165E-04,
     0.82871615165407083021E-03,
     0.13207748375849572564E-01,
     0.19845733620194439894E+00,
     0.29853823231804706294E+00,
     0.64676372828356211712E+00,
     0.10697266613089193593E+01,
     0.38831308000420560970E+00,
     0.74854243745107710333E+00,
     0.84664854642567359993E+00,
     0.58385732464244384564E+00,
     0.80600584524215772824E+00,
     0.53880362132692947616E+00,
     0.72175037834698998506E+00,
     0.58007844794544189900E+00,
     0.60151910385440804463E+00,
     0.70611511147286827018E+00,
     0.61631110327201338454E+00,
     0.62778480765443656489E+00 };

  static double x_vec[N_MAX] = {
        0.0019531250E+00,
       -0.0078125000E+00,
        0.0625000000E+00,
       -0.2500000000E+00,
        1.0000000000E+00,
        1.2500000000E+00,
        2.0000000000E+00,
       -4.0000000000E+00,
        7.5000000000E+00,
       11.0000000000E+00,
       11.5000000000E+00,
      -16.0000000000E+00,
       20.0000000000E+00,
       25.0000000000E+00,
      -30.0000000000E+00,
       50.0000000000E+00,
       75.0000000000E+00,
      -80.0000000000E+00,
      100.0000000000E+00,
     -125.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void struve_l0_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_L0_VALUES returns some values of the Struve L0 function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      StruveL[0,x]
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
      0.12433985199262820188E-02,
     -0.19896526647882937004E-01,
      0.79715713253115014945E-01,
     -0.32724069939418078025E+00,
      0.71024318593789088874E+00,
      0.19374337579914456612E+01,
     -0.11131050203248583431E+02,
      0.16850062034703267148E+03,
     -0.28156522493745948555E+04,
      0.89344618796978400815E+06,
      0.11382025002851451057E+07,
     -0.23549701855860190304E+07,
      0.43558282527641046718E+08,
      0.49993516476037957165E+09,
     -0.57745606064408041689E+10,
      0.78167229782395624524E+12,
     -0.14894774793419899908E+17,
      0.29325537838493363267E+21,
      0.58940770556098011683E+25,
     -0.12015889579125463605E+30 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
      -0.0312500000E+00,
       0.1250000000E+00,
      -0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
      -4.0000000000E+00,
       7.0000000000E+00,
     -10.0000000000E+00,
      16.0000000000E+00,
      16.2500000000E+00,
     -17.0000000000E+00,
      20.0000000000E+00,
      22.5000000000E+00,
     -25.0000000000E+00,
      30.0000000000E+00,
     -40.0000000000E+00,
      50.0000000000E+00,
      60.0000000000E+00,
     -70.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void struve_l1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STRUVE_L1_VALUES returns some values of the Struve L1 function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      StruveL[1,x]
//
//    The data was reported by McLeod.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.80950410749865126939E-06,
     0.20724649092571514607E-03,
     0.33191834066894516744E-02,
     0.53942182623522663292E-01,
     0.22676438105580863683E+00,
     0.11027597873677158176E+01,
     0.91692778117386847344E+01,
     0.15541656652426660966E+03,
     0.26703582852084829694E+04,
     0.86505880175304633906E+06,
     0.11026046613094942620E+07,
     0.22846209494153934787E+07,
     0.42454972750111979449E+08,
     0.48869614587997695539E+09,
     0.56578651292431051863E+10,
     0.76853203893832108948E+12,
     0.14707396163259352103E+17,
     0.29030785901035567967E+21,
     0.58447515883904682813E+25,
     0.11929750788892311875E+30 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
      -0.0312500000E+00,
       0.1250000000E+00,
      -0.5000000000E+00,
       1.0000000000E+00,
       2.0000000000E+00,
      -4.0000000000E+00,
       7.0000000000E+00,
     -10.0000000000E+00,
      16.0000000000E+00,
      16.2500000000E+00,
     -17.0000000000E+00,
      20.0000000000E+00,
      22.5000000000E+00,
     -25.0000000000E+00,
      30.0000000000E+00,
     -40.0000000000E+00,
      50.0000000000E+00,
      60.0000000000E+00,
     -70.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void student_cdf_values ( int &n_data, double &c, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_CDF_VALUES returns some values of the Student CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = StudentTDistribution [ c ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 November 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &C, is usually called the number of
//    degrees of freedom of the distribution.  C is typically an
//    integer, but that is not essential.  It is required that
//    C be strictly positive.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 13

  static double c_vec[N_MAX] = {
    1.0, 2.0, 3.0, 4.0,
    5.0, 2.0, 5.0, 2.0,
    5.0, 2.0, 3.0, 4.0,
    5.0 };

  static double fx_vec[N_MAX] = {
     0.6000231200328521,
     0.6001080279134390,
     0.6001150934648930,
     0.6000995134721354,
     0.5999341989834830,
     0.7498859393137811,
     0.7500879487671045,
     0.9500004222186464,
     0.9499969138365968,
     0.9900012348724744,
     0.9900017619355059,
     0.9900004567580596,
     0.9900007637471291 };

  static double x_vec[N_MAX] = {
     0.325,
     0.289,
     0.277,
     0.271,
     0.267,
     0.816,
     0.727,
     2.920,
     2.015,
     6.965,
     4.541,
     3.747,
     3.365 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    c = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    c = c_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void student_noncentral_cdf_values ( int &n_data, int &df, double &lambda,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    STUDENT_NONCENTRAL_CDF_VALUES returns values of the noncentral Student CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = NoncentralStudentTDistribution [ df, lambda ]
//      CDF [ dist, x ]
//
//    Mathematica seems to have some difficulty computing this function
//    to the desired number of digits.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &DF, double &LAMBDA, the parameters of the
//    function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 30

  static int df_vec[N_MAX] = {
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
     1,  2,  3,
    15, 20, 25,
     1,  2,  3,
    10, 10, 10,
    10, 10, 10,
    10, 10, 10 };

  static double fx_vec[N_MAX] = {
     0.8975836176504333E+00,
     0.9522670169E+00,
     0.9711655571887813E+00,
     0.8231218864E+00,
     0.9049021510E+00,
     0.9363471834E+00,
     0.7301025986E+00,
     0.8335594263E+00,
     0.8774010255E+00,
     0.5248571617E+00,
     0.6293856597E+00,
     0.6800271741E+00,
     0.20590131975E+00,
     0.2112148916E+00,
     0.2074730718E+00,
     0.9981130072E+00,
     0.9994873850E+00,
     0.9998391562E+00,
     0.168610566972E+00,
     0.16967950985E+00,
     0.1701041003E+00,
     0.9247683363E+00,
     0.7483139269E+00,
     0.4659802096E+00,
     0.9761872541E+00,
     0.8979689357E+00,
     0.7181904627E+00,
     0.9923658945E+00,
     0.9610341649E+00,
     0.8688007350E+00 };

  static double lambda_vec[N_MAX] = {
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.5E+00,
     0.5E+00,
     0.5E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     2.0E+00,
     2.0E+00,
     2.0E+00,
     4.0E+00,
     4.0E+00,
     4.0E+00,
     7.0E+00,
     7.0E+00,
     7.0E+00,
     1.0E+00,
     1.0E+00,
     1.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00,
     2.0E+00,
     3.0E+00,
     4.0E+00 };

  static double x_vec[N_MAX] = {
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
      3.00E+00,
     15.00E+00,
     15.00E+00,
     15.00E+00,
      0.05E+00,
      0.05E+00,
      0.05E+00,
      4.00E+00,
      4.00E+00,
      4.00E+00,
      5.00E+00,
      5.00E+00,
      5.00E+00,
      6.00E+00,
      6.00E+00,
      6.00E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    df = 0;
    lambda = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    df = df_vec[n_data-1];
    lambda = lambda_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void subfactorial_values ( int &n_data, int &n, int &fn )

//****************************************************************************80
//
//  Purpose:
//
//    SUBFACTORIAL_VALUES returns values of the subfactorial function.
//
//  Discussion:
//
//    The subfactorial function Subfactorial(N) counts the number of
//    permutations of N objects which leave no object unchanged.
//
//    Such a permutation is known as a derangement.
//
//    In Mathematica, the function can be evaluated by:
//
//      << DiscreteMath`CombinatorialFunctions`
//      Subfactorial[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 March 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the function.
//
//    Output, int &FN, the value of the function.
//
{
# define N_MAX 13

  static int fn_vec[N_MAX] = {
            1,
            0,
            1,
            2,
            9,
           44,
          265,
         1854,
        14833,
       133496,
      1334961,
     14684570,
    176214841 };

  static int n_vec[N_MAX] = {
     0,  1,  2,  3,
     4,  5,  6,  7,
     8,  9, 10, 11,
    12 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    fn = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    fn = fn_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void surten_values ( int &n_data, double &tc, double &sigma )

//****************************************************************************80
//
//  Purpose:
//
//    SURTEN_VALUES returns some values of the surface tension.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, pages 267.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double *TC, the temperature, in degrees Celsius.
//
//    Output, double *SIGMA, the surface tension,
//    in Pascal * m = Newton / m.
//
{
# define N_MAX 14

  static double sigma_vec[N_MAX] = {
     74.22E+00,
     72.74E+00,
     71.20E+00,
     69.60E+00,
     67.95E+00,
     58.92E+00,
     48.75E+00,
     37.68E+00,
     26.05E+00,
     14.37E+00,
      8.78E+00,
      3.67E+00,
      0.40E+00,
      0.00E+00 };

  static double tc_vec[N_MAX] = {
      10.000E+00,
      20.000E+00,
      30.000E+00,
      40.000E+00,
      50.000E+00,
     100.000E+00,
     150.000E+00,
     200.000E+00,
     250.000E+00,
     300.000E+00,
     325.000E+00,
     350.000E+00,
     370.000E+00,
     373.976E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    sigma = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    sigma = sigma_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void synch1_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SYNCH1_VALUES returns some values of the synchrotron radiation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      SYNCH1(x) = x * integral ( x <= t < +oo ) K(5/3)(t) dt
//
//    where K(5/3) is a modified Bessel function of order 5/3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
       0.26514864547487397044E+00,
       0.62050129979079045645E+00,
       0.85112572132368011206E+00,
       0.87081914687546885094E+00,
       0.65142281535536396975E+00,
       0.45064040920322354579E+00,
       0.30163590285073940285E+00,
       0.19814490804441305867E+00,
       0.12856571000906381300E+00,
       0.52827396697866818297E-01,
       0.42139298471720305542E-01,
       0.21248129774981984268E-01,
       0.13400258907505536491E-01,
       0.84260797314108699935E-02,
       0.12884516186754671469E-02,
       0.19223826430086897418E-03,
       0.28221070834007689394E-04,
       0.15548757973038189372E-05,
       0.11968634456097453636E-07,
       0.89564246772237127742E-10 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      25.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void synch2_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    SYNCH2_VALUES returns some values of the synchrotron radiation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      SYNCH2(x) = x * K(2/3)(x)
//
//    where K(2/3) is a modified Bessel function of order 2/3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.13430727275667378338E+00,
     0.33485265272424176976E+00,
     0.50404224110911078651E+00,
     0.60296523236016785113E+00,
     0.49447506210420826699E+00,
     0.36036067860473360389E+00,
     0.24967785497625662113E+00,
     0.16813830542905833533E+00,
     0.11117122348556549832E+00,
     0.46923205826101330711E-01,
     0.37624545861980001482E-01,
     0.19222123172484106436E-01,
     0.12209535343654701398E-01,
     0.77249644268525771866E-02,
     0.12029044213679269639E-02,
     0.18161187569530204281E-03,
     0.26884338006629353506E-04,
     0.14942212731345828759E-05,
     0.11607696854385161390E-07,
     0.87362343746221526073E-10 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      12.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      25.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tan_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TAN_VALUES returns some values of the tangent function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Tan[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 15

  static double fx_vec[N_MAX] = {
     0.00000000000000000000,
     0.26794919243112270647,
     0.54630248984379051326,
     0.57735026918962576451,
     1.0000000000000000000,
     1.5574077246549022305,
     1.7320508075688772935,
     3.7320508075688772935,
     7.5957541127251504405,
    15.257051688265539110,
    -2.1850398632615189916,
    -0.14254654307427780530,
     0.0000000000000000000,
     1.1578212823495775831,
    -3.3805150062465856370 };

  static double x_vec[N_MAX] = {
    0.00000000000000000000,
    0.26179938779914943654,
    0.50000000000000000000,
    0.52359877559829887308,
    0.78539816339744830962,
    1.0000000000000000000,
    1.0471975511965977462,
    1.3089969389957471827,
    1.4398966328953219010,
    1.5053464798451092601,
    2.0000000000000000000,
    3.0000000000000000000,
    3.1415926535897932385,
    4.0000000000000000000,
    5.0000000000000000000 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tanh_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TANH_VALUES returns some values of the hyperbolic tangent function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Tanh[x]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 June 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 18

  static double fx_vec[N_MAX] = {
   -0.99990920426259513121,
   -0.76159415595576488812,
    0.00000000000000000000,
    0.099667994624955817118,
    0.19737532022490400074,
    0.29131261245159090582,
    0.37994896225522488527,
    0.46211715726000975850,
    0.53704956699803528586,
    0.60436777711716349631,
    0.66403677026784896368,
    0.71629787019902442081,
    0.76159415595576488812,
    0.96402758007581688395,
    0.99505475368673045133,
    0.99932929973906704379,
    0.99990920426259513121,
    0.99999999587769276362 };

  static double x_vec[N_MAX] = {
   -5.0,
   -1.0,
    0.0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
   10.0 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tau_values ( int &n_data, int &n, int &c )

//****************************************************************************80
//
//  Purpose:
//
//    TAU_VALUES returns some values of the Tau function.
//
//  Discussion:
//
//    TAU(N) is the number of divisors of N, including 1 and N.
//
//    In Mathematica, the function can be evaluated by:
//
//      DivisorSigma[1,n]
//
//  First values:
//
//     N   TAU(N)
//
//     1    1
//     2    2
//     3    2
//     4    3
//     5    2
//     6    4
//     7    2
//     8    4
//     9    3
//    10    4
//    11    2
//    12    6
//    13    2
//    14    4
//    15    4
//    16    5
//    17    2
//    18    6
//    19    2
//    20    6
//
//  Formula:
//
//    If the prime factorization of N is
//
//      N = P1**E1 * P2**E2 * ... * PM**EM,
//
//    then
//
//      TAU(N) = ( E1 + 1 ) * ( E2 + 1 ) * ... * ( EM + 1 ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the Tau function.
//
//    Output, int &C, the value of the Tau function.
//
{
# define N_MAX 20

  static int c_vec[N_MAX] = {
    1,  2,  2,  3,  2,  4,  2,  4,  3,  4,
    2, 12, 12,  4, 18, 24,  2,  8, 14, 28 };

  static int n_vec[N_MAX] = {
      1,   2,   3,   4,   5,   6,   7,   8,   9,  10,
     23,  72, 126, 226, 300, 480, 521, 610, 832, 960 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    c = 0;
  }
  else
  {
    n = n_vec[n_data-1];
    c = c_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void thercon_values ( int &n_data, double &tc, double &p, double &lambda )

//****************************************************************************80
//
//  Purpose:
//
//    THERCON_VALUES returns some values of the thermal conductivity.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, page 264.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &LAMBDA, the thermal conductivity, in
//    mW/(m degrees Kelvin).
//
{
# define N_MAX 35

  static double lambda_vec[N_MAX] = {
     561.00E+00,
     561.30E+00,
     561.50E+00,
     562.40E+00,
     563.70E+00,
     565.10E+00,
     566.50E+00,
     567.90E+00,
     569.30E+00,
     570.60E+00,
     572.00E+00,
     573.40E+00,
     574.80E+00,
     576.10E+00,
     577.50E+00,
     580.20E+00,
     582.90E+00,
     585.50E+00,
     588.10E+00,
     590.70E+00,
     593.30E+00,
     595.80E+00,
     598.30E+00,
     603.10E+00,
     607.80E+00,
     612.20E+00,
     607.20E+00,
     643.60E+00,
     666.80E+00,
      25.08E+00,
      28.85E+00,
      33.28E+00,
      54.76E+00,
      79.89E+00,
     107.30E+00 };

  static double p_vec[N_MAX] = {
        1.0E+00,
        5.0E+00,
       10.0E+00,
       25.0E+00,
       50.0E+00,
       75.0E+00,
      100.0E+00,
      125.0E+00,
      150.0E+00,
      175.0E+00,
      200.0E+00,
      225.0E+00,
      250.0E+00,
      275.0E+00,
      300.0E+00,
      350.0E+00,
      400.0E+00,
      450.0E+00,
      500.0E+00,
      550.0E+00,
      600.0E+00,
      650.0E+00,
      700.0E+00,
      800.0E+00,
      900.0E+00,
     1000.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00 };

  static double tc_vec[N_MAX] = {
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
      25.0E+00,
      50.0E+00,
      75.0E+00,
     100.0E+00,
     150.0E+00,
     200.0E+00,
     400.0E+00,
     600.0E+00,
     800.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    lambda = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    lambda = lambda_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void three_j_values ( int &n_data, double &j1, double &j2, double &j3,
  double &m1, double &m2, double &m3, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    THREE_J_VALUES returns some values of the Wigner 3J function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      ThreeJSymbol[{j1,m1},{j2,m2},{j3,m3}]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &J1, &J2, &J3, &M1, &M2, &M3, the arguments
//    of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 8

  static double fx_vec[N_MAX] = {
     0.2788866755113585,
    -0.09534625892455923,
    -0.06741998624632421,
     0.1533110351679666,
    -0.1564465546936860,
     0.1099450412156551,
    -0.05536235693131719,
     0.01799835451137786 };
  static double j1_vec[N_MAX] = {
    1.0,
    2.0,
    3.0,
    4.0,
    5.0,
    6.0,
    7.0,
    8.0 };
  static double j2_vec[N_MAX] = {
    4.5,
    4.5,
    4.5,
    4.5,
    4.5,
    4.5,
    4.5,
    4.5 };
  static double j3_vec[N_MAX] = {
    3.5,
    3.5,
    3.5,
    3.5,
    3.5,
    3.5,
    3.5,
    3.5 };
  static double m1_vec[N_MAX] = {
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0 };
  static double m2_vec[N_MAX] = {
    -3.5,
    -3.5,
    -3.5,
    -3.5,
    -3.5,
    -3.5,
    -3.5,
    -3.5 };
  static double m3_vec[N_MAX] = {
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5,
    2.5 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    j1 = 0.0;
    j2 = 0.0;
    j3 = 0.0;
    m1 = 0.0;
    m2 = 0.0;
    m3 = 0.0;
    fx = 0.0;
  }
  else
  {
    j1 = j1_vec[n_data-1];
    j2 = j2_vec[n_data-1];
    j3 = j3_vec[n_data-1];
    m1 = m1_vec[n_data-1];
    m2 = m2_vec[n_data-1];
    m3 = m3_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void tran02_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN02_VALUES returns some values of the order 2 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN02(x) = integral ( 0 <= t <= x ) t^2 exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.19531247930394515480E-02,
     0.31249152314331109004E-01,
     0.12494577194783451032E+00,
     0.49655363615640595865E+00,
     0.97303256135517012845E+00,
     0.14121978695932525805E+01,
     0.18017185674405776809E+01,
     0.21350385339277043015E+01,
     0.24110500490169534620E+01,
     0.28066664045631179931E+01,
     0.28777421863296234131E+01,
     0.30391706043438554330E+01,
     0.31125074928667355940E+01,
     0.31656687817738577185E+01,
     0.32623520367816009184E+01,
     0.32843291144979517358E+01,
     0.32897895167775788137E+01,
     0.32898672226665499687E+01,
     0.32898681336064325400E+01,
     0.32898681336964528724E+01 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran03_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN03_VALUES returns some values of the order 3 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN03(x) = integral ( 0 <= t <= x ) t^3 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.19073483296476379584E-05,
     0.48826138243180786081E-03,
     0.78074163848431205820E-02,
     0.12370868718812031049E+00,
     0.47984100657241749994E+00,
     0.10269431622039754738E+01,
     0.17063547219458658863E+01,
     0.24539217444475937661E+01,
     0.32106046629422467723E+01,
     0.45792174372291563703E+01,
     0.48722022832940370805E+01,
     0.56143866138422732286E+01,
     0.59984455864575470009E+01,
     0.63033953673480961120E+01,
     0.69579908688361166266E+01,
     0.71503227120085929750E+01,
     0.72110731475871876393E+01,
     0.72123221966388461839E+01,
     0.72123414161609465119E+01,
     0.72123414189575656868E+01 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran04_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN04_VALUES returns some values of the order 4 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN04(x) = integral ( 0 <= t <= x ) t^4 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.24835263919461834041E-08,
     0.10172029353616724881E-04,
     0.65053332405940765479E-03,
     0.41150448004155727767E-01,
     0.31724404523442648241E+00,
     0.10079442901142373591E+01,
     0.22010881024333408363E+01,
     0.38846508619156545210E+01,
     0.59648223973714765245E+01,
     0.10731932392998622219E+02,
     0.11940028876819364777E+02,
     0.15359784316882182982E+02,
     0.17372587633093742893E+02,
     0.19122976016053166969E+02,
     0.23583979156921941515E+02,
     0.25273667677030441733E+02,
     0.25955198214572256372E+02,
     0.25975350935212241910E+02,
     0.25975757522084093747E+02,
     0.25975757609067315288E+02 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran05_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN05_VALUES returns some values of the order 5 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN05(x) = integral ( 0 <= t <= x ) t^5 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.36379780361036116971E-11,
     0.23840564453948442379E-06,
     0.60982205372226969189E-04,
     0.15410004586376649337E-01,
     0.23661587923909478926E+00,
     0.11198756851307629651E+01,
     0.32292901663684049171E+01,
     0.70362973105160654056E+01,
     0.12770557691044159511E+02,
     0.29488339015245845447E+02,
     0.34471340540362254586E+02,
     0.50263092218175187785E+02,
     0.60819909101127165207E+02,
     0.70873334429213460498E+02,
     0.10147781242977788097E+03,
     0.11638074540242071077E+03,
     0.12409623901262967878E+03,
     0.12442270155632550228E+03,
     0.12443132790838589548E+03,
     0.12443133061720432435E+03 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran06_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN06_VALUES returns some values of the order 6 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN06(x) = integral ( 0 <= t <= x ) t^6 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.56843405953641209574E-14,
     0.59601180165247401484E-08,
     0.60978424397580572815E-05,
     0.61578909866319494394E-02,
     0.18854360275680840514E+00,
     0.13319251347921659134E+01,
     0.50857202271697616755E+01,
     0.13729222365466557122E+02,
     0.29579592481641441292E+02,
     0.88600835706899853768E+02,
     0.10916037113373004909E+03,
     0.18224323749575359518E+03,
     0.23765383125586756031E+03,
     0.29543246745959381136E+03,
     0.50681244381280455592E+03,
     0.63878231134946125623E+03,
     0.72699203556994876111E+03,
     0.73230331643146851717E+03,
     0.73248692015882096369E+03,
     0.73248700462879996604E+03 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran07_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN07_VALUES returns some values of the order 7 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN07(x) = integral ( 0 <= t <= x ) t^7 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.92518563327283409427E-17,
     0.15521095556949867541E-09,
     0.63516238373841716290E-06,
     0.25638801246626135714E-02,
     0.15665328993811649746E+00,
     0.16538225039181097423E+01,
     0.83763085709508211054E+01,
     0.28078570717830763747E+02,
     0.72009676046751991365E+02,
     0.28174905701691911450E+03,
     0.36660227975327792529E+03,
     0.70556067982603601123E+03,
     0.99661927562755629434E+03,
     0.13288914430417403901E+04,
     0.27987640273169129925E+04,
     0.39721376409416504325E+04,
     0.49913492839319899726E+04,
     0.50781562639825019000E+04,
     0.50820777202028708434E+04,
     0.50820803580047164618E+04 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran08_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN08_VALUES returns some values of the order 8 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN08(x) = integral ( 0 <= t <= x ) t^8 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.15488598634539359463E-19,
     0.41574269117845953797E-11,
     0.68050651245227411689E-07,
     0.10981703519563009836E-02,
     0.13396432776187883834E+00,
     0.21153387806998617182E+01,
     0.14227877028750735641E+02,
     0.59312061431647843226E+02,
     0.18139614577043147745E+03,
     0.93148001928992220863E+03,
     0.12817928112604611804E+04,
     0.28572838386329242218E+04,
     0.43872971687877730010E+04,
     0.62993229139406657611E+04,
     0.16589426277154888511E+05,
     0.27064780798797398935E+05,
     0.38974556062543661284E+05,
     0.40400240716905025786E+05,
     0.40484316504120655568E+05,
     0.40484399001892184901E+05 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tran09_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRAN09_VALUES returns some values of the order 9 transportation function.
//
//  Discussion:
//
//    The function is defined by:
//
//      TRAN09(x) = integral ( 0 <= t <= x ) t^9 * exp(t) / ( exp(t) - 1 )^2 dt
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Allan McLeod,
//    Algorithm 757:
//    MISCFUN: A software package to compute uncommon special functions,
//    ACM Transactions on Mathematical Software,
//    Volume 22, Number 3, September 1996, pages 288-301.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 20

  static double fx_vec[N_MAX] = {
     0.26469772870084897671E-22,
     0.11367943653594246210E-12,
     0.74428246255329800255E-08,
     0.48022728485415366194E-03,
     0.11700243014358676725E+00,
     0.27648973910899914391E+01,
     0.24716631405829192997E+02,
     0.12827119828849828583E+03,
     0.46842894800662208986E+03,
     0.31673967371627895718E+04,
     0.46140886546630195390E+04,
     0.11952718545392302185E+05,
     0.20001612666477027728E+05,
     0.31011073271851366554E+05,
     0.10352949905541130133E+06,
     0.19743173017140591390E+06,
     0.33826030414658460679E+06,
     0.36179607036750755227E+06,
     0.36360622124777561525E+06,
     0.36360880558827162725E+06 };

  static double x_vec[N_MAX] = {
       0.0019531250E+00,
       0.0312500000E+00,
       0.1250000000E+00,
       0.5000000000E+00,
       1.0000000000E+00,
       1.5000000000E+00,
       2.0000000000E+00,
       2.5000000000E+00,
       3.0000000000E+00,
       4.0000000000E+00,
       4.2500000000E+00,
       5.0000000000E+00,
       5.5000000000E+00,
       6.0000000000E+00,
       8.0000000000E+00,
      10.0000000000E+00,
      15.0000000000E+00,
      20.0000000000E+00,
      30.0000000000E+00,
      50.0000000000E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x  = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void trigamma_values ( int &n_data, double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    TRIGAMMA_VALUES returns some values of the TriGamma function.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      PolyGamma[1,x]
//
//    TriGamma(X) = d^2 ln ( Gamma ( X ) ) / d X^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 11

  static double fx_vec[N_MAX] = {
    0.1644934066848226E+01,
    0.1433299150792759E+01,
    0.1267377205423779E+01,
    0.1134253434996619E+01,
    0.1025356590529597E+01,
    0.9348022005446793E+00,
    0.8584318931245799E+00,
    0.7932328301639984E+00,
    0.7369741375017002E+00,
    0.6879720582426356E+00,
    0.6449340668482264E+00 };

  static double x_vec[N_MAX] = {
     1.0E+00,
     1.1E+00,
     1.2E+00,
     1.3E+00,
     1.4E+00,
     1.5E+00,
     1.6E+00,
     1.7E+00,
     1.8E+00,
     1.9E+00,
     2.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void tsat_values ( int &n_data, double &p, double &tc )

//****************************************************************************80
//
//  Purpose:
//
//    TSAT_VALUES returns some values of the saturation temperature.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, pages 16-22.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &TC, the saturation temperature, in
//    degrees Celsius.
//
{
# define N_MAX 20

  static double p_vec[N_MAX] = {
       0.0061173E+00,
       0.012E+00,
       0.025E+00,
       0.055E+00,
       0.080E+00,
       0.110E+00,
       0.160E+00,
       0.250E+00,
       0.500E+00,
       0.750E+00,
       1.000E+00,
       1.500E+00,
       2.000E+00,
       5.000E+00,
      10.000E+00,
      20.000E+00,
      50.000E+00,
     100.000E+00,
     200.000E+00,
     220.550E+00 };

  static double tc_vec[N_MAX] = {
       0.010E+00,
       9.655E+00,
      21.080E+00,
      34.589E+00,
      41.518E+00,
      47.695E+00,
      55.327E+00,
      64.980E+00,
      81.339E+00,
      91.783E+00,
      99.632E+00,
     111.378E+00,
     120.443E+00,
     151.866E+00,
     179.916E+00,
     212.417E+00,
     263.977E+00,
     311.031E+00,
     365.800E+00,
     373.976E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    p = 0.0;
    tc = 0.0;
  }
  else
  {
    p = p_vec[n_data-1];
    tc = tc_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void van_der_corput_values ( int &n_data, int &base, int &seed, double &value )

//****************************************************************************80
//
//  Purpose:
//
//    VAN_DER_CORPUT_VALUES returns some values of the van der Corput sequence.
//
//  Discussion:
//
//    The van der Corput sequence is often used to generate a "subrandom"
//    sequence of points which have a better covering property
//    than pseudorandom points.
//
//    The van der Corput sequence generates a sequence of points in [0,1]
//    which (theoretically) never repeats.  Except for SEED = 0, the
//    elements of the van der Corput sequence are strictly between 0 and 1.
//
//    The van der Corput sequence writes an int *in a given base B,
//    and then its digits are "reflected" about the decimal point.
//    This maps the numbers from 1 to N into a set of numbers in [0,1],
//    which are especially nicely distributed if N is one less
//    than a power of the base.
//
//    Hammersley suggested generating a set of N nicely distributed
//    points in two dimensions by setting the first component of the
//    Ith point to I/N, and the second to the van der Corput
//    value of I in base 2.
//
//    Halton suggested that in many cases, you might not know the number
//    of points you were generating, so Hammersley's formulation was
//    not ideal.  Instead, he suggested that to generate a nicely
//    distributed sequence of points in M dimensions, you simply
//    choose the first M primes, P(1:M), and then for the J-th component of
//    the I-th point in the sequence, you compute the van der Corput
//    value of I in base P(J).
//
//    Thus, to generate a Halton sequence in a 2 dimensional space,
//    it is typical practice to generate a pair of van der Corput sequences,
//    the first with prime base 2, the second with prime base 3.
//    Similarly, by using the first K primes, a suitable sequence
//    in K-dimensional space can be generated.
//
//    The generation is quite simple.  Given an int *SEED, the expansion
//    of SEED in base BASE is generated.  Then, essentially, the result R
//    is generated by writing a decimal point followed by the digits of
//    the expansion of SEED, in reverse order.  This decimal value is actually
//    still in base BASE, so it must be properly interpreted to generate
//    a usable value.
//
//  Example:
//
//    BASE = 2
//
//    SEED     SEED      van der Corput
//    decimal  binary    binary   decimal
//    -------  ------    ------   -------
//        0  =     0  =>  .0     = 0.0;
//        1  =     1  =>  .1     = 0.5
//        2  =    10  =>  .01    = 0.25
//        3  =    11  =>  .11    = 0.75
//        4  =   100  =>  .001   = 0.125
//        5  =   101  =>  .101   = 0.625
//        6  =   110  =>  .011   = 0.375
//        7  =   111  =>  .111   = 0.875
//        8  =  1000  =>  .0001  = 0.0625
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    28 June 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
//
//    John Hammersley,
//    Monte Carlo methods for solving multivariable problems,
//    Proceedings of the New York Academy of Science,
//    Volume 86, pages 844-874, 1960.
//
//    J G van der Corput,
//    Verteilungsfunktionen,
//    Proc Akad Amsterdam,
//    Volume 38, 1935,
//    Volume 39, 1936.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &BASE, the base of the sequence.
//
//    Output, int &SEED, the index of the element of the sequence.
//
//    Output, double &VALUE, the value of the SEED-th element of the
//    van der Corput sequence in base BASE.
//
{
# define N_MAX 75

  static int base_vec[N_MAX] = {
     2,   2,   2,   2,   2,
     2,   2,   2,   2,   3,
     3,   3,   3,   3,   3,
     3,   3,   3,   4,   4,
     4,   4,   4,   4,   4,
     4,   4,   2,   3,   4,
     5,   7,  11,  13,   2,
     3,   4,   5,   7,  11,
    13,   2,   3,   4,   5,
     7,  11,  13,   2,   3,
     4,   5,   7,  11,  13,
    29,  29,  29,  29,  29,
    71,  71,  71,  71,  71,
   173, 173, 173, 173, 173,
   409, 409, 409, 409, 409 };

  static int seed_vec[N_MAX] = {
        0,     1,     2,     3,     4,
        5,     6,     7,     8,     0,
        1,     2,     3,     4,     5,
        6,     7,     8,     0,     1,
        2,     3,     4,     5,     6,
        7,     8,    10,    10,    10,
       10,    10,    10,    10,   100,
      100,   100,   100,   100,   100,
      100,  1000,  1000,  1000,  1000,
     1000,  1000,  1000, 10000, 10000,
    10000, 10000, 10000, 10000, 10000,
     1000,  1001,  1002,  1003,  1004,
     1000,  1001,  1002,  1003,  1004,
     1000,  1001,  1002,  1003,  1004,
     1000,  1001,  1002,  1003,  1004 };

  static double value_vec[N_MAX] = {
     0.0000000000000000E+00,
     0.5000000000000000E+00,
     0.2500000000000000E+00,
     0.7500000000000000E+00,
     0.1250000000000000E+00,
     0.6250000000000000E+00,
     0.3750000000000000E+00,
     0.8750000000000000E+00,
     0.0625000000000000E+00,
     0.0000000000000000E+00,
     0.3333333333333333E+00,
     0.6666666666666666E+00,
     0.1111111111111111E+00,
     0.4444444444444444E+00,
     0.7777777777777777E+00,
     0.2222222222222222E+00,
     0.5555555555555556E+00,
     0.8888888888888888E+00,
     0.0000000000000000E+00,
     0.2500000000000000E+00,
     0.5000000000000000E+00,
     0.7500000000000000E+00,
     0.0625000000000000E+00,
     0.3125000000000000E+00,
     0.5625000000000000E+00,
     0.8125000000000000E+00,
     0.1250000000000000E+00,
     0.3125000000000000E+00,
     0.3703703703703703E+00,
     0.6250000000000000E+00,
     0.0800000000000000E+00,
     0.4489795918367347E+00,
     0.9090909090909092E+00,
     0.7692307692307693E+00,
     0.1484375000000000E+00,
     0.4115226337448559E+00,
     0.0976562500000000E+00,
     0.0320000000000000E+00,
     0.2915451895043731E+00,
     0.1652892561983471E+00,
     0.7337278106508875E+00,
     0.0927734375000000E+00,
     0.3475080018289895E+00,
     0.1708984375000000E+00,
     0.0051200000000000E+00,
     0.9162848812994586E+00,
     0.9316303531179565E+00,
     0.9904415111515704E+00,
     0.0347290039062500E+00,
     0.3861200020322105E+00,
     0.0189208984375000E+00,
     0.0005120000000000E+00,
     0.5749985125245433E+00,
     0.1529950140017758E+00,
     0.2459297643639929E+00,
     0.4887449259912255E+00,
     0.5232276846119153E+00,
     0.5577104432326049E+00,
     0.5921932018532945E+00,
     0.6266759604739842E+00,
     0.0872842689942472E+00,
     0.1013687760365007E+00,
     0.1154532830787542E+00,
     0.1295377901210077E+00,
     0.1436222971632613E+00,
     0.7805138828560928E+00,
     0.7862942296769020E+00,
     0.7920745764977113E+00,
     0.7978549233185205E+00,
     0.8036352701393298E+00,
     0.4449997309915651E+00,
     0.4474447187666262E+00,
     0.4498897065416874E+00,
     0.4523346943167484E+00,
     0.4547796820918096E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    base = 0;
    seed = 0;
    value = 0.0;
  }
  else
  {
    base = base_vec[n_data-1];
    seed = seed_vec[n_data-1];
    value = value_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void viscosity_values ( int &n_data, double &tc, double &p, double &eta )

//****************************************************************************80
//
//  Purpose:
//
//    VISCOSITY_VALUES returns some values of the viscosity function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 February 2002
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Lester Haar, John Gallagher and George Kell,
//    NBS/NRC Steam Tables:
//    Thermodynamic and Transport Properties and Computer Programs
//    for Vapor and Liquid States of Water in SI Units,
//    Hemisphere Publishing Corporation, Washington, 1984,
//    TJ270.H3, page 263.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &TC, the temperature, in degrees Celsius.
//
//    Output, double &P, the pressure, in bar.
//
//    Output, double &ETA, the viscosity, in MegaPascal seconds.
//
{
# define N_MAX 34

  static double eta_vec[N_MAX] = {
     1792.0E+00,
     1791.0E+00,
     1790.0E+00,
     1786.0E+00,
     1780.0E+00,
     1775.0E+00,
     1769.0E+00,
     1764.0E+00,
     1759.0E+00,
     1754.0E+00,
     1749.0E+00,
     1744.0E+00,
     1739.0E+00,
     1735.0E+00,
     1731.0E+00,
     1722.0E+00,
     1714.0E+00,
     1707.0E+00,
     1700.0E+00,
     1694.0E+00,
     1687.0E+00,
     1682.0E+00,
     1676.0E+00,
     1667.0E+00,
     1659.0E+00,
     1653.0E+00,
      890.8E+00,
      547.1E+00,
      378.4E+00,
      12.28E+00,
      16.18E+00,
      24.45E+00,
      32.61E+00,
      40.38E+00 };

  static double p_vec[N_MAX] = {
        1.0E+00,
        5.0E+00,
       10.0E+00,
       25.0E+00,
       50.0E+00,
       75.0E+00,
      100.0E+00,
      125.0E+00,
      150.0E+00,
      175.0E+00,
      200.0E+00,
      225.0E+00,
      250.0E+00,
      275.0E+00,
      300.0E+00,
      350.0E+00,
      400.0E+00,
      450.0E+00,
      500.0E+00,
      550.0E+00,
      600.0E+00,
      650.0E+00,
      700.0E+00,
      800.0E+00,
      900.0E+00,
     1000.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00,
        1.0E+00 };

  static double tc_vec[N_MAX] = {
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
       0.0E+00,
      25.0E+00,
      50.0E+00,
      75.0E+00,
     100.0E+00,
     200.0E+00,
     400.0E+00,
     600.0E+00,
     800.0E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    tc = 0.0;
    p = 0.0;
    eta = 0.0;
  }
  else
  {
    tc = tc_vec[n_data-1];
    p = p_vec[n_data-1];
    eta = eta_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void von_mises_cdf_values ( int &n_data, double &a, double &b, double &x,
  double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    VON_MISES_CDF_VALUES returns some values of the von Mises CDF.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 December 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Kanti Mardia, Peter Jupp,
//    Directional Statistics,
//    Wiley, 2000, QA276.M335
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &A, &B, the parameters of the function.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 23

  static double a_vec[N_MAX] = {
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
     0.1E+01,
    -0.2E+01,
    -0.1E+01,
     0.0E+01,
     0.1E+01,
     0.2E+01,
     0.3E+01,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00,
     0.0E+00 };

  static double b_vec[N_MAX] = {
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.1E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.2E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.3E+01,
      0.0E+00,
      0.1E+01,
      0.2E+01,
      0.3E+01,
      0.4E+01,
      0.5E+01 };

  static double fx_vec[N_MAX] = {
     0.2535089956281180E-01,
     0.1097539041177346E+00,
     0.5000000000000000E+00,
     0.8043381312498558E+00,
     0.9417460124555197E+00,
     0.5000000000000000E+00,
     0.6018204118446155E+00,
     0.6959356933122230E+00,
     0.7765935901304593E+00,
     0.8410725934916615E+00,
     0.8895777369550366E+00,
     0.9960322705517925E+00,
     0.9404336090170247E+00,
     0.5000000000000000E+00,
     0.5956639098297530E-01,
     0.3967729448207649E-02,
     0.2321953958111930E-03,
     0.6250000000000000E+00,
     0.7438406999109122E+00,
     0.8369224904294019E+00,
     0.8941711407897124E+00,
     0.9291058600568743E+00,
     0.9514289900655436E+00 };

  static double x_vec[N_MAX] = {
     -0.2617993977991494E+01,
     -0.1570796326794897E+01,
      0.0000000000000000E+00,
      0.1047197551196598E+01,
      0.2094395102393195E+01,
      0.1000000000000000E+01,
      0.1200000000000000E+01,
      0.1400000000000000E+01,
      0.1600000000000000E+01,
      0.1800000000000000E+01,
      0.2000000000000000E+01,
      0.0000000000000000E+00,
      0.0000000000000000E+00,
      0.0000000000000000E+00,
      0.0000000000000000E+00,
      0.0000000000000000E+00,
      0.0000000000000000E+00,
      0.7853981633974483E+00,
      0.7853981633974483E+00,
      0.7853981633974483E+00,
      0.7853981633974483E+00,
      0.7853981633974483E+00,
      0.7853981633974483E+00 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    a = 0.0;
    b = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    a = a_vec[n_data-1];
    b = b_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void weekday_values ( int &n_data, int &y, int &m, int &d, int &w )

//****************************************************************************80
//
//  Purpose:
//
//    WEEKDAY_VALUES returns the day of the week for various dates.
//
//  Discussion:
//
//    The CE or Common Era calendar is used, under the
//    hybrid Julian/Gregorian Calendar, with a transition from Julian
//    to Gregorian.  The day after 04 October 1582 was 15 October 1582.
//
//    The year before 1 AD or CE is 1 BC or BCE.  In this data set,
//    years BC/BCE are indicated by a negative year value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Edward Reingold, Nachum Dershowitz,
//    Calendrical Calculations: The Millennium Edition,
//    Cambridge University Press, 2001,
//    ISBN: 0 521 77752 6
//    LC: CE12.R45.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0
//    before the first call.  On each call, the routine increments N_DATA by 1,
//    and returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &Y, &M, &D, the Common Era date.
//
//    Output, int &W, the day of the week.  Sunday = 1.
//
{
# define N_MAX 34

  static int d_vec[N_MAX] = {
    30,
     8,
    26,
     3,
     7,
    18,
     7,
    19,
    14,
    18,
    16,
     3,
    26,
    20,
     4,
    25,
    31,
     9,
    24,
    10,
    30,
    24,
    19,
     2,
    27,
    19,
    25,
    29,
    19,
     7,
    17,
    25,
    10,
    18 };
  static int m_vec[N_MAX] = {
     7,
    12,
     9,
    10,
     1,
     5,
    11,
     4,
    10,
     5,
     3,
     3,
     3,
     4,
     6,
     1,
     3,
     9,
     2,
     6,
     6,
     7,
     6,
     8,
     3,
     4,
     8,
     9,
     4,
    10,
     3,
     2,
    11,
     7 };
  static int w_vec[N_MAX] = {
    1,
    4,
    4,
    1,
    4,
    2,
    7,
    1,
    7,
    1,
    6,
    7,
    6,
    1,
    1,
    4,
    7,
    7,
    7,
    4,
    1,
    6,
    1,
    2,
    4,
    1,
    1,
    2,
    2,
    5,
    3,
    1,
    4,
    1 };
  static int y_vec[N_MAX] = {
    - 587,
    - 169,
       70,
      135,
      470,
      576,
      694,
     1013,
     1066,
     1096,
     1190,
     1240,
     1288,
     1298,
     1391,
     1436,
     1492,
     1553,
     1560,
     1648,
     1680,
     1716,
     1768,
     1819,
     1839,
     1903,
     1929,
     1941,
     1943,
     1943,
     1992,
     1996,
     2038,
     2094 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    y = 0;
    m = 0;
    d = 0;
    w = 0;
  }
  else
  {
    y = y_vec[n_data-1];
    m = m_vec[n_data-1];
    d = d_vec[n_data-1];
    w = w_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void weibull_cdf_values ( int &n_data, double &alpha, double &beta,
  double &x, double &fx )

//****************************************************************************80
//
//  Purpose:
//
//    WEIBULL_CDF_VALUES returns some values of the Weibull CDF.
//
//  Discussion:
//
//    In Mathematica, the function can be evaluated by:
//
//      Needs["Statistics`ContinuousDistributions`"]
//      dist = WeibullDistribution [ alpha, beta ]
//      CDF [ dist, x ]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, double &ALPHA, the first parameter of the distribution.
//
//    Output, double &BETA, the second parameter of the distribution.
//
//    Output, double &X, the argument of the function.
//
//    Output, double &FX, the value of the function.
//
{
# define N_MAX 12

  static double alpha_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01 };

  static double beta_vec[N_MAX] = {
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.5000000000000000E+00,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.5000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01 };

  static double fx_vec[N_MAX] = {
     0.8646647167633873E+00,
     0.9816843611112658E+00,
     0.9975212478233336E+00,
     0.9996645373720975E+00,
     0.6321205588285577E+00,
     0.4865828809674080E+00,
     0.3934693402873666E+00,
     0.3296799539643607E+00,
     0.8946007754381357E+00,
     0.9657818816883340E+00,
     0.9936702845725143E+00,
     0.9994964109502630E+00 };

  static double x_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.4000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.2000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01,
     0.3000000000000000E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    alpha = 0.0;
    beta = 0.0;
    x = 0.0;
    fx = 0.0;
  }
  else
  {
    alpha = alpha_vec[n_data-1];
    beta = beta_vec[n_data-1];
    x = x_vec[n_data-1];
    fx = fx_vec[n_data-1];
  }

  return;
# undef N_MAX
}
//****************************************************************************80

void zeta_values ( int &n_data, int &n, double &zeta )

//****************************************************************************80
//
//  Purpose:
//
//    ZETA_VALUES returns some values of the Riemann Zeta function.
//
//  Discussion:
//
//    ZETA(N) = sum ( 1 <= I < +oo ) 1 / I^N
//
//    In Mathematica, the function can be evaluated by:
//
//      Zeta[n]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    12 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz, Irene Stegun,
//    Handbook of Mathematical Functions,
//    National Bureau of Standards, 1964,
//    ISBN: 0-486-61272-4,
//    LC: QA47.A34.
//
//    Stephen Wolfram,
//    The Mathematica Book,
//    Fourth Edition,
//    Cambridge University Press, 1999,
//    ISBN: 0-521-64314-7,
//    LC: QA76.95.W65.
//
//  Parameters:
//
//    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
//    first call.  On each call, the routine increments N_DATA by 1, and
//    returns the corresponding data; when there is no more data, the
//    output value of N_DATA will be 0 again.
//
//    Output, int &N, the argument of the Zeta function.
//
//    Output, double &ZETA, the value of the Zeta function.
//
{
# define N_MAX 15

  static int n_vec[N_MAX] = {
     2,
     3,
     4,
     5,
     6,
     7,
     8,
     9,
    10,
    11,
    12,
    16,
    20,
    30,
    40 };

  static double zeta_vec[N_MAX] = {
     0.164493406684822643647E+01,
     0.120205690315959428540E+01,
     0.108232323371113819152E+01,
     0.103692775514336992633E+01,
     0.101734306198444913971E+01,
     0.100834927738192282684E+01,
     0.100407735619794433939E+01,
     0.100200839292608221442E+01,
     0.100099457512781808534E+01,
     0.100049418860411946456E+01,
     0.100024608655330804830E+01,
     0.100001528225940865187E+01,
     0.100000095396203387280E+01,
     0.100000000093132743242E+01,
     0.100000000000090949478E+01 };

  if ( n_data < 0 )
  {
    n_data = 0;
  }

  n_data = n_data + 1;

  if ( N_MAX < n_data )
  {
    n_data = 0;
    n = 0;
    zeta = 0.0;
  }
  else
  {
    n = n_vec[n_data-1];
    zeta = zeta_vec[n_data-1];
  }

  return;
# undef N_MAX
}
