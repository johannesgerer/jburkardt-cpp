# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <stdint.h>

using namespace std;

# include "ziggurat_inline.hpp"

static float fe[256];
static float fn[128];
static int32_t hz;
static uint32_t iz;
static uint32_t jcong = 234567891;
static uint32_t jsr = 123456789;
static uint32_t jz;
static uint32_t ke[256];
static uint32_t kn[128];
static uint32_t w = 345678912;
static float we[256];
static float wn[128];
static uint32_t z = 456789123;
//
//  The original SHR3 random number generator was replaced by
//  KISS, a combination of MWC, CONG and SHR3 as suggested in
//  the Leong reference.
//
//  Thanks to Dirk Eddelbuettel, 04 October 2013.
//
# define znew ( z = 36969 * ( z & 65535 ) + ( z >> 16 ) )
# define wnew ( w = 18000 * ( w & 65535 ) + ( w >> 16 ) )
# define MWC ( ( znew << 16 ) + wnew )
# define SHR3 ( jz = jsr, jsr ^= ( jsr << 13 ), jsr ^= ( jsr >> 17 ), jsr ^= ( jsr << 5 ), jz + jsr )
# define CONG ( jcong = 69069 * jcong + 1234567 )
# define KISS ( ( MWC ^ CONG ) + SHR3 )

# define UNI ( 0.5 + ( signed ) KISS * 0.2328306e-09 )
# define IUNI KISS
# define RNOR ( hz = KISS, iz = hz & 127, ( fabs ( hz ) < kn[iz] ) ? hz * wn[iz] : nfix() )
# define REXP ( jz = KISS, iz = jz & 255, (        jz   < ke[iz] ) ? jz * we[iz] : efix() )

//****************************************************************************80

uint32_t cong_seeded ( uint32_t &jcong )

//****************************************************************************80
//
//  Purpose:
//
//    CONG_SEEDED evaluates the CONG generator for uint32_t's.
//
//  Discussion:
//
//    This function requires the user to supply the seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    George Marsaglia, Wai Wan Tsang.
//    Modifications by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JCONG, the seed.
//
//    Output, uint32_t CONG_SEEDED, the randomly chosen value.
//
{
  jcong = 69069 * ( jcong ) + 1234567;

  return jcong;
}
//****************************************************************************80

uint32_t cong_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    CONG_VALUE evaluates the CONG generator for uint32_t's.
//
//  Discussion:
//
//    This function uses an internally managed seed, which can be viewed by
//    ZIGGET or set by ZIGSET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    George Marsaglia, Wai Wan Tsang.
//    Modifications by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t CONG_VALUE, the randomly chosen value.
//
{
  return CONG;
}
//****************************************************************************80

double cpu_time ( )

//****************************************************************************80
//
//  Purpose:
//
//    CPU_TIME returns the current reading on the CPU clock.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 June 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double CPU_TIME, the current reading of the CPU clock, in seconds.
//
{
  double value;

  value = ( double ) clock ( ) / ( double ) CLOCKS_PER_SEC;

  return value;
}
//****************************************************************************80

float efix ( )

//****************************************************************************80
//
//  Purpose:
//
//    EFIX generates variates when rejection occurs in the exponential code. 
//
//  Discussion:
//
//    This routine is NOT designed for direct calls by the user!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2008
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, float EFIX, an exponential variate.
//
{ 
  float x;

  for ( ; ; ) 
  {
//
//  IZ = 0.
//
    if ( iz == 0 )
    {
      return ( 7.69711 - log ( UNI ) );
    }

    x = jz * we[iz];
    if ( fe[iz] + UNI * ( fe[iz-1] - fe[iz] ) < exp ( - x ) ) 
    {
      return x;
    }
// 
//  Initiate, try to exit the loop.
//
    jz = KISS;
    iz = ( jz & 255 );
    if ( jz < ke[iz] ) 
    {
      return ( ( float ) ( jz * we[iz] ) );
    }
  }
}
//****************************************************************************80

uint32_t kiss_seeded ( uint32_t &jcong, uint32_t &jsr, uint32_t &w, uint32_t &z )

//****************************************************************************80
//
//  Purpose:
//
//    KISS_SEEDED evaluates the KISS random number generator.
//
//  Discussion:
//
//    This function requires the user to supply the seeds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JCONG, uint32_t &JSR, uint32_t &W, uint32_t &Z, 
//    the seeds, which are updated on each call.
//
//    Output, uint32_t KISS_SEEDED, the new value.
//
{
  uint32_t value;

  value = ( mwc_seeded ( w, z ) ^ cong_seeded ( jcong ) ) + shr3_seeded ( jsr );

  return value;
}
//****************************************************************************80

uint32_t kiss_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    KISS_VALUE evaluates the KISS generator for uint32_t's.
//
//  Discussion:
//
//    This provides a functional interface to the inline KISS generator.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2013
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t KISS_VALUE, the value of the KISS generator.
//
{
  return KISS;
}
//*****************************************************************************/

uint32_t mwc_seeded ( uint32_t &w, uint32_t &z )

//*****************************************************************************/
//
//  Purpose:
//
//    MWC_SEEDED evaluates the MWC multiply-with-carry random number generator.
//
//  Discussion:
//
//    This function requires the user to supply the seeds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &W, uint32_t &Z, the seeds, which are updated 
//    on each call.
//
//    Output, uint32_t MWC_SEEDED, the new value.
//
{
  uint32_t value;

  z = 36969 * ( z & 65535 ) + ( z >> 16 );
  w = 18000 * ( w & 65535 ) + ( w >> 16 );

  value = ( z << 16 ) + w;

  return value;
}
//****************************************************************************80

uint32_t mwc_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    MWC_VALUE evaluates the MWC generator for uint32_t's.
//
//  Discussion:
//
//    This function uses internally managed seeds, which can be viewed by
//    ZIGGET or set by ZIGSET.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    George Marsaglia, Wai Wan Tsang.
//    Modifications by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t MWC_VALUE, the randomly chosen value.
//
{
  return MWC;
}
//****************************************************************************80

float nfix ( )

//****************************************************************************80
//
//  Purpose:
// 
//    NFIX generates variates when rejection occurs in the normal code.
//
//  Discussion:
//
//    This routine is NOT designed for direct calls by the user!
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, float NFIX, a normal variate.
//
{
  const float r = 3.442620;
  static float x;
  static float y;

  for ( ; ; )
  {
//
//  IZ = 0 handles the base strip.
//
    x = ( float ) ( hz * wn[iz] );
    if ( iz == 0 )
    { 
      do
      {
        x = - log ( UNI ) * 0.2904764; 
        y = - log ( UNI );
      }
      while ( y + y < x * x );

      return ( 0 < hz ) ? r + x : - r - x;
    }
// 
//  0 < IZ, handle the wedges of other strips.
//
    if ( fn[iz] + UNI * ( fn[iz-1] - fn[iz] ) < exp ( - 0.5 * x * x ) ) 
    {
      return x;
    }
// 
//  Initiate, try to exit the loop.
//
    hz = SHR3;
    iz = ( hz & 127 );
    if ( fabs ( hz ) < kn[iz] )
    {
      return ( ( float ) ( hz * wn[iz] ) );
    }
  }
}
//****************************************************************************80

void r4_exp_setup ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EXP_SETUP sets data needed by R4_EXP.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Global, uint32_t KE[256], data needed by R4_EXP.
//
//    Global, float FE[256], WE[256], data needed by R4_EXP.
//
{
  double de = 7.697117470131487;
  int i;
  const double m2 = 4294967296.0;
  double q;
  double te = 7.697117470131487;
  double ve = 3.949659822581572e-03;

  q = ve / exp ( - de );

  ke[0] = ( uint32_t ) ( ( de / q ) * m2 );
  ke[1] = 0;

  we[0] = ( float ) q / m2;
  we[255] = ( float ) de / m2;

  fe[0] = 1.0;
  fe[255] = ( float ) exp ( - de );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    ke[i+1] = ( uint32_t ) ( ( de / te ) * m2 );
    te = de;
    fe[i] = ( float ) exp ( - de );
    we[i] = ( float ) ( de / m2 );
  }
  return;
}
//****************************************************************************80

float r4_exp_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_EXP_VALUE returns an exponentially distributed float.
//
//  Discussion:
//
//    The value returned is generated from a distribution on [0,+00) with density 
//    exp(-x).
//
//    The underlying algorithm is the ziggurat method.
//
//    Before calling this function, the user must call ZIGSET, supplying
//    a nonzero seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2008
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, float R4_EXP_VALUE, an exponentially distributed random value.
//
{
  return REXP;
}
//****************************************************************************80

void r4_nor_setup ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR_SETUP sets data needed by R4_NOR.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Global, uint32_t KN[128], data needed by R4_NOR.
//
//    Global, float FN[128], WN[128], data needed by R4_NOR.
//
{
  double dn = 3.442619855899;
  int i;
  const double m1 = 2147483648.0;
  double q;
  double tn = 3.442619855899;
  double vn = 9.91256303526217e-03;
//
//  Set up the tables for the normal random number generator.
//
  q = vn / exp ( - 0.5 * dn * dn );

  kn[0] = ( uint32_t ) ( ( dn / q ) * m1 );
  kn[1] = 0;

  wn[0] = ( float ) ( q / m1 );
  wn[127] = ( float ) ( dn / m1 );

  fn[0] = 1.0;
  fn[127] = ( float ) exp ( - 0.5 * dn * dn );

  for ( i = 126; 1 <= i; i-- )
  {
    dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
    kn[i+1] = ( uint32_t ) ( ( dn / tn ) * m1 );
    tn = dn;
    fn[i] = ( float ) exp ( - 0.5 * dn * dn );
    wn[i] = ( float ) ( dn / m1 );
  }
  return;
}
//****************************************************************************80

float r4_nor_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_NOR_VALUE returns a normally distributed float.
//
//  Discussion:
//
//    The value returned is generated from a distribution with mean 0 and 
//    variance 1.
//
//    The underlying algorithm is the ziggurat method.
//
//    Before calling this function, the user must call ZIGSET, supplying
//    a nonzero seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2008
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, float R4_NOR_VALUE, a normally distributed random value.
//
{
  return RNOR;
}
//****************************************************************************80

float r4_uni_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    R4_UNI_VALUE returns a uniformly distributed float.
//
//  Discussion:
//
//    Before calling this function, the user may call ZIGSET, supplying
//    a nonzero seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 May 2008
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, float R4_UNI_VALUE, a uniformly distributed random value in
//    the range [0,1].
//
{
  return UNI;
}
//****************************************************************************80

uint32_t shr3_seeded ( uint32_t &jsr )

//****************************************************************************80
//
//  Purpose:
//
//    SHR3_SEEDED evaluates the SHR3 generator for unsigned 32 bit integers.
//
//  Discussion:
//
//    This function requires the user to supply the seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2013
//
//  Author:
//
//    George Marsaglia, Wai Wan Tsang.
//    Modifications by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input/output, uint32_t &JSR, the seed, which is updated 
//    on each call.
//
//    Output, uint32_t SHR3_SEEDED, the value of the SHR3 generator.
//
{
  uint32_t value;

  value = jsr;

  jsr = ( jsr ^ ( jsr <<   13 ) );
  jsr = ( jsr ^ ( jsr >>   17 ) );
  jsr = ( jsr ^ ( jsr <<    5 ) );

  value = value + jsr;

  return value;
}
//****************************************************************************80

uint32_t shr3_value ( )

//****************************************************************************80
//
//  Purpose:
//
//    SHR3_VALUE evaluates the SHR3 generator for unsigned 32 bit integers.
//
//  Discussion:
//
//    Before calling this function, the user may call ZIGSET, supplying
//    a nonzero seed.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2008
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t SHR3_VALUE, the value of the SHR3 generator.
//
{
  return SHR3;
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
//    31 May 2001 09:45:54 AM
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

void zigget ( uint32_t &jsr_value, uint32_t &jcong_value,
  uint32_t &w_value, uint32_t &z_value )

//****************************************************************************80
//
//  Purpose:
//
//    ZIGGET gets the seeds.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2013
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Output, uint32_t &JSR_VALUE, the seed for the SHR3 generator.
//
//    Output, uint32_t &JCONG_VALUE, the seed for the congruential generator.
//
//    Output, uint32_t &W_VALUE, the seed for the first MWC generator.
//
//    Output, uint32_t &Z_VALUE, the seed for the second MWC generator.
//
{
  jsr_value = jsr;
  jcong_value = jcong;
  w_value = w;
  z_value = z;

  return;
}
//****************************************************************************80

void zigset ( uint32_t jsr_value, uint32_t jcong_value, uint32_t w_value, 
  uint32_t z_value )

//****************************************************************************80
//
//  Purpose:
//
//    ZIGSET sets the seeds and creates the tables for the Ziggurat method.
//
//  Discussion:
//
//    ZIGSET must be called before the exponential and normal random number
//    generators are called.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 October 2013
//
//  Author:
//
//    Original C version by George Marsaglia, Wai Wan Tsang.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Philip Leong, Guanglie Zhang, Dong-U Lee, Wayne Luk, John Villasenor,
//    A comment on the implementation of the ziggurat method,
//    Journal of Statistical Software,
//    Volume 12, Number 7, February 2005.
//
//    George Marsaglia, Wai Wan Tsang,
//    The Ziggurat Method for Generating Random Variables,
//    Journal of Statistical Software,
//    Volume 5, Number 8, October 2000, seven pages.
//
//  Parameters:
//
//    Input, uint32_t JSR_VALUE, the seed for the SHR3 generator.
//
//    Input, uint32_t JCONG_VALUE, the seed for the congruential generator.
//
//    Input, uint32_t W_VALUE, the seed for the first MWC generator.
//
//    Input, uint32_t Z_VALUE, the seed for the second MWC generator.
//
{
  jsr = jsr_value;
  jcong = jcong_value;
  w = w_value;
  z = z_value;
//
//  Set up the tables for the normal random number generator.
//
  r4_nor_setup ( );
//
//  Set up tables for the exponential random number generator.
//
  r4_exp_setup ( );

  return;
}
