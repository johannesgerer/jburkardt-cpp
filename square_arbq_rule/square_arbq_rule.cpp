# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "square_arbq_rule.hpp"

//****************************************************************************80

double *lege2eva ( int degree, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGE2EVA evaluates orthogonal polynomials on the symmetric square.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the maximum degree of the polynomials.
//
//    Input, double Z[2], the evaluation point.
//
//    Output, double LEGE2EVA[(DEGREE+1)*(DEGREE+2)/2], the orthogonal
//    polynomials evaluated at Z.
//
{
  double *f1;
  double *f2;
  int kk;
  int m;
  int n;
  int npols;
  double *pols;
  double scale;

  npols = ( ( degree + 1 ) * ( degree + 2 ) ) / 2;
  pols = new double[npols];

  f1 = llegepols1 ( degree, z[0] );
  f2 = llegepols1 ( degree, z[1] );

  kk = 0;
  for ( m = 0; m <= degree; m++ )
  {
    for ( n = 0; n <= m; n++ )
    {
      pols[kk] = f1[m-n] * f2[n];
      scale = ( double ) ( ( 1 + 2 * n ) * ( 1 + 2 * ( m - n ) ) );
      scale = 0.5 * sqrt ( scale );
      pols[kk] = pols[kk] * scale;
      kk = kk + 1;
    }
  }

  delete [] f1;
  delete [] f2;

  return pols;
}
//****************************************************************************80

double *llegepols1 ( int degree, double x )

//****************************************************************************80
//
//  Purpose:
//
//    LLEGEPOLS1 evaluates orthogonal polynomials on the symmetric interval.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the maximum degree.
//
//    Input, double X, the evaluation point.
//
//    Output, double LLEGEPOLS1[DEGREE+1], the orthogonal
//    polynomials evaluated at X.
//
{
  int k;
  double pk;
  double pkm1;
  double pkp1;
  double *pols;

  pols = new double[degree+1];
  pkp1 = 1.0;
  pols[0] = pkp1;

  if ( degree == 0 )
  {
    return pols;
  }

  pk = pkp1;
  pkp1 = x;
  pols[1] = pkp1;

  if ( degree == 1 )
  {
    return pols;
  }

  for ( k = 1; k <= degree - 1; k++ )
  {
    pkm1 = pk;
    pk = pkp1;
    pkp1 = ( ( double ) ( 2 * k + 1 ) * x * pk 
           - ( double ) (     k     ) * pkm1 ) 
           / ( double ) (     k + 1 );

    pols[k+1] = pkp1;
  }

  return pols;
}
//****************************************************************************80

void r8mat_row_copy ( int m, int n, int i, double v[], double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_ROW_COPY copies a vector into a row of an R8MAT.
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
//    30 June 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the order of the matrix.
//
//    Input, int I, the index of the row.
//    0 <= I <= M-1.
//
//    Input, double V[N], the row to be copied.
//
//    Input/output, double A[M*N], the matrix into which
//    the row is to be copied.
//
{
  int j;

  for ( j = 0; j < n; j++ )
  {
    a[i+j*m] = v[j];
  }
  return;
}
//****************************************************************************80

void r8vec_copy ( int n, double a1[], double a2[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_COPY copies an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], the vector to be copied.
//
//    Input, double A2[N], the copy of A1.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
//****************************************************************************80

double r8vec_sum ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 August 2008
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a[i];
  }

  return value;
}
//****************************************************************************80

void rule01 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE01 returns the rule of degree 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[1] = {
       0.00000000000000000E+00 };
  double ys[1] = {
       0.00000000000000000E+00 };
  double ws[1] = {
       0.2828427124746189E+01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule02 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE02 returns the rule of degree 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.6700319885140564E+00,0.6424528854269665E+00, 
    -.5079273596590297E+00 };
  static double ys[] = {
    -.8727274074699508E+00,0.8751842913892396E+00, 
    -.8014238374817481E-02 };
  static double ws[] = {
    0.6106555690526828E+00,0.6235399890121793E+00, 
    0.1594231566681328E+01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule03 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE03 returns the rule of degree 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.5773502691896256E+00,0.5773502691896260E+00, 
    0.5773502691896256E+00,-.5773502691896260E+00 };
  static double ys[] = {
    -.5773502691896260E+00,-.5773502691896256E+00, 
    0.5773502691896260E+00,0.5773502691896256E+00 };
  static double ws[] = {
    0.7071067811865476E+00,0.7071067811865476E+00, 
    0.7071067811865476E+00,0.7071067811865476E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule04 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE04 returns the rule of degree 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.8434867168489557E+00,0.7481765370372371E+00, 
    -.4061851500656107E+00,-.4581090172172534E+00, 
    0.1816993641713940E+00,-.9077196977637252E+00 };
  static double ys[] = {
    0.7332250861965538E+00,-.6280294280975105E+00, 
    -.7973798546121016E+00,0.8743017248509551E+00, 
    0.1628466016041256E+00,0.8506794801388022E-02 };
  static double ws[] = {
    0.2754640609309160E+00,0.4372667066134153E+00, 
    0.4966805368802857E+00,0.3707670373943532E+00, 
    0.8675752571634261E+00,0.3806735257637942E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule05 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE05 returns the rule of degree 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.1775868202077551E-01,-.1775868202077539E-01, 
    0.7788710544649639E+00,-.7788710544649639E+00, 
    -.7703781288541645E+00,0.7703781288541645E+00, 
    -.7490353914168658E-33 };
  static double ys[] = {
    -.9659285494001192E+00,0.9659285494001192E+00, 
    -.5715708301251639E+00,0.5715708301251639E+00, 
    -.5829672991828014E+00,0.5829672991828014E+00, 
    0.1356144833394667E-33 };
  static double ws[] = {
    0.2246199725165690E+00,0.2246199725165690E+00, 
    0.3901817339168917E+00,0.3901817339168917E+00, 
    0.3953508381187504E+00,0.3953508381187504E+00, 
    0.8081220356417684E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule06 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE06 returns the rule of degree 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.7347550761673839E+00,0.8662152988634034E+00, 
    0.1596873653424614E+00,-.8905137714296896E+00, 
    -.1469707846748791E+00,0.9240009259977663E+00, 
    -.8463324986375500E+00,-.4086308482879689E+00, 
    0.5175294652720337E+00,0.4801002492857063E+00 };
  static double ys[] = {
    0.8933891941643415E+00,-.7037359670513631E+00, 
    -.9085856749287847E+00,0.1644347368502312E+00, 
    0.5352177835541986E+00,0.4879643750888035E+00, 
    -.8394767448218339E+00,-.4262330870004397E+00, 
    0.9176357850707917E+00,-.1009764561823168E+00 };
  static double ws[] = {
    0.1541850714382379E+00,0.1900556513689156E+00, 
    0.2246942645703077E+00,0.2465847648329768E+00, 
    0.5062382287542438E+00,0.1829226437278864E+00, 
    0.1373586623279704E+00,0.4754388545735908E+00, 
    0.1856913242244974E+00,0.5252576589275637E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule07 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE07 returns the rule of degree 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.4595981103653579E-16,0.9258200997725515E+00, 
    0.6742045114073804E-16,-.9258200997725515E+00, 
    -.3805544332083157E+00,0.3805544332083157E+00, 
    0.3805544332083157E+00,-.3805544332083157E+00, 
    -.8059797829185990E+00,0.8059797829185988E+00, 
    0.8059797829185990E+00,-.8059797829185988E+00 };
  static double ys[] = {
    -.9258200997725515E+00,-.1073032005210112E-16, 
    0.9258200997725515E+00,0.1241105822293750E-15, 
    -.3805544332083157E+00,-.3805544332083157E+00, 
    0.3805544332083157E+00,0.3805544332083157E+00, 
    -.8059797829185988E+00,-.8059797829185990E+00, 
    0.8059797829185988E+00,0.8059797829185990E+00 };
  static double ws[] = {
    0.1711023816204485E+00,0.1711023816204485E+00, 
    0.1711023816204485E+00,0.1711023816204485E+00, 
    0.3681147816131979E+00,0.3681147816131979E+00, 
    0.3681147816131979E+00,0.3681147816131979E+00, 
    0.1678896179529011E+00,0.1678896179529011E+00, 
    0.1678896179529011E+00,0.1678896179529011E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule08 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE08 returns the rule of degree 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.2272218649369121E+00,0.2786782798124801E+00, 
    0.9215721988395638E+00,-.5229427015551803E+00, 
    0.8309170589376613E+00,-.6080254018462903E+00, 
    -.9822549066167084E+00,0.4959470731361600E-01, 
    0.5910013957537859E+00,0.3626589212754838E+00, 
    -.9369162594801185E+00,-.8850131220663160E+00, 
    -.1934658240272289E+00,0.5772453681919104E+00, 
    0.9213070164035271E+00,-.7176037958340967E+00 };
  static double ys[] = {
    0.8703146041404044E+00,0.9856262640199153E+00, 
    0.2224095500358621E+00,-.9282264259882677E+00, 
    0.8435111761265234E+00,0.5825946042673711E+00, 
    -.8211266831948021E+00,-.6917239446781449E+00, 
    -.2614406969784849E+00,0.5198121135620160E+00, 
    0.2153771996329335E+00,0.9090384216207131E+00, 
    0.3526321874643216E-01,-.9622555555961493E+00, 
    -.7082682674817122E+00,-.4130619139730907E+00 };
  static double ws[] = {
    0.1444235515837947E+00,0.5206905878850610E-01, 
    0.1365819925705312E+00,0.1136963049256808E+00, 
    0.1156201396846171E+00,0.2194056396025883E+00, 
    0.4570142629159132E-01,0.3040158377300561E+00, 
    0.3227772111095287E+00,0.3341175763908440E+00, 
    0.1202823186503543E+00,0.6155232134515501E-01, 
    0.4037250536437860E+00,0.8510021531985533E-01, 
    0.1026971066172272E+00,0.2666613704920739E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule09 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE09 returns the rule of degree 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.7502770999789001E+00,-.7502770999789001E+00, 
    -.9279616459595694E+00,0.9279616459595694E+00, 
    0.6306801197316682E+00,-.6306801197316682E+00, 
    0.9688499663619775E+00,-.9688499663619775E+00, 
    -.7620832819261721E-01,0.7620832819261719E-01, 
    0.8526157293336627E+00,-.8526157293336627E+00, 
    0.4533398211356476E+00,-.5237358202144297E+00, 
    0.5237358202144297E+00,-.4533398211356476E+00, 
    -.7154960467453349E-17 };
  static double ys[] = {
    -.9279616459595700E+00,0.9279616459595700E+00, 
    -.7502770999789009E+00,0.7502770999789010E+00, 
    0.9688499663619778E+00,-.9688499663619778E+00, 
    -.6306801197316696E+00,0.6306801197316697E+00, 
    0.8526157293336619E+00,-.8526157293336618E+00, 
    0.7620832819261708E-01,-.7620832819261709E-01, 
    0.5237358202144290E+00,0.4533398211356468E+00, 
    -.4533398211356468E+00,-.5237358202144290E+00, 
    -.1536427274631298E-17 };
  static double ws[] = {
    0.7926638883415150E-01,0.7926638883415155E-01, 
    0.7926638883415171E-01,0.7926638883415174E-01, 
    0.6284721101179108E-01,0.6284721101179105E-01, 
    0.6284721101179126E-01,0.6284721101179125E-01, 
    0.1902480253324007E+00,0.1902480253324007E+00, 
    0.1902480253324001E+00,0.1902480253324001E+00, 
    0.2816282136297291E+00,0.2816282136297291E+00, 
    0.2816282136297291E+00,0.2816282136297291E+00, 
    0.3724677695139016E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule10 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE10 returns the rule of degree 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
     -0.9853833119314600E+00, 
     -0.9262105001258388E+00, 
     -0.9119588710357346E+00, 
     -0.8792348043990323E+00, 
      0.7551269206143556E+00, 
     -0.3188453596839296E+00, 
     -0.6474946981752547E+00, 
      0.9492191314088700E+00, 
     -0.6188661913929927E+00, 
     -0.9215290755789827E+00, 
     -0.1043123255663635E+00, 
      0.9707183739677747E+00, 
      0.9246684242905355E+00, 
     -0.1655785251003832E+00, 
     -0.1844717206212201E+00, 
      0.6666473305982110E+00, 
     -0.6115967830349248E+00, 
      0.3187935759364068E+00, 
      0.3663139167806795E+00, 
      0.2105273891482153E+00, 
      0.7631114939243835E+00, 
      0.6258661935323967E+00 };
  static double ys[] = {
      0.6240243795898477E+00, 
     -0.9577495916000753E+00, 
      0.2065013461988724E+00, 
      0.9225481682574119E+00, 
     -0.9538019223425510E+00, 
      0.9406185571992116E+00, 
      0.6476354842626755E+00, 
      0.1142951736422384E+00, 
     -0.1073322786510873E+00, 
     -0.5089131904296068E+00, 
     -0.9663420836873585E+00, 
     -0.7529656324799600E+00, 
      0.8117151060164874E+00, 
     -0.4732489884927658E+00, 
      0.3507267260891900E+00, 
      0.4711392149070170E+00, 
     -0.8103749226019182E+00, 
     -0.3211002312038632E-01, 
     -0.7605065507139738E+00, 
      0.7788254159831852E+00, 
     -0.3924748753960959E+00, 
      0.9817119264047970E+00 };
  static double ws[] = {
      0.2360826549495683E-01, 
      0.2437222107469133E-01, 
      0.9518348391179431E-01, 
      0.4422100229516362E-01, 
      0.5077443794275623E-01, 
      0.8474079765068482E-01, 
      0.1630376427353708E+00, 
      0.8387536029255430E-01, 
      0.2227307039000182E+00, 
      0.9431200171676321E-01, 
      0.6925434184055250E-01, 
      0.4179721982807350E-01, 
      0.6271061498530471E-01, 
      0.2628440705130901E+00, 
      0.2744383764593098E+00, 
      0.2007889259004832E+00, 
      0.1415509025856471E+00, 
      0.2886706704866435E+00, 
      0.1841475599924024E+00, 
      0.1879992823185209E+00, 
      0.1826311110112870E+00, 
      0.4473813181012295E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule11 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE11 returns the rule of degree 11.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.1885861387186414E+00,0.9535395282015320E+00, 
    -.1885861387186413E+00,-.9535395282015320E+00, 
    -.6980761045495679E+00,0.9826392235408555E+00, 
    0.6980761045495681E+00,-.9826392235408555E+00, 
    -.9394863828167370E+00,0.8257758359029639E+00, 
    0.9394863828167370E+00,-.8257758359029637E+00, 
    -.7120019130753364E+00,0.5253202503645475E+00, 
    0.7120019130753364E+00,-.5253202503645475E+00, 
    -.3156234329152547E+00,0.8125205483048131E+00, 
    0.3156234329152548E+00,-.8125205483048131E+00, 
    -.4248472488486694E+00,0.4165807191202203E-01, 
    0.4248472488486694E+00,-.4165807191202197E-01 };
  static double ys[] = {
    -.9535395282015320E+00,0.1885861387186414E+00, 
    0.9535395282015320E+00,-.1885861387186413E+00, 
    -.9826392235408555E+00,-.6980761045495680E+00, 
    0.9826392235408555E+00,0.6980761045495683E+00, 
    -.8257758359029640E+00,-.9394863828167370E+00, 
    0.8257758359029638E+00,0.9394863828167370E+00, 
    -.5253202503645475E+00,-.7120019130753364E+00, 
    0.5253202503645475E+00,0.7120019130753364E+00, 
    -.8125205483048131E+00,-.3156234329152547E+00, 
    0.8125205483048131E+00,0.3156234329152549E+00, 
    -.4165807191202205E-01,-.4248472488486694E+00, 
    0.4165807191202200E-01,0.4248472488486694E+00 };
  static double ws[] = {
    0.6886285066821880E-01,0.6886285066821880E-01, 
    0.6886285066821880E-01,0.6886285066821880E-01, 
    0.3395580740305121E-01,0.3395580740305121E-01, 
    0.3395580740305121E-01,0.3395580740305121E-01, 
    0.4671948489426224E-01,0.4671948489426224E-01, 
    0.4671948489426224E-01,0.4671948489426224E-01, 
    0.1595417182608939E+00,0.1595417182608939E+00, 
    0.1595417182608939E+00,0.1595417182608939E+00, 
    0.1497202089079448E+00,0.1497202089079448E+00, 
    0.1497202089079448E+00,0.1497202089079448E+00, 
    0.2483067110521767E+00,0.2483067110521767E+00, 
    0.2483067110521767E+00,0.2483067110521767E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule12 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE12 returns the rule of degree 12.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.9711185107918885E+00,0.6489450771045480E+00, 
    -.9547543104262661E+00,-.9065777988000044E+00, 
    0.9288045791287373E+00,-.9425162358139516E+00, 
    -.9438108523148829E+00,-.6477285885089000E+00, 
    0.9399983037047607E+00,-.6866282782659429E+00, 
    0.7379913501268124E+00,-.3293152288819712E+00, 
    0.6556399582616308E+00,0.4257309111871534E+00, 
    0.9692829476897494E+00,-.8505721097622355E+00, 
    0.2079264382173936E+00,-.8025201782903676E+00, 
    -.5197237466355563E+00,-.2734035281447398E-02, 
    -.3699658428845123E+00,-.6558970744607242E+00, 
    0.3202734978144128E+00,0.8244469498554706E+00, 
    0.2278925123542080E-01,0.5651896934359196E-01, 
    0.4889968338954821E+00,-.3555976156822369E+00, 
    0.7575512967066254E+00,-.1870234315112276E+00, 
    -.9967741042631649E+00 };
  static double ys[] = {
    -.8672105154213969E+00,0.9928922644702000E+00, 
    0.2857493181383339E+00,0.9656011790176721E+00, 
    0.8921207951072256E+00,-.9100219543607504E+00, 
    0.7000393501436398E+00,-.9664634507775797E+00, 
    0.4769006675305678E+00,0.4257467094739614E+00, 
    -.9615153562096814E+00,0.9693604253119810E+00, 
    0.7113042249747283E+00,-.7943285461026974E+00, 
    -.1992840398255900E+00,-.7990773209000775E-01, 
    0.8876704665740045E+00,-.6649760891057823E+00, 
    -.3390779542043381E+00,-.5354471390418425E+00, 
    0.1502820099360215E+00,0.8334029046137713E+00, 
    0.3881559105148546E+00,-.5879856922234445E+00, 
    -.2000991752759776E-01,-.9646721637922943E+00, 
    -.2761642039851812E+00,-.8128294162538594E+00, 
    0.1215344546399007E+00,0.6390274229299343E+00, 
    -.4400036004541968E+00 };
  static double ws[] = {
    0.2294319212922989E-01,0.2269718640167010E-01, 
    0.3814000586151853E-01,0.1921567026521910E-01, 
    0.3433859117158319E-01,0.2503589002871782E-01, 
    0.4151906822977771E-01,0.3384747145223248E-01, 
    0.5960510578836526E-01,0.1273691684426847E+00, 
    0.3629732156183973E-01,0.4288352023218015E-01, 
    0.1186836445978463E+00,0.1223937757234154E+00, 
    0.4775972626669994E-01,0.1037607311404478E+00, 
    0.1017344934330748E+00,0.9441812422392200E-01, 
    0.1662942328844954E+00,0.1752158503094866E+00, 
    0.1535404788337684E+00,0.8331450401650711E-01, 
    0.1951461758691787E+00,0.1055576202902579E+00, 
    0.1749572560557213E+00,0.5131431669876880E-01, 
    0.1804321296492865E+00,0.1127530309084298E+00, 
    0.1400307997981144E+00,0.1721261711453589E+00, 
    0.2510187133639127E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule13 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE13 returns the rule of degree 13.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.9572976997863074E+00,0.8595560056416388E+00, 
    0.9572976997863074E+00,-.8595560056416386E+00, 
    -.7788097115544195E+00,0.9834866824398722E+00, 
    0.7788097115544197E+00,-.9834866824398722E+00, 
    -.4758086252182752E+00,0.8500766736997490E+00, 
    0.4758086252182753E+00,-.8500766736997490E+00, 
    0.3907362161294613E+00,0.9413272258729251E+00, 
    -.3907362161294612E+00,-.9413272258729251E+00, 
    -.1381834598624646E+00,0.9589251702875351E+00, 
    0.1381834598624647E+00,-.9589251702875351E+00, 
    0.6478216371870111E+00,0.7558053565720809E+00, 
    -.6478216371870111E+00,-.7558053565720809E+00, 
    0.7074150899644462E-01,0.6962500784917495E+00, 
    -.7074150899644453E-01,-.6962500784917495E+00, 
    -.3427165560404070E+00,0.4093045616940387E+00, 
    0.3427165560404070E+00,-.4093045616940387E+00, 
    -.7375869198366919E-30 };
  static double ys[] = {
    -.8595560056416389E+00,-.9572976997863074E+00, 
    0.8595560056416387E+00,0.9572976997863074E+00, 
    -.9834866824398722E+00,-.7788097115544196E+00, 
    0.9834866824398722E+00,0.7788097115544198E+00, 
    -.8500766736997490E+00,-.4758086252182752E+00, 
    0.8500766736997490E+00,0.4758086252182753E+00, 
    -.9413272258729251E+00,0.3907362161294612E+00, 
    0.9413272258729251E+00,-.3907362161294611E+00, 
    -.9589251702875351E+00,-.1381834598624647E+00, 
    0.9589251702875351E+00,0.1381834598624648E+00, 
    -.7558053565720809E+00,0.6478216371870111E+00, 
    0.7558053565720809E+00,-.6478216371870111E+00, 
    -.6962500784917495E+00,0.7074150899644457E-01, 
    0.6962500784917495E+00,-.7074150899644449E-01, 
    -.4093045616940387E+00,-.3427165560404070E+00, 
    0.4093045616940387E+00,0.3427165560404070E+00, 
    -.6522588594679827E-30 };
  static double ws[] = {
    0.2699339218118215E-01,0.2699339218118215E-01, 
    0.2699339218118215E-01,0.2699339218118215E-01, 
    0.2120743264134157E-01,0.2120743264134157E-01, 
    0.2120743264134157E-01,0.2120743264134157E-01, 
    0.8403587015611026E-01,0.8403587015611026E-01, 
    0.8403587015611026E-01,0.8403587015611026E-01, 
    0.5479564090947502E-01,0.5479564090947502E-01, 
    0.5479564090947502E-01,0.5479564090947502E-01, 
    0.4272687338421139E-01,0.4272687338421139E-01, 
    0.4272687338421139E-01,0.4272687338421139E-01, 
    0.9175668641747110E-01,0.9175668641747110E-01, 
    0.9175668641747110E-01,0.9175668641747110E-01, 
    0.1508552789574409E+00,0.1508552789574409E+00, 
    0.1508552789574409E+00,0.1508552789574409E+00, 
    0.1816350488471704E+00,0.1816350488471704E+00, 
    0.1816350488471704E+00,0.1816350488471704E+00, 
    0.2124022307685795E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule14 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE14 returns the rule of degree 14.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.7971028933442358E-01,0.5471061974611687E+00, 
    0.9884165047157199E-02,-.3207474755493107E+00, 
    -.9553717233098673E-01,0.4082649649717013E+00, 
    0.7282812616257787E+00,-.3472810864715026E+00, 
    -.7915820141518555E+00,0.2785829529075714E+00, 
    0.2099338375706851E+00,-.3770688843422310E+00, 
    0.1546298932128112E+00,0.9706353820201650E+00, 
    -.9710188308157890E+00,-.8089422169597223E-01, 
    -.6907103686706517E+00,-.3766500271820025E+00, 
    -.5730050261123937E+00,0.2628896117741670E+00, 
    0.5452072294160648E+00,-.8660081844945786E+00, 
    0.5293743439059322E+00,0.8326824008480398E+00, 
    -.9838433394477736E+00,0.9724693366798177E+00, 
    0.7239917862034132E+00,-.8848792748044129E+00, 
    0.9060729882068399E+00,0.9965542304973446E+00, 
    -.8318061454455522E+00,-.3342315451147794E+00, 
    0.7013190684905666E+00,0.9589226582728134E+00, 
    -.6067867745349950E+00,0.9099368375229098E+00, 
    -.9529480307145739E+00,-.9784887714287833E+00, 
    -.6288709123199404E+00,-.8984845926717366E+00, 
    0.7834593049247998E+00 };
  static double ys[] = {
    0.7997982556617036E+00,0.8939858342635323E+00, 
    0.4772113760639611E+00,0.7590210924614017E+00, 
    -.1641336121692613E+00,0.6627639082081634E+00, 
    0.9910176287834889E+00,0.1897335228068912E+00, 
    0.1071320357310403E+00,0.1628717488606178E+00, 
    -.4874120187172733E+00,0.9815580199193836E+00, 
    0.9428958078193180E+00,-.9568920474226420E+00, 
    0.2935671028498448E+00,-.8120950161302742E+00, 
    0.8737574543065562E+00,-.5687517308443202E+00, 
    0.5033945159864817E+00,-.9476918771401073E+00, 
    -.1446903403457750E+00,0.6085144943023661E+00, 
    -.7379657281768257E+00,0.7553490523545380E+00, 
    -.7955014974234428E+00,0.5032746878924270E+00, 
    -.9513994070448035E+00,-.9641836389975580E+00, 
    -.7647494365215827E+00,-.4197122679990931E+00, 
    -.5991295673421837E+00,-.9870420901808716E+00, 
    0.3662527369607260E+00,0.9255947941322170E+00, 
    -.2432439926029902E+00,0.2516141162390789E-01, 
    -.2593532979101687E+00,0.8203127557579633E+00, 
    -.8613154366604666E+00,0.9766163862390111E+00, 
    -.4322894664059549E+00 };
  static double ws[] = {
    0.5306119137159240E-01,0.4422789838333461E-01, 
    0.1273937944867034E+00,0.7175684820689888E-01, 
    0.1568138398057906E+00,0.1082893191419578E+00, 
    0.1425381117124736E-01,0.1362424642572790E+00, 
    0.9621728531430029E-01,0.1420579284045241E+00, 
    0.1461917185380448E+00,0.2659569714822968E-01, 
    0.5071663615284024E-01,0.8265110949266337E-02, 
    0.3089045691835443E-01,0.9808292446464849E-01, 
    0.5780479268214644E-01,0.1229597898540734E+00, 
    0.1040216750542093E+00,0.4885605487211036E-01, 
    0.1309068037017440E+00,0.5947734706318328E-01, 
    0.9812697893659589E-01,0.5941729254373463E-01, 
    0.1721117182291158E-01,0.3143615055818434E-01, 
    0.3320596149477782E-01,0.1938593289821144E-01, 
    0.4495089530188040E-01,0.1917919328230456E-01, 
    0.7662013919695503E-01,0.2566629926215708E-01, 
    0.1191871942506123E+00,0.1558050723920070E-01, 
    0.1194301874701882E+00,0.7623831743082932E-01, 
    0.4693631023011145E-01,0.1603574648435874E-01, 
    0.7006443647504916E-01,0.1296824822458430E-01, 
    0.9170277370106396E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule15 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE15 returns the rule of degree 15.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.7749527857778351E+00,0.9885448991378063E+00, 
    -.7749527857778349E+00,-.9885448991378063E+00, 
    -.9070374303651182E+00,0.9571446613308432E+00, 
    0.9070374303651184E+00,-.9571446613308430E+00, 
    -.4303978306869286E+00,0.9769578054468787E+00, 
    0.4303978306869287E+00,-.9769578054468787E+00, 
    -.9756646723906326E+00,0.1107064048513496E+00, 
    0.9756646723906326E+00,-.1107064048513495E+00, 
    -.7388921437312957E+00,0.7868610204187212E+00, 
    0.7388921437312957E+00,-.7868610204187212E+00, 
    0.1995220876718269E+00,0.6659287668239546E+00, 
    -.1995220876718268E+00,-.6659287668239546E+00, 
    -.1934983412061240E+00,0.8412271039808018E+00, 
    0.1934983412061241E+00,-.8412271039808018E+00, 
    0.4882189227791580E+00,0.8922368778153702E+00, 
    -.4882189227791579E+00,-.8922368778153702E+00, 
    -.5772265461040059E+00,0.9526539504944950E+00, 
    0.5772265461040061E+00,-.9526539504944950E+00, 
    -.4474426063114782E+00,0.5675455860909890E+00, 
    0.4474426063114783E+00,-.5675455860909890E+00, 
    -.7044956995149931E-01,0.3256679896817100E+00, 
    0.7044956995149934E-01,-.3256679896817100E+00 };
  static double ys[] = {
    -.9885448991378063E+00,0.7749527857778350E+00, 
    0.9885448991378063E+00,-.7749527857778348E+00, 
    -.9571446613308433E+00,-.9070374303651183E+00, 
    0.9571446613308431E+00,0.9070374303651185E+00, 
    -.9769578054468787E+00,-.4303978306869286E+00, 
    0.9769578054468787E+00,0.4303978306869287E+00, 
    -.1107064048513496E+00,-.9756646723906326E+00, 
    0.1107064048513495E+00,0.9756646723906326E+00, 
    -.7868610204187212E+00,-.7388921437312957E+00, 
    0.7868610204187212E+00,0.7388921437312957E+00, 
    -.6659287668239546E+00,0.1995220876718268E+00, 
    0.6659287668239546E+00,-.1995220876718268E+00, 
    -.8412271039808018E+00,-.1934983412061240E+00, 
    0.8412271039808018E+00,0.1934983412061241E+00, 
    -.8922368778153702E+00,0.4882189227791580E+00, 
    0.8922368778153702E+00,-.4882189227791578E+00, 
    -.9526539504944950E+00,-.5772265461040060E+00, 
    0.9526539504944950E+00,0.5772265461040063E+00, 
    -.5675455860909890E+00,-.4474426063114783E+00, 
    0.5675455860909890E+00,0.4474426063114784E+00, 
    -.3256679896817100E+00,-.7044956995149933E-01, 
    0.3256679896817100E+00,0.7044956995149936E-01 };
  static double ws[] = {
    0.1443015463807196E-01,0.1443015463807196E-01, 
    0.1443015463807196E-01,0.1443015463807196E-01, 
    0.1816242330920956E-01,0.1816242330920956E-01, 
    0.1816242330920956E-01,0.1816242330920956E-01, 
    0.1290815898308381E-01,0.1290815898308381E-01, 
    0.1290815898308381E-01,0.1290815898308381E-01, 
    0.3010764365372140E-01,0.3010764365372140E-01, 
    0.3010764365372140E-01,0.3010764365372140E-01, 
    0.6540469907131932E-01,0.6540469907131932E-01, 
    0.6540469907131932E-01,0.6540469907131932E-01, 
    0.1197895531736646E+00,0.1197895531736646E+00, 
    0.1197895531736646E+00,0.1197895531736646E+00, 
    0.8473841548096289E-01,0.8473841548096289E-01, 
    0.8473841548096289E-01,0.8473841548096289E-01, 
    0.6453833756714425E-01,0.6453833756714425E-01, 
    0.6453833756714425E-01,0.6453833756714425E-01, 
    0.2403055376316494E-01,0.2403055376316494E-01, 
    0.2403055376316494E-01,0.2403055376316494E-01, 
    0.1196130510491228E+00,0.1196130510491228E+00, 
    0.1196130510491228E+00,0.1196130510491228E+00, 
    0.1533837904970821E+00,0.1533837904970821E+00, 
    0.1533837904970821E+00,0.1533837904970821E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule16 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE16 returns the rule of degree 16.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.4455077315117606E+00,-.2848639978436560E+00, 
    0.3281612731066613E+00,0.8918323090897292E+00, 
    0.5080739455207711E+00,-.1604014838336734E+00, 
    -.5084786013975685E+00,0.6830501661618820E+00, 
    0.3958609078545774E-01,-.9505884228486957E+00, 
    -.9618117119546932E+00,-.9824876605856204E+00, 
    0.7294909884375445E+00,0.1663512176814138E+00, 
    -.8401346406026955E+00,0.7196828159773350E+00, 
    0.2530272226250717E+00,0.6403128558574714E-01, 
    0.9867649245062882E+00,-.2009570267598735E+00, 
    -.4385685045405642E+00,0.9942720499027967E+00, 
    -.8306577686854565E+00,0.1540093368841131E+00, 
    0.5127230397362322E+00,-.9896697019211779E+00, 
    0.9721677355950014E+00,-.7432189863550824E+00, 
    0.7845206176595362E+00,0.9792382924978336E+00, 
    -.7071472873150695E+00,-.1416891896032315E+00, 
    0.8901118844470474E+00,0.8991291622644930E+00, 
    -.8948787586840297E+00,-.1685236423690618E+00, 
    0.8814614259143865E+00,-.9727794481100781E+00, 
    -.2335640652756550E+00,0.9317484014546712E+00, 
    0.4462548323794501E+00,-.4956609448877936E+00, 
    -.7464462325957708E+00,-.6209288757442605E+00, 
    -.9729904625306640E+00,0.6851689211463186E+00, 
    -.6950237982778372E+00,0.1406973647399049E+00, 
    -.8903519744155600E+00,-.4494459013630722E+00, 
    0.6393855383975146E+00,0.4185604687784284E+00 };
  static double ys[] = {
    -.4133544718603057E+00,-.7739045957079853E+00, 
    -.9418618376442869E+00,-.5606355473579230E+00, 
    0.9698717423735355E+00,-.4892197731244818E+00, 
    -.9117050044897801E+00,-.5146278686624721E+00, 
    -.8810600362311165E+00,-.9545711191201147E+00, 
    -.5466258413042109E+00,-.8038321990753307E+00, 
    -.1738438525238659E+00,-.6521123405187266E+00, 
    -.3199762775663758E+00,0.8583352444264194E+00, 
    0.9559713160961977E+00,0.8638442022353128E+00, 
    -.7782646945578277E+00,0.6577321912732575E+00, 
    -.2429524974731521E+00,0.3876054361677954E+00, 
    -.7951217128043849E+00,-.2163105603296342E+00, 
    -.7522598930908062E+00,0.5914422527084512E+00, 
    -.2390927720573675E+00,-.9749994483165604E+00, 
    -.8444622014821822E+00,0.8810775159513397E+00, 
    0.6173187727164633E+00,0.7292324652124864E-01, 
    0.1389731666217061E+00,0.6539121663076461E+00, 
    0.8152409648202001E+00,-.9939582277051608E+00, 
    0.9818676808679303E+00,0.9582010007957961E+00, 
    0.9824162852107247E+00,-.9611203598796459E+00, 
    0.1014796681115945E+00,0.8630083224142983E+00, 
    0.9693553936129727E+00,-.5969027498090264E+00, 
    -.3789880325419188E-01,0.4272704958652744E+00, 
    0.4313302857584844E-01,0.4101941952194764E+00, 
    0.3277950128817075E+00,0.3540266909627898E+00, 
    -.9786311482371555E+00,0.6960968466161754E+00 };
  static double ws[] = {
    0.7760560264177564E-01,0.6557384620049388E-01, 
    0.3492505367961147E-01,0.4902394774171961E-01, 
    0.2071780611220039E-01,0.9663554422001587E-01, 
    0.4275428371393922E-01,0.4251351999693562E-01, 
    0.4614134112945410E-01,0.1008804262208153E-01, 
    0.2597909269494542E-01,0.7243920128199269E-02, 
    0.8689727864752500E-01,0.7770693383996630E-01, 
    0.6824785727785596E-01,0.4942860806128915E-01, 
    0.1939352812224098E-01,0.6287664325857020E-01, 
    0.1228661441972498E-01,0.1029914685037788E+00, 
    0.1129049857606719E+00,0.1469344234054041E-01, 
    0.4431063501362172E-01,0.1183113207778475E+00, 
    0.6355581532812171E-01,0.1527602237522992E-01, 
    0.2841229216361006E-01,0.1472432061291935E-01, 
    0.4129113754659276E-01,0.1149253014935389E-01, 
    0.7883123749996855E-01,0.1339728049062325E+00, 
    0.6215863520566767E-01,0.4742159015437373E-01, 
    0.3530277887828205E-01,0.1663690106625590E-01, 
    0.1074607484017106E-01,0.6944512701761030E-02, 
    0.2298755190493018E-01,0.1230850655669064E-01, 
    0.1267827808948161E+00,0.6295607382580996E-01, 
    0.2066068139152711E-01,0.8653018393026629E-01, 
    0.2938486912464724E-01,0.9600456501593493E-01, 
    0.9622610521466571E-01,0.1308302132754554E+00, 
    0.6115363632570159E-01,0.1144365242523608E+00, 
    0.1628247685682751E-01,0.9586498584301183E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule17 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE17 returns the rule of degree 17.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.7710386602263628E+00,0.7710386602263630E+00, 
    0.9803457456469387E+00,-.9803457456469384E+00, 
    -.2292639639675523E+00,0.2292639639675524E+00, 
    0.4847176019505991E-03,-.4847176019504780E-03, 
    -.6189416389750175E+00,0.6189416389750177E+00, 
    0.9587315519802511E+00,-.9587315519802511E+00, 
    0.8409306922533750E+00,-.8409306922533748E+00, 
    -.4308042054877432E+00,0.4308042054877433E+00, 
    0.4761431266211590E+00,-.4761431266211589E+00, 
    0.8651144531733139E+00,-.8651144531733139E+00, 
    0.9846617345267017E+00,-.9846617345267017E+00, 
    -.7981639404863030E+00,0.7981639404863030E+00, 
    0.6877591943414725E+00,-.6877591943414725E+00, 
    -.3038305486106544E+00,0.3038305486106544E+00, 
    0.9852576255116258E+00,-.9852576255116258E+00, 
    0.9853756930046446E+00,-.9853756930046446E+00, 
    0.7024672194580522E+00,-.7024672194580522E+00, 
    0.4589513024499020E+00,-.4589513024499019E+00, 
    -.5838938372432102E+00,0.5838938372432102E+00, 
    0.4855363777625971E+00,-.4855363777625971E+00, 
    0.1909552287968119E+00,-.1909552287968118E+00, 
    0.1970910744873101E+00,-.1970910744873101E+00, 
    0.9070140000742543E+00,-.9070140000742543E+00, 
    -.9370706813548184E+00,0.9370706813548186E+00, 
    -.1024098809482286E+00,0.1024098809482287E+00, 
    0.9018657853563646E+00,-.9018657853563646E+00, 
    0.7422255782894629E+00,-.7422255782894629E+00, 
    -.1975779250586182E-19 };
  static double ys[] = {
    -.9187170657318696E+00,0.9187170657318696E+00, 
    -.9679135253250817E+00,0.9679135253250819E+00, 
    -.9437800394025085E+00,0.9437800394025085E+00, 
    -.9886578344699537E+00,0.9886578344699537E+00, 
    -.9803491213417113E+00,0.9803491213417113E+00, 
    -.8226737868824753E+00,0.8226737868824755E+00, 
    -.9649601466712245E+00,0.9649601466712245E+00, 
    -.8370492275539414E+00,0.8370492275539414E+00, 
    -.9716943047473653E+00,0.9716943047473653E+00, 
    -.6326447362896030E+00,0.6326447362896030E+00, 
    0.2029425559112923E+00,-.2029425559112922E+00, 
    -.7906135688735062E+00,0.7906135688735062E+00, 
    -.8442560578129694E+00,0.8442560578129694E+00, 
    -.3117615836793495E+00,0.3117615836793495E+00, 
    0.7701659795648228E+00,-.7701659795648226E+00, 
    -.4379432170880169E+00,0.4379432170880170E+00, 
    -.3820619012323893E+00,0.3820619012323894E+00, 
    -.6514286057161101E+00,0.6514286057161101E+00, 
    -.5711068454496305E+00,0.5711068454496305E+00, 
    -.8072896746317025E-01,0.8072896746317031E-01, 
    -.8630149364726712E+00,0.8630149364726712E+00, 
    -.3872678175415290E+00,0.3872678175415290E+00, 
    0.5103334842355030E+00,-.5103334842355027E+00, 
    -.9584329986119476E+00,0.9584329986119474E+00, 
    -.6619201369182062E+00,0.6619201369182062E+00, 
    -.1238115372273944E+00,0.1238115372273945E+00, 
    0.2071876599633523E+00,-.2071876599633522E+00, 
    0.5346688849930886E-20 };
  static double ws[] = {
    0.1261638293838951E-01,0.1261638293838951E-01, 
    0.3408339905429266E-02,0.3408339905429266E-02, 
    0.2796862081921473E-01,0.2796862081921473E-01, 
    0.1252812914329644E-01,0.1252812914329644E-01, 
    0.1635296523785200E-01,0.1635296523785200E-01, 
    0.1720881227455075E-01,0.1720881227455075E-01, 
    0.1523407270818440E-01,0.1523407270818440E-01, 
    0.5600796522816800E-01,0.5600796522816800E-01, 
    0.2382823797668716E-01,0.2382823797668716E-01, 
    0.4513279974663867E-01,0.4513279974663867E-01, 
    0.1931215256841166E-01,0.1931215256841166E-01, 
    0.4158804216001467E-01,0.4158804216001467E-01, 
    0.4685849665862760E-01,0.4685849665862760E-01, 
    0.1200522449400290E+00,0.1200522449400290E+00, 
    0.1238565802221201E-01,0.1238565802221201E-01, 
    0.1760077392303538E-01,0.1760077392303538E-01, 
    0.8264937698824523E-01,0.8264937698824523E-01, 
    0.8629133710270168E-01,0.8629133710270168E-01, 
    0.8660536182880048E-01,0.8660536182880048E-01, 
    0.1134857467272575E+00,0.1134857467272575E+00, 
    0.6518861145910534E-01,0.6518861145910534E-01, 
    0.1184802238173896E+00,0.1184802238173896E+00, 
    0.4767526390300979E-01,0.4767526390300979E-01, 
    0.1203076112968188E-01,0.1203076112968188E-01, 
    0.1010849820160845E+00,0.1010849820160845E+00, 
    0.5753445241741756E-01,0.5753445241741756E-01, 
    0.8946701652955226E-01,0.8946701652955226E-01, 
    0.1312734684062163E+00 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule18 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE18 returns the rule of degree 18.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.9824292972819758E+00,-.9803800150830472E+00, 
    0.9786873970954008E+00,0.9783064717068459E+00, 
    0.9721036514650782E+00,0.5934299192489508E+00, 
    0.5976341012197494E+00,-.3692961886011435E+00, 
    0.7049479055795181E+00,-.7373468008004338E+00, 
    0.9861989064660380E+00,-.8465158209874328E-01, 
    0.8850139974482023E+00,-.9313926150284875E+00, 
    -.7089981957319805E+00,0.8732614092173598E+00, 
    -.8419399215622921E+00,0.3224279462028189E+00, 
    0.2244640124522700E+00,-.8347180004932071E+00, 
    -.8279536304784142E+00,-.6568125381001307E+00, 
    -.8460042964230963E+00,-.9242174294702851E+00, 
    -.5414368844899894E+00,0.9265304841706828E+00, 
    -.8981156173043389E+00,-.4115732120271069E+00, 
    0.6047453622518654E+00,-.5983900835098156E+00, 
    0.4320397790648277E+00,0.4499667234152572E+00, 
    0.9513764564603604E+00,-.6795979700947666E+00, 
    -.6564652118524049E+00,0.8378047852517292E+00, 
    -.1113556772390856E+00,-.3608746962710803E+00, 
    0.9862753904173321E+00,-.2331918978886594E+00, 
    -.9870492227974196E+00,-.9894160854750882E+00, 
    -.9836502238795436E+00,0.8002432442611336E+00, 
    0.6492404239850669E+00,0.8701160725699196E+00, 
    0.7389321385106327E+00,-.4715855030569622E+00, 
    0.7436165358202136E+00,0.3534834684769252E+00, 
    -.3650849452075164E+00,-.2283506208713122E+00, 
    0.1740274102467612E+00,0.8995839783500573E+00, 
    0.3405925655026003E+00,0.4295328792252409E-01, 
    0.1448653016063912E+00,-.1312827821241097E+00, 
    0.5004638295686703E+00,0.6006609259220026E-01, 
    0.1552257232885370E+00,-.9395864632354990E+00, 
    -.4791507758959901E+00,-.8137902397325638E+00 };
  static double ys[] = {
    -.9492304540136947E+00,0.9469609026131530E+00, 
    0.7622186219256711E+00,-.9196983748762729E+00, 
    0.9764183470032382E+00,-.9616466720526158E+00, 
    0.8343589711117383E+00,0.9800036938552533E+00, 
    0.9870747330066914E+00,-.8333460482416958E+00, 
    -.3278737740897394E+00,0.8908466843089227E+00, 
    -.9855803626832146E+00,-.8212845135474046E+00, 
    0.6822143912173052E+00,0.9031148761319726E+00, 
    0.9885786707367459E+00,0.5365154393395204E-01, 
    0.7350437600270140E+00,0.5084609108470370E+00, 
    -.6132315128115761E+00,0.2290780703350274E+00, 
    -.9786418155875946E+00,-.3102577845761006E+00, 
    -.6275098391746888E+00,0.1407483577400233E-01, 
    0.8295680022967102E+00,0.7753196502785407E+00, 
    -.6526316599009200E+00,-.9333674622588783E+00, 
    0.9243432783386737E+00,0.7005646989935736E-01, 
    -.6925969506257750E+00,-.3641945758998874E+00, 
    0.9193665018045939E+00,-.4639915525654867E+00, 
    -.9041936920710849E+00,-.7808865818750257E+00, 
    0.2974467982752778E+00,-.4343746051515523E+00, 
    -.5898322990766285E+00,-.1604383366455782E-01, 
    0.6393477839238390E+00,-.8481629731994191E+00, 
    -.2054507717193144E+00,-.1017355764563722E+00, 
    0.7159568279790096E+00,-.1245817088489083E+00, 
    0.2542532076804419E+00,-.8415199898068800E+00, 
    -.9897219383699781E+00,0.1556384545656453E+00, 
    -.9791074324031704E+00,0.5295267411341480E+00, 
    -.4165391069331053E+00,-.1509802365687459E+00, 
    0.3347376662738560E+00,0.5759504751299614E+00, 
    0.5086985983019631E+00,-.6661951769904182E+00, 
    0.9847264496502689E+00,0.3003342534564235E+00, 
    0.4492326697835803E+00,-.6261874651130582E-02 };
  static double ws[] = {
    0.4880182878194577E-02,0.5801719791846293E-02, 
    0.1363202620128554E-01,0.7756485561516194E-02, 
    0.4249395346845679E-02,0.2449228653823736E-01, 
    0.2751148878355652E-01,0.1865262619502059E-01, 
    0.1160827955745764E-01,0.2909451783554787E-01, 
    0.1614009453230349E-01,0.5169172255655628E-01, 
    0.7388591103278081E-02,0.2021700469490708E-01, 
    0.4636019735255641E-01,0.2315518049910366E-01, 
    0.7815468109974862E-02,0.3769198715796499E-01, 
    0.7881010252197586E-01,0.3988856078693273E-01, 
    0.4079805534552380E-01,0.6645673137795649E-01, 
    0.1072546708888146E-01,0.3868425989902273E-01, 
    0.5522501374016728E-01,0.3253997663365805E-01, 
    0.2752684329481517E-01,0.6117844466269161E-01, 
    0.7168779734796255E-01,0.2620423047578233E-01, 
    0.3392741843404458E-01,0.7458228761237336E-01, 
    0.2478342723311915E-01,0.6447492653722675E-01, 
    0.3284650513806075E-01,0.5368779603413000E-01, 
    0.4545538155398856E-01,0.5397016565409892E-01, 
    0.1519721989624682E-01,0.1031663549936232E+00, 
    0.1231416020827792E-01,0.1349682430871745E-01, 
    0.1395336957267540E-01,0.3753801076980155E-01, 
    0.8258232039188609E-01,0.2489123627268648E-01, 
    0.4546454333835988E-01,0.9452263223987649E-01, 
    0.7842163445187621E-01,0.5876926206730445E-01, 
    0.1249000951778932E-01,0.1088923354481117E+00, 
    0.1951366247998040E-01,0.4179385297531394E-01, 
    0.1033273611701436E+00,0.1162365062083885E+00, 
    0.1023818705273872E+00,0.8804476252208088E-01, 
    0.8910543971624120E-01,0.8874335328605459E-01, 
    0.1602722062379702E-01,0.3378674145303614E-01, 
    0.7595232871119993E-01,0.6022146552676879E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule19 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE19 returns the rule of degree 19.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.9734386316165470E+00,0.9744990929036832E+00, 
    0.9734386316165472E+00,-.9744990929036830E+00, 
    -.3841574585766744E+00,0.9670641778942685E+00, 
    0.3841574585766745E+00,-.9670641778942685E+00, 
    0.2986734938364671E+00,0.9905525689050123E+00, 
    -.2986734938364670E+00,-.9905525689050123E+00, 
    -.7396581737067777E+00,0.9869464369033261E+00, 
    0.7396581737067779E+00,-.9869464369033261E+00, 
    -.1425244970455050E+00,0.9733021904515969E+00, 
    0.1425244970455051E+00,-.9733021904515969E+00, 
    0.7650240374639232E+00,0.9804863471920530E+00, 
    -.7650240374639230E+00,-.9804863471920530E+00, 
    -.7599006633708002E+00,0.7279453517455540E+00, 
    0.7599006633708002E+00,-.7279453517455540E+00, 
    -.1192987760526789E+00,-.2637912058730560E-02, 
    0.1192987760526789E+00,0.2637912058730575E-02, 
    -.8850504442537889E+00,0.9022234232868145E+00, 
    0.8850504442537889E+00,-.9022234232868145E+00, 
    0.5304174652462883E+00,0.9125489607085608E+00, 
    -.5304174652462881E+00,-.9125489607085608E+00, 
    -.2858528945041368E+00,0.2941600854694212E+00, 
    0.2858528945041368E+00,-.2941600854694212E+00, 
    -.5671850101113227E+00,0.8836081660895880E+00, 
    0.5671850101113227E+00,-.8836081660895880E+00, 
    0.3174295148500719E+00,0.7293427112089215E+00, 
    -.3174295148500718E+00,-.7293427112089215E+00, 
    -.2492430513869149E+00,0.7672563284436533E+00, 
    0.2492430513869150E+00,-.7672563284436533E+00, 
    -.5087793568494521E+00,0.5623738439118215E+00, 
    0.5087793568494521E+00,-.5623738439118215E+00, 
    0.7335719396414396E-01,0.8930925855397183E+00, 
    -.7335719396414385E-01,-.8930925855397183E+00, 
    0.8350775723842838E-02,0.5392457387102469E+00, 
    -.8350775723842772E-02,-.5392457387102469E+00 };
  static double ys[] = {
    -.9744990929036833E+00,-.9734386316165471E+00, 
    0.9744990929036831E+00,0.9734386316165473E+00, 
    -.9670641778942685E+00,-.3841574585766744E+00, 
    0.9670641778942685E+00,0.3841574585766745E+00, 
    -.9905525689050123E+00,0.2986734938364670E+00, 
    0.9905525689050123E+00,-.2986734938364669E+00, 
    -.9869464369033261E+00,-.7396581737067778E+00, 
    0.9869464369033261E+00,0.7396581737067780E+00, 
    -.9733021904515969E+00,-.1425244970455050E+00, 
    0.9733021904515969E+00,0.1425244970455051E+00, 
    -.9804863471920530E+00,0.7650240374639231E+00, 
    0.9804863471920530E+00,-.7650240374639229E+00, 
    -.7279453517455540E+00,-.7599006633708002E+00, 
    0.7279453517455540E+00,0.7599006633708002E+00, 
    0.2637912058730553E-02,-.1192987760526789E+00, 
    -.2637912058730568E-02,0.1192987760526789E+00, 
    -.9022234232868145E+00,-.8850504442537889E+00, 
    0.9022234232868145E+00,0.8850504442537889E+00, 
    -.9125489607085608E+00,0.5304174652462882E+00, 
    0.9125489607085608E+00,-.5304174652462880E+00, 
    -.2941600854694212E+00,-.2858528945041368E+00, 
    0.2941600854694212E+00,0.2858528945041368E+00, 
    -.8836081660895880E+00,-.5671850101113227E+00, 
    0.8836081660895880E+00,0.5671850101113227E+00, 
    -.7293427112089215E+00,0.3174295148500719E+00, 
    0.7293427112089215E+00,-.3174295148500718E+00, 
    -.7672563284436533E+00,-.2492430513869149E+00, 
    0.7672563284436533E+00,0.2492430513869150E+00, 
    -.5623738439118215E+00,-.5087793568494521E+00, 
    0.5623738439118215E+00,0.5087793568494521E+00, 
    -.8930925855397183E+00,0.7335719396414390E-01, 
    0.8930925855397183E+00,-.7335719396414379E-01, 
    -.5392457387102469E+00,0.8350775723842805E-02, 
    0.5392457387102469E+00,-.8350775723842739E-02 };
  static double ws[] = {
    0.4076118519980060E-02,0.4076118519980060E-02, 
    0.4076118519980060E-02,0.4076118519980060E-02, 
    0.1627326938099484E-01,0.1627326938099484E-01, 
    0.1627326938099484E-01,0.1627326938099484E-01, 
    0.1254157952509427E-01,0.1254157952509427E-01, 
    0.1254157952509427E-01,0.1254157952509427E-01, 
    0.1028929333936017E-01,0.1028929333936017E-01, 
    0.1028929333936017E-01,0.1028929333936017E-01, 
    0.1475928282295525E-01,0.1475928282295525E-01, 
    0.1475928282295525E-01,0.1475928282295525E-01, 
    0.1207323692393111E-01,0.1207323692393111E-01, 
    0.1207323692393111E-01,0.1207323692393111E-01, 
    0.4619184040692218E-01,0.4619184040692218E-01, 
    0.4619184040692218E-01,0.4619184040692218E-01, 
    0.3696173437828049E-01,0.3696173437828049E-01, 
    0.3696173437828049E-01,0.3696173437828049E-01, 
    0.2018069481193246E-01,0.2018069481193246E-01, 
    0.2018069481193246E-01,0.2018069481193246E-01, 
    0.3738944032940469E-01,0.3738944032940469E-01, 
    0.3738944032940469E-01,0.3738944032940469E-01, 
    0.9821701539315209E-01,0.9821701539315209E-01, 
    0.9821701539315209E-01,0.9821701539315209E-01, 
    0.3844110871724747E-01,0.3844110871724747E-01, 
    0.3844110871724747E-01,0.3844110871724747E-01, 
    0.7127049386881731E-01,0.7127049386881731E-01, 
    0.7127049386881731E-01,0.7127049386881731E-01, 
    0.6966749913838975E-01,0.6966749913838975E-01, 
    0.6966749913838975E-01,0.6966749913838975E-01, 
    0.7715964130310782E-01,0.7715964130310782E-01, 
    0.7715964130310782E-01,0.7715964130310782E-01, 
    0.4598470092336809E-01,0.4598470092336809E-01, 
    0.4598470092336809E-01,0.4598470092336809E-01, 
    0.9562983140360957E-01,0.9562983140360957E-01, 
    0.9562983140360957E-01,0.9562983140360957E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule20 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE20 returns the rule of degree 20.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.7903555901298286E+00,0.2252273789362463E+00, 
    -.1371556203352229E+00,0.9516912467651903E+00, 
    0.8981438613797404E+00,0.3511327336932318E+00, 
    0.9759068343062448E-01,-.1920781702210826E+00, 
    0.6695436306885109E+00,0.5830289505702392E+00, 
    -.4278190227205511E+00,0.5727845316573248E-01, 
    0.7351507200926066E+00,0.7898934964234535E+00, 
    -.3594399777482915E+00,0.3376839183112454E+00, 
    0.9886396714914109E+00,-.3907721596661082E+00, 
    -.6046899139754487E+00,0.9173007153030285E+00, 
    -.4141396903659944E-01,-.6385982895810176E+00, 
    0.6013396134438070E+00,0.9161661661998538E+00, 
    -.3278639859221721E+00,0.8802460227981840E+00, 
    -.2034581078060147E+00,0.6969166399839279E+00, 
    0.5351593293803870E+00,-.2787747819837111E+00, 
    0.9910314479371701E+00,-.8252889264259545E+00, 
    -.9599491206891601E+00,0.4944916835714113E+00, 
    -.7703896716180660E+00,0.9129604837711394E+00, 
    0.8500373866499284E+00,-.3755453354695931E-01, 
    -.9734574105217517E+00,0.2897016770313429E+00, 
    0.7688438290660360E+00,-.8388618516990838E+00, 
    -.5341286913684105E+00,-.9174806831430735E+00, 
    -.9474346838846874E+00,-.6131936263177837E-01, 
    0.7605933106317154E+00,-.9848886462575595E+00, 
    -.7572987820459326E+00,0.5080745682281066E-01, 
    -.4636550029355001E+00,0.9705319244227489E+00, 
    -.8914079947263784E+00,-.9873056954830064E+00, 
    -.9758313372889643E+00,-.8819170869186215E+00, 
    -.9709635895278415E+00,0.9822722262702339E+00, 
    -.8873503942891041E+00,-.6569247082044029E+00, 
    -.9870905373622706E+00,-.5734674551654146E+00, 
    0.1507888619083396E-01,0.5506050560553492E+00, 
    -.7366839885772312E+00,-.7117026556332434E+00, 
    0.2333965922395885E+00,-.2878681152474931E+00, 
    0.4361746990352844E+00,0.9814633895599361E+00, 
    -.5319601634991685E+00,0.9850987522022563E+00, 
    0.1957972824970721E+00,0.2084217400195901E+00, 
    0.4755027720742199E+00,0.9147494843686150E+00, 
    -.8645689328165518E+00,0.6772212070927294E+00 };
  static double ys[] = {
    -.8898248168845113E-01,-.7280953904485143E-01, 
    -.8653524201560835E+00,-.2902346585197535E+00, 
    0.1218660048097990E+00,-.2471221282702101E+00, 
    -.7467915296200941E+00,-.5589830789081398E+00, 
    -.2090651941599068E+00,-.4355704406099477E+00, 
    -.7264771605644266E+00,-.1363590388950549E+00, 
    -.9878512802669168E+00,0.9209737732078764E+00, 
    -.9535217391501140E+00,-.5953361238026516E+00, 
    -.7515863705429210E-01,0.9865491452054407E+00, 
    -.9940939029178233E+00,-.9840832904192385E+00, 
    0.8284898564271251E-01,-.8660720584668250E+00, 
    0.9809450037856990E+00,0.8048533499129369E+00, 
    0.2790289740151233E+00,-.4082743843038831E+00, 
    0.9740226654958868E+00,0.3035697582018394E+00, 
    -.9862725985120795E+00,0.7185672952170555E+00, 
    -.9193256145694509E+00,-.9544384245126734E+00, 
    -.9898146349839960E+00,0.7338947136259082E-01, 
    -.2075330463780984E+00,0.9907318998968393E+00, 
    0.5387986220469478E+00,0.8834131451429665E+00, 
    -.8934522638909749E+00,-.9207656972246965E+00, 
    -.6367076114754525E+00,0.9978913196317913E+00, 
    0.5300943949356414E+00,-.4262998259673880E+00, 
    0.1841144567160324E+00,0.5090691042712456E+00, 
    -.9214501612581178E+00,-.6400027811642439E+00, 
    0.3209231678914174E+00,-.3846636405176737E+00, 
    0.8697161095697501E+00,0.3470098428369617E+00, 
    0.5711239311935090E+00,0.4067912914297930E+00, 
    0.7624786524458222E+00,-.7770143745259210E+00, 
    0.9584356802550765E+00,0.6551973976459680E+00, 
    0.2854578772720695E-01,0.9551619780680761E+00, 
    -.1675185677029804E+00,0.4240021284329500E-01, 
    -.9864545039606490E+00,-.7989550278861819E+00, 
    -.6071839996966222E+00,0.7388672054082017E+00, 
    0.3140775120692750E+00,-.1793826353802805E+00, 
    0.8772336982049433E+00,-.6033076898494673E+00, 
    -.4039805189289538E+00,0.9282498625490752E+00, 
    0.7185098783764428E+00,0.9758090760986288E+00, 
    0.5388255455349226E+00,-.8020212431755418E+00, 
    0.8814335976449187E+00,0.7393267487868044E+00 };
   static double ws[] = {
    0.4220810999582407E-01,0.4210468769855168E-01, 
    0.4111495137274749E-01,0.1602722818539338E-01, 
    0.3619967055440460E-01,0.5719761023753370E-01, 
    0.5553932136722374E-01,0.7027947643456225E-01, 
    0.3877224057962559E-01,0.5615678549912440E-01, 
    0.5631946977958934E-01,0.1879677326120675E-01, 
    0.1313057162810730E-02,0.2109221815204442E-01, 
    0.2451480331747192E-01,0.6556422061834499E-01, 
    0.1119235292343633E-01,0.7115563056913567E-02, 
    0.6714100241623257E-02,0.5885829448742507E-02, 
    0.8309083447844690E-01,0.3590868298232087E-01, 
    0.1246127664791048E-01,0.2018878910488524E-01, 
    0.8451815870413146E-01,0.3248692905689807E-01, 
    0.1341510929402505E-01,0.6472935261954060E-01, 
    0.1145237221800998E-01,0.5872712613961156E-01, 
    0.4508959666015173E-02,0.1493403238757605E-01, 
    0.2665017934074124E-02,0.7704806248110788E-01, 
    0.5726954359349042E-01,0.4528203295472384E-02, 
    0.4279817302411395E-01,0.4285481643886387E-01, 
    0.8072115177434423E-02,0.3701323444432086E-01, 
    0.4613328158120365E-01,0.3966194341947075E-02, 
    0.6788885180206457E-01,0.3480913280397030E-01, 
    0.1958947629336374E-01,0.7993874305725486E-01, 
    0.2433364291226502E-01,0.1062990044699051E-01, 
    0.6170560315322090E-01,0.7521981238448959E-01, 
    0.3621544028408890E-01,0.1839043728583574E-01, 
    0.3635800635820290E-01,0.1087561990443862E-01, 
    0.1235528041521356E-01,0.2670676504697193E-01, 
    0.5756245567725824E-02,0.1016410330068744E-01, 
    0.3155355506307831E-01,0.1843695052008212E-01, 
    0.1362846461516491E-01,0.8050841953401032E-01, 
    0.1413883594050338E-01,0.5019314017226575E-01, 
    0.5096315858878126E-01,0.4491888769330514E-01, 
    0.8622086129893808E-01,0.9032569802582428E-01, 
    0.4214427133348429E-01,0.1315668261444553E-01, 
    0.7308371765738744E-01,0.5122200740309473E-02, 
    0.6652468674436235E-01,0.1888655008314831E-01, 
    0.7163064878569746E-01,0.2352000357823148E-01, 
    0.2167515013056250E-01,0.4797944511124803E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule21 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE21 returns the rule of degree 21.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.9754554632015029E+00,-.9754554632015029E+00, 
    0.9807628529972264E+00,-.9807628529972262E+00, 
    0.9864691810349792E+00,-.9864691810349792E+00, 
    0.9817752275809579E+00,-.9817752275809579E+00, 
    0.8942647503363536E+00,-.8942647503363533E+00, 
    -.7926169940479519E+00,0.7926169940479522E+00, 
    0.9673238098903247E+00,-.9673238098903247E+00, 
    -.6233424472175393E+00,0.6233424472175395E+00, 
    0.8075494298208362E+00,-.8075494298208362E+00, 
    -.9112339448388380E+00,0.9112339448388382E+00, 
    0.5294700780841287E+00,-.5294700780841285E+00, 
    0.9126054184220483E+00,-.9126054184220483E+00, 
    -.8005037324172778E+00,0.8005037324172778E+00, 
    -.9786316733709013E+00,0.9786316733709015E+00, 
    -.1931251782159229E+00,0.1931251782159230E+00, 
    0.7448428411561702E+00,-.7448428411561699E+00, 
    0.9910640341773119E+00,-.9910640341773119E+00, 
    0.4487876769216255E+00,-.4487876769216254E+00, 
    0.7730731265491211E+00,-.7730731265491211E+00, 
    -.8538350485874476E-02,0.8538350485874577E-02, 
    0.9305948052973293E+00,-.9305948052973293E+00, 
    -.4410686164106586E+00,0.4410686164106588E+00, 
    0.5726168020738448E+00,-.5726168020738448E+00, 
    0.9094921563319007E+00,-.9094921563319007E+00, 
    0.6650799110223121E+00,-.6650799110223121E+00, 
    0.9845379327215906E+00,-.9845379327215906E+00, 
    -.1489897109425285E-01,0.1489897109425291E-01, 
    0.4777858927462423E+00,-.4777858927462423E+00, 
    0.2346150739656185E-01,-.2346150739656185E-01, 
    -.6464750210052006E+00,0.6464750210052006E+00, 
    -.3898509033065152E+00,0.3898509033065153E+00, 
    -.2596817297614839E+00,0.2596817297614840E+00, 
    0.8930982971146899E-01,-.8930982971146888E-01, 
    0.2378419730518350E+00,-.2378419730518349E+00, 
    -.4923226235686091E+00,0.4923226235686091E+00, 
    0.3260553152624549E+00,-.3260553152624548E+00, 
    0.6900640280303905E+00,-.6900640280303905E+00, 
    0.8399949350854392E+00,-.8399949350854392E+00, 
    0.2175201355296100E+00,-.2175201355296099E+00, 
    -.2552655348509120E+00,0.2552655348509120E+00, 
    0.8967910069881992E+00,-.8967910069881992E+00 };
  static double ys[] = {
    -.8684124523298049E+00,0.8684124523298051E+00, 
    -.9445108390325675E+00,0.9445108390325677E+00, 
    0.8224833612237958E+00,-.8224833612237956E+00, 
    -.6674381398519379E+00,0.6674381398519381E+00, 
    -.9903377411510502E+00,0.9903377411510502E+00, 
    -.9907543940322484E+00,0.9907543940322484E+00, 
    0.4473528589231669E-01,-.4473528589231656E-01, 
    -.9343866955469575E+00,0.9343866955469575E+00, 
    0.4520870573336636E+00,-.4520870573336634E+00, 
    -.9274277947725086E+00,0.9274277947725084E+00, 
    -.9866712872285782E+00,0.9866712872285782E+00, 
    -.4550564012095777E+00,0.4550564012095778E+00, 
    -.8127883212479330E+00,0.8127883212479330E+00, 
    -.9813307175904386E+00,0.9813307175904383E+00, 
    -.9279142876740704E+00,0.9279142876740704E+00, 
    -.9320431789499202E+00,0.9320431789499202E+00, 
    0.4282792320995888E+00,-.4282792320995887E+00, 
    -.5362236701446004E+00,0.5362236701446004E+00, 
    -.6399087385202717E+00,0.6399087385202717E+00, 
    -.8135057649888232E+00,0.8135057649888232E+00, 
    0.6470101042020234E+00,-.6470101042020232E+00, 
    -.8155924092874611E+00,0.8155924092874611E+00, 
    -.7941248493230142E+00,0.7941248493230142E+00, 
    0.2473440793805384E+00,-.2473440793805383E+00, 
    -.3454192145795777E+00,0.3454192145795778E+00, 
    -.2668503935770928E+00,0.2668503935770929E+00, 
    -.4625176827521818E+00,0.4625176827521818E+00, 
    -.6895677594851285E-01,0.6895677594851290E-01, 
    -.2404342948898286E-01,0.2404342948898286E-01, 
    -.6484663233565804E+00,0.6484663233565804E+00, 
    -.9898736812731077E+00,0.9898736812731077E+00, 
    -.6343263441026139E+00,0.6343263441026139E+00, 
    -.9814223809013130E+00,0.9814223809013130E+00, 
    -.2765567197102583E+00,0.2765567197102583E+00, 
    -.4080894058196405E+00,0.4080894058196404E+00, 
    -.9078967063955937E+00,0.9078967063955937E+00, 
    0.1497779073446914E+00,-.1497779073446913E+00, 
    -.1279062882959751E+00,0.1279062882959752E+00, 
    -.6902255847834162E+00,0.6902255847834162E+00, 
    -.2059314662901914E+00,0.2059314662901914E+00, 
    -.8238753864861620E+00,0.8238753864861620E+00 };
  static double ws[] = {
    0.2118008413970087E-02,0.2118008413970087E-02, 
    0.4561991935236101E-02,0.4561991935236101E-02, 
    0.6876351365235698E-02,0.6876351365235698E-02, 
    0.1124466702102733E-01,0.1124466702102733E-01, 
    0.4818715192537238E-02,0.4818715192537238E-02, 
    0.5989564896092799E-02,0.5989564896092799E-02, 
    0.1767795981651036E-01,0.1767795981651036E-01, 
    0.2380089975465993E-01,0.2380089975465993E-01, 
    0.4148355130225962E-01,0.4148355130225962E-01, 
    0.1184358597738201E-01,0.1184358597738201E-01, 
    0.1089514188613057E-01,0.1089514188613057E-01, 
    0.3129460714416973E-01,0.3129460714416973E-01, 
    0.2951898806186879E-01,0.2951898806186879E-01, 
    0.2370685766690177E-02,0.2370685766690177E-02, 
    0.2822396585947564E-01,0.2822396585947564E-01, 
    0.2219335711268346E-01,0.2219335711268346E-01, 
    0.9583664348129897E-02,0.9583664348129897E-02, 
    0.6372673897970250E-01,0.6372673897970250E-01, 
    0.4375827379330900E-01,0.4375827379330900E-01, 
    0.4429543031127132E-01,0.4429543031127132E-01, 
    0.2471410210516637E-01,0.2471410210516637E-01, 
    0.4249825844225831E-01,0.4249825844225831E-01, 
    0.4558876680953739E-01,0.4558876680953739E-01, 
    0.2756490338260340E-01,0.2756490338260340E-01, 
    0.6243396099268773E-01,0.6243396099268773E-01, 
    0.1168182946833649E-01,0.1168182946833649E-01, 
    0.7382665975220684E-01,0.7382665975220684E-01, 
    0.7666679963909392E-01,0.7666679963909392E-01, 
    0.4391713647273730E-01,0.4391713647273730E-01, 
    0.4705013653249818E-01,0.4705013653249818E-01, 
    0.9574427334796009E-02,0.9574427334796009E-02, 
    0.6262807887841917E-01,0.6262807887841917E-01, 
    0.1507382504443397E-01,0.1507382504443397E-01, 
    0.7853972113541645E-01,0.7853972113541645E-01, 
    0.6890694586576196E-01,0.6890694586576196E-01, 
    0.3574201277291281E-01,0.3574201277291281E-01, 
    0.6408284590143187E-01,0.6408284590143187E-01, 
    0.4908916812527610E-01,0.4908916812527610E-01, 
    0.5502691558122864E-01,0.5502691558122864E-01, 
    0.8126618865329004E-01,0.8126618865329004E-01, 
    0.2206473054465979E-01,0.2206473054465979E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule22 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE22 returns the rule of degree 22.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.1168773459913643E-01,-.7371881387645758E+00, 
    0.8294952934036814E+00,-.5954751704041336E+00, 
    -.4324369794733402E+00,-.7688482950823662E+00, 
    0.9204074310294709E+00,0.9345989584398445E-01, 
    -.8936155650763420E+00,-.7702459791051814E+00, 
    -.1373629144933700E+00,0.5509811434354096E+00, 
    -.1491245367450688E+00,0.3561825582103059E+00, 
    0.9353031268387029E+00,-.2924575431987692E+00, 
    0.2443965863661181E+00,-.6390847658587071E+00, 
    -.4984730279242101E+00,0.7392293471332223E+00, 
    -.9874262069353944E+00,-.9183316919543670E-01, 
    -.7612182335622527E+00,-.8617096626977842E-01, 
    0.9888069308129792E+00,-.5588987313328233E+00, 
    -.5751189401814027E+00,-.9225108518720977E+00, 
    0.8984615209750768E+00,0.8382620472105732E+00, 
    0.1626926025129410E+00,-.9967788148796793E+00, 
    -.6486037112886791E+00,-.5668454455159495E+00, 
    0.9856879280293260E+00,0.9717573267273887E+00, 
    -.3076008058514624E+00,0.3587016306051349E+00, 
    -.7623612106456181E+00,0.6388680957517785E+00, 
    0.5510537420201546E+00,0.1053488710493425E+00, 
    0.9916832391481870E+00,-.8977327473577641E+00, 
    0.4016301282586624E+00,-.9838441232429354E+00, 
    0.7905399896596540E+00,-.9824504337674188E+00, 
    -.5775436902276508E+00,-.9772275458992810E+00, 
    -.8927278831009643E+00,0.6072501413880697E+00, 
    0.3391264841215570E+00,0.8802736916527715E+00, 
    -.3614572700854332E+00,0.7628120974888042E+00, 
    0.9410464685954009E+00,0.5766667841609323E+00, 
    0.3785886954980848E+00,-.7336046871225165E+00, 
    -.3595240670259577E+00,0.1395498733753600E+00, 
    -.8379800610987231E+00,0.9026858439122091E+00, 
    -.7544128178700894E+00,-.1371254063139990E+00, 
    0.7326981395585584E+00,0.1008288194176049E+00, 
    -.1402954353279801E+00,0.5600957930752593E+00, 
    -.9835815408757524E+00,0.9794922442499654E+00, 
    -.3422476201633468E+00,0.8847947144492797E+00, 
    0.3480028949057162E+00,0.1399751797907734E+00, 
    0.7410269267329871E+00,-.8707213419949226E+00, 
    -.7829703967945006E+00,-.7518531217472861E-01, 
    0.1797607493737034E+00,0.5514938721272684E+00, 
    0.9189774855743247E+00,-.9478502721386363E+00, 
    -.3753442057328584E+00,0.9824594434447735E+00, 
    0.7794560599617187E+00,0.5868549502954560E+00, 
    -.9830132030284936E+00,-.9109869789024120E+00, 
    -.3369368343786870E+00,-.9656919342504539E+00, 
    0.9878086002165142E+00 };
  static double ys[] = {
    0.2190948368619670E+00,-.5638633812668349E+00, 
    0.7772605528528743E+00,-.2566398183828403E+00, 
    -.5129047910078139E+00,-.4506194723257167E+00, 
    0.8425861006849568E+00,-.4518983245417754E+00, 
    0.3543337144777458E+00,0.4959761932603439E+00, 
    0.2485911682035064E+00,-.3132773446075675E+00, 
    -.9937282156984992E+00,-.2003997888710157E+00, 
    -.5607980017506478E+00,-.3522379616257777E+00, 
    -.4426431371959611E+00,-.3692415241854752E-01, 
    0.7496615691984209E+00,0.6926590385285363E+00, 
    -.1850687765410226E+00,-.1925693267588198E+00, 
    0.1363788459324560E+00,0.9515651662120697E+00, 
    0.5579343853339161E-01,0.9517652105508706E+00, 
    -.9389370283909926E+00,-.3771353916795439E+00, 
    -.9947202210423953E+00,-.3609118075153294E+00, 
    0.9902936741983037E+00,0.3133668776125615E+00, 
    0.6199365920821208E+00,0.2873358944700919E+00, 
    -.7434073625585887E+00,0.3666474918307532E+00, 
    0.8685050837394444E+00,0.9342984101597417E+00, 
    -.8496275328566278E+00,-.4918800631346358E+00, 
    0.5207829647449653E+00,0.4875959314271490E+00, 
    -.3759522691536638E+00,-.9452616179455003E+00, 
    -.6399489008638399E+00,-.5698761741539897E+00, 
    -.6821182736060648E+00,-.9838318115418717E+00, 
    -.7110797215818025E+00,-.8424531888684956E+00, 
    -.7019404243765398E+00,-.8133442851406898E+00, 
    0.6937136042433247E+00,0.1403628667969644E+00, 
    0.4853142503301350E+00,-.9404688641752270E+00, 
    -.1594253939471121E+00,0.8410558281897980E+00, 
    -.9090435768164216E+00,0.8618858293632953E+00, 
    -.9801201774367985E+00,0.2925291147308246E-01, 
    -.1389383474162093E+00,-.8488076014783108E+00, 
    -.9950422147374585E+00,0.6774278489603331E+00, 
    0.3281092020661301E+00,0.8336091093259060E+00, 
    -.6447342914933101E+00,-.9855619048067824E+00, 
    0.9863388067410731E+00,-.9474857272173681E+00, 
    -.8291238758195407E+00,0.5386555528488844E+00, 
    0.2895743833834052E+00,-.7848534988754099E+00, 
    -.1043902046106473E+00,0.7337865418484663E+00, 
    0.9901322140668845E+00,-.9181369104835962E+00, 
    -.9807391493006852E+00,0.7602577556121397E-01, 
    0.9902833594797992E+00,0.7639725935225382E-01, 
    0.1089930212643242E-01,0.6870529955030933E+00, 
    0.9401212847078045E+00,0.9891928053902059E+00, 
    0.8336461919417253E+00,0.9356252744608189E+00, 
    0.9934890766781104E+00,0.5866648232658751E+00, 
    0.9371299738758988E+00 };
  static double ws[] = {
    0.6067121956035168E-02,0.2527104513859865E-01, 
    0.9025528955311097E-02,0.5803636803454140E-01, 
    0.4942739001800841E-01,0.2794072665175644E-01, 
    0.1603496923754518E-01,0.6013880064040703E-01, 
    0.2941264773858220E-01,0.3334592290253230E-01, 
    0.7897021060987033E-01,0.3685303544908184E-01, 
    0.3948104138816246E-02,0.6679001222692597E-01, 
    0.2111954975158539E-01,0.4792141258829016E-01, 
    0.2695976162319685E-01,0.2470945430647790E-01, 
    0.3947338905928970E-01,0.3294627756714954E-01, 
    0.9846011394782925E-02,0.7297001782497028E-01, 
    0.4106989048022942E-01,0.2161795857158166E-01, 
    0.8365413396512252E-02,0.1967362232030983E-01, 
    0.2188213465301344E-01,0.2494392526364261E-01, 
    0.3519731383753954E-02,0.3478325192787569E-01, 
    0.8190571090769309E-02,0.6108434551087310E-02, 
    0.3434968062753964E-01,0.6007998437323292E-01, 
    0.7504666111448662E-02,0.1603184108436172E-01, 
    0.3437516964723778E-01,0.2571624349884583E-01, 
    0.2694235332130491E-01,0.4466450734173113E-01, 
    0.5500707054662403E-01,0.7222938823304295E-01, 
    0.7486350326225666E-02,0.1143517520249437E-01, 
    0.5817171888321045E-01,0.9994978126228900E-02, 
    0.3467709360229788E-01,0.1964349392895841E-02, 
    0.4240321318174279E-01,0.8728495324106815E-02, 
    0.2591557111759762E-01,0.3720861411904331E-01, 
    0.5505688737445376E-01,0.3631603425342831E-01, 
    0.6567401226660792E-01,0.1795583567053723E-01, 
    0.2289367999547543E-01,0.3547488447682710E-01, 
    0.2995254108777568E-01,0.2686044680493277E-01, 
    0.1079041750322994E-01,0.7877736606773408E-01, 
    0.3733692480386966E-01,0.1821095142444074E-01, 
    0.5325368320113208E-02,0.5922160899036626E-01, 
    0.5076619476679856E-01,0.4372734278286988E-01, 
    0.6332661510009212E-01,0.9790273018996520E-02, 
    0.1776555836345017E-02,0.4813300061774173E-02, 
    0.4258882597057490E-01,0.3052111793035846E-01, 
    0.7003826037949971E-01,0.5120830680793959E-01, 
    0.4944440063240135E-01,0.2605683245528804E-01, 
    0.6201489406026509E-02,0.3089929178556832E-01, 
    0.1316878814412241E-01,0.6382065038385205E-01, 
    0.3990900163624796E-02,0.2334151205865433E-01, 
    0.7860422025755921E-01,0.1000548940130294E-01, 
    0.1798134856368640E-01,0.8411973937140507E-02, 
    0.7382035812019881E-02,0.1170576178999751E-01, 
    0.6920800809156403E-02,0.1595811031270932E-01, 
    0.3880611624295888E-02 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule23 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE23 returns the rule of degree 23.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.1922689542773817E+00,-.1922689542773816E+00, 
    -.9880860844879094E+00,0.9880860844879096E+00, 
    0.4588599108929147E+00,-.4588599108929146E+00, 
    0.2847661416086631E+00,-.2847661416086630E+00, 
    0.8023730521007825E+00,-.8023730521007822E+00, 
    -.4008300617025742E+00,0.4008300617025743E+00, 
    0.9128169992246806E+00,-.9128169992246804E+00, 
    0.6402796135160949E+00,-.6402796135160946E+00, 
    0.6316970739891037E+00,-.6316970739891037E+00, 
    -.5737921022421535E+00,0.5737921022421535E+00, 
    0.4346256934364524E+00,-.4346256934364522E+00, 
    0.7764676048715878E+00,-.7764676048715878E+00, 
    0.7244917270725689E+00,-.7244917270725689E+00, 
    -.1888739484686455E+00,0.1888739484686456E+00, 
    0.9943536134460622E+00,-.9943536134460622E+00, 
    -.7313948142734760E+00,0.7313948142734760E+00, 
    0.8943097577057932E+00,-.8943097577057932E+00, 
    0.8436479267198643E+00,-.8436479267198643E+00, 
    0.7791831851449014E+00,-.7791831851449014E+00, 
    0.8641774713700604E+00,-.8641774713700604E+00, 
    0.9617690298704565E+00,-.9617690298704563E+00, 
    -.4195385320317492E+00,0.4195385320317493E+00, 
    0.2226040785176633E+00,-.2226040785176632E+00, 
    0.9898168214058006E+00,-.9898168214058006E+00, 
    -.7489474270848241E+00,0.7489474270848243E+00, 
    0.1917103338448669E-01,-.1917103338448658E-01, 
    0.9493597561722253E+00,-.9493597561722253E+00, 
    0.9904394141829281E+00,-.9904394141829279E+00, 
    -.5829692270288439E+00,0.5829692270288439E+00, 
    -.6147395289066456E+00,0.6147395289066458E+00, 
    0.9535616627723794E+00,-.9535616627723796E+00, 
    -.2159492946526814E+00,0.2159492946526815E+00, 
    0.6262127585905694E+00,-.6262127585905694E+00, 
    0.9907200353310519E+00,-.9907200353310519E+00, 
    0.9843113976598121E+00,-.9843113976598121E+00, 
    -.8786718625785459E+00,0.8786718625785461E+00, 
    0.9324996820585237E+00,-.9324996820585237E+00, 
    0.8556365285580105E+00,-.8556365285580105E+00, 
    0.6098296498831587E+00,-.6098296498831587E+00, 
    0.2401978383668566E-02,-.2401978383668536E-02, 
    -.2160362298752296E+00,0.2160362298752297E+00, 
    0.9419047852611938E+00,-.9419047852611938E+00, 
    0.9784325270256451E+00,-.9784325270256451E+00, 
    0.4245835794119392E+00,-.4245835794119392E+00, 
    0.2180806198402934E+00,-.2180806198402933E+00, 
    0.4211715480246619E+00,-.4211715480246619E+00, 
    0.2159184187190165E+00,-.2159184187190165E+00, 
    0.9145012870119151E+00,-.9145012870119151E+00, 
    -.6659715296954374E-02,0.6659715296954456E-02 };
  static double ys[] = {
    -.8413204499428454E+00,0.8413204499428454E+00, 
    -.9814588472800125E+00,0.9814588472800123E+00, 
    -.6684346870698402E+00,0.6684346870698402E+00, 
    -.7992729842209235E+00,0.7992729842209235E+00, 
    -.9416742417176776E+00,0.9416742417176776E+00, 
    -.6678440121187212E+00,0.6678440121187212E+00, 
    -.9891421484914774E+00,0.9891421484914776E+00, 
    -.9915720334011255E+00,0.9915720334011255E+00, 
    -.4898989236348730E+00,0.4898989236348731E+00, 
    -.4662751769077947E+00,0.4662751769077946E+00, 
    -.9416451658784137E+00,0.9416451658784137E+00, 
    -.2882717499568009E+00,0.2882717499568010E+00, 
    0.2486649083667506E+00,-.2486649083667505E+00, 
    -.9882692965051619E+00,0.9882692965051619E+00, 
    0.7885280263196212E+00,-.7885280263196209E+00, 
    -.6589298745979582E+00,0.6589298745979582E+00, 
    -.5251907426680982E+00,0.5251907426680982E+00, 
    0.2698952252390472E-01,-.2698952252390462E-01, 
    -.7058148204376444E+00,0.7058148204376444E+00, 
    0.8039684450793269E+00,-.8039684450793269E+00, 
    -.9795669903911515E+00,0.9795669903911517E+00, 
    -.9436139885736695E+00,0.9436139885736695E+00, 
    -.9897986315498021E+00,0.9897986315498021E+00, 
    -.1515697280553813E-02,0.1515697280553934E-02, 
    -.9288927030420643E+00,0.9288927030420643E+00, 
    -.9319407769628315E+00,0.9319407769628315E+00, 
    0.6363710989431733E+00,-.6363710989431731E+00, 
    -.9179467622904829E+00,0.9179467622904831E+00, 
    -.8263113277740309E+00,0.8263113277740309E+00, 
    -.9928811295871781E+00,0.9928811295871781E+00, 
    0.9091589955140115E+00,-.9091589955140112E+00, 
    -.4548720285260120E+00,0.4548720285260120E+00, 
    -.8475542857772346E+00,0.8475542857772346E+00, 
    0.4247422802343108E+00,-.4247422802343107E+00, 
    -.3947150748420805E+00,0.3947150748420806E+00, 
    -.9832744227413978E+00,0.9832744227413978E+00, 
    -.1923063082978693E+00,0.1923063082978694E+00, 
    0.4580266855460115E+00,-.4580266855460114E+00, 
    -.4698541060861779E-01,0.4698541060861786E-01, 
    -.2437626817883838E+00,0.2437626817883838E+00, 
    -.8348003321715464E+00,0.8348003321715464E+00, 
    0.2268917425003757E+00,-.2268917425003756E+00, 
    -.6927513253690014E+00,0.6927513253690016E+00, 
    -.2736597506928621E+00,0.2736597506928621E+00, 
    -.4885665960318681E+00,0.4885665960318681E+00, 
    0.2093478053463650E+00,-.2093478053463650E+00, 
    -.1584479979856305E-01,0.1584479979856307E-01, 
    -.8398490539050634E+00,0.8398490539050636E+00, 
    -.6725993072235357E+00,0.6725993072235357E+00 };
  static double ws[] = {
    0.2087996398690324E-01,0.2087996398690324E-01, 
    0.1454948755827531E-02,0.1454948755827531E-02, 
    0.4748321475897626E-01,0.4748321475897626E-01, 
    0.2581005409389433E-01,0.2581005409389433E-01, 
    0.1501004260111713E-01,0.1501004260111713E-01, 
    0.4842495146008856E-01,0.4842495146008856E-01, 
    0.2506826413864942E-02,0.2506826413864942E-02, 
    0.6535276655159778E-02,0.6535276655159778E-02, 
    0.4894280385061923E-01,0.4894280385061923E-01, 
    0.4975905503606815E-01,0.4975905503606815E-01, 
    0.2214413184536724E-01,0.2214413184536724E-01, 
    0.4313124475180520E-01,0.4313124475180520E-01, 
    0.4433897643935200E-01,0.4433897643935200E-01, 
    0.9046276329783944E-02,0.9046276329783944E-02, 
    0.3911096929389304E-02,0.3911096929389304E-02, 
    0.3599231426419872E-01,0.3599231426419872E-01, 
    0.2838276711743061E-01,0.2838276711743061E-01, 
    0.3373992192798946E-01,0.3373992192798946E-01, 
    0.3336454855497337E-01,0.3336454855497337E-01, 
    0.2095059245259030E-01,0.2095059245259030E-01, 
    0.2046044703097752E-02,0.2046044703097752E-02, 
    0.2164205919445887E-01,0.2164205919445887E-01, 
    0.8285628342573098E-02,0.8285628342573098E-02, 
    0.8397251850551729E-02,0.8397251850551729E-02, 
    0.1731140501927820E-01,0.1731140501927820E-01, 
    0.2467601783099452E-01,0.2467601783099452E-01, 
    0.1682157937509966E-01,0.1682157937509966E-01, 
    0.3331814718792869E-02,0.3331814718792869E-02, 
    0.3259337068821334E-01,0.3259337068821334E-01, 
    0.6005493842446140E-02,0.6005493842446140E-02, 
    0.8210801440654199E-02,0.8210801440654199E-02, 
    0.6548837153964436E-01,0.6548837153964436E-01, 
    0.3088139176161778E-01,0.3088139176161778E-01, 
    0.7183392776926659E-02,0.7183392776926659E-02, 
    0.9748230121097153E-02,0.9748230121097153E-02, 
    0.5383514720217621E-02,0.5383514720217621E-02, 
    0.2224440941283463E-01,0.2224440941283463E-01, 
    0.3213763319312607E-01,0.3213763319312607E-01, 
    0.5959500864875663E-01,0.5959500864875663E-01, 
    0.7353356562572214E-01,0.7353356562572214E-01, 
    0.3920268743920687E-01,0.3920268743920687E-01, 
    0.2191097090277860E-01,0.2191097090277860E-01, 
    0.1035031173873475E-01,0.1035031173873475E-01, 
    0.6659548623039759E-01,0.6659548623039759E-01, 
    0.6553251667930458E-01,0.6553251667930458E-01, 
    0.6644084983327379E-01,0.6644084983327379E-01, 
    0.7384544963870446E-01,0.7384544963870446E-01, 
    0.1649581640678158E-01,0.1649581640678158E-01, 
    0.5651348047241081E-01,0.5651348047241081E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule24 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE24 returns the rule of degree 24.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.9947148928250171E+00,0.9910990446809216E+00, 
    0.9693866534373358E+00,0.8366669887301708E+00, 
    0.9976180969133432E+00,-.9626395872771291E+00, 
    -.9956012694874316E+00,0.9827883755011416E+00, 
    0.9899696481565038E+00,-.1027795858264751E+00, 
    -.5372013809570906E+00,-.1302736589054738E+00, 
    -.9966597496356319E+00,-.8416889461091814E+00, 
    -.9888026187105032E+00,0.5314747340473951E+00, 
    -.8162827902273060E+00,0.9739798586001176E+00, 
    -.9411304883795972E+00,-.9821743630414852E+00, 
    -.8146483736791913E+00,0.2216111751930999E+00, 
    0.9919273511735743E+00,0.8515268909036164E+00, 
    -.6344208939912110E-01,0.8124697219833941E+00, 
    -.9159342331978603E+00,0.1939576002609792E+00, 
    0.9879450437516667E+00,0.2535126040886103E+00, 
    -.9628662071921108E+00,0.8784885445188073E+00, 
    -.8980693800663294E+00,-.4876187319409098E+00, 
    0.8244493216582833E+00,-.7771226957924775E-01, 
    0.9133806640305893E+00,-.2863497657037777E+00, 
    -.7746256411943806E+00,0.8788745310361548E+00, 
    -.9610911899226917E+00,0.9600320888811048E+00, 
    0.9425138210931967E+00,0.7626320556994604E+00, 
    -.3308772778351785E+00,0.1279919952613696E+00, 
    0.5294732769614892E+00,0.3555627902477361E+00, 
    -.9279977222115954E+00,0.9931449449560651E+00, 
    0.5102490315575532E-01,-.7525118688957001E+00, 
    -.4719724839934526E+00,0.1520260254758981E+00, 
    0.9489604126336816E+00,-.2169062988895815E+00, 
    -.8687157102597637E+00,-.7360447596642431E+00, 
    0.8497144016937899E+00,-.6403949304059380E+00, 
    0.9160879280172893E+00,0.5974122417304695E+00, 
    0.9356564573866695E+00,0.2602110178808541E-02, 
    0.2185370159809754E+00,-.4174601636720920E+00, 
    -.6704915695516799E+00,0.7098684394657311E+00, 
    -.6662903518807954E+00,0.6842750704411169E+00, 
    0.6940401495592319E+00,0.2055436622280245E+00, 
    0.5409016726579078E+00,-.2111056222157408E+00, 
    0.4249715565953610E+00,-.4579852398972861E+00, 
    -.4762834069005017E+00,-.3112968492969387E+00, 
    -.9861402159543649E+00,0.5994434867749180E+00, 
    -.6329536320758256E+00,0.4071516862564607E+00, 
    -.5899049368040069E+00,-.8873886456704851E+00, 
    -.9429884035167735E+00,-.8103411028772554E+00, 
    -.6172197178932676E+00,0.4403002993242678E+00, 
    -.5964555068203565E+00,0.4214363261995724E+00, 
    -.2167017340827659E+00,0.3600751655325609E+00, 
    0.9573652308066442E-02,-.1986519498186896E+00, 
    0.5940360363603530E+00,0.7449094683537840E+00, 
    0.3446198859685110E+00,0.7550496004921502E+00, 
    0.2034447821293972E+00,-.8764429892755836E+00, 
    -.9821148921468710E-02,-.4283906719131795E+00, 
    -.2565291139819166E+00,-.4253124084387095E+00, 
    -.4107973035474266E+00,0.6107899647074595E+00, 
    -.7740269973316680E+00,-.9874686104628365E+00, 
    -.1780069523272345E-01 };
  static double ys[] = {
    0.9921816068213789E+00,0.9927512676582938E+00, 
    -.9807210913639813E+00,-.9849335979706825E+00, 
    -.9031203174435843E+00,-.9807185407289059E+00, 
    0.1826280218638680E+00,0.1985977578306220E-03, 
    -.1343242025779884E+00,0.9836572598301550E+00, 
    -.9945228775385402E+00,-.9932466930618153E+00, 
    -.9021729201562282E+00,0.9919450680865269E+00, 
    -.5469411165135486E+00,-.9946152107188805E+00, 
    -.9833999582719447E+00,-.7694911408818028E+00, 
    0.9454206624475765E+00,0.5820485751363867E+00, 
    0.8726697552774600E+00,0.9884035519660979E+00, 
    0.3773859079404019E+00,0.9854732000616159E+00, 
    0.7430686212004588E+00,-.8172733392382462E+00, 
    0.7368090773498677E+00,0.8191066542200217E+00, 
    -.5129857572745989E+00,0.7339438859190308E+00, 
    0.3121292180404563E+00,0.4238556817241846E+00, 
    -.9023577502860924E+00,-.8635258163522489E+00, 
    -.5019428343401774E+00,-.9100031022978671E+00, 
    -.6538056836085887E+00,0.9419753943798090E+00, 
    -.8083978001656728E+00,0.7971819583886433E+00, 
    -.7612052337333153E+00,0.6275466264499783E+00, 
    0.1573227985990733E+00,0.6249214522371203E+00, 
    -.9587696584061124E+00,-.9822260622287698E+00, 
    -.8459845466706856E+00,-.6962190414609735E+00, 
    -.3797057196057647E+00,0.8208860694583049E+00, 
    0.9117911472030502E+00,0.2466478148353290E+00, 
    0.2531249337761478E+00,-.8311014632352622E+00, 
    0.9327963019854183E+00,0.3901208978319467E+00, 
    -.6258430618900800E+00,-.4854699629184363E+00, 
    -.8013355833035751E-01,0.7663811452715913E+00, 
    -.9064248594697549E+00,0.9816002629218432E+00, 
    -.3156827868701995E+00,0.5928915583696132E+00, 
    -.2758048916460264E-01,0.6590354293018823E+00, 
    -.9381955060600012E+00,-.2998249948498765E+00, 
    0.9581605160533719E+00,-.6881340269166623E+00, 
    -.9406665893940601E+00,0.3996485712287530E+00, 
    -.5101678225103823E+00,-.4421541038328301E+00, 
    0.5912997929743822E+00,0.9932623278482910E+00, 
    0.8781241163707275E+00,0.5625436230792860E+00, 
    -.1854744570391582E+00,0.7725191379684401E+00, 
    -.2272561659397919E-01,-.2598606537214366E+00, 
    0.4549553364199058E+00,0.4412559753061466E+00, 
    -.7387531735800691E-02,-.1993131251585695E+00, 
    -.7001646090032370E+00,0.8894305723468252E+00, 
    -.3651421534696686E+00,0.2021031040499760E+00, 
    -.4759573989947036E-01,0.9560701833067996E+00, 
    -.2523062776553245E+00,0.8029027500701331E+00, 
    -.3052326412673496E-01,0.9085765488533346E+00, 
    -.9429472498836997E+00,0.1985140340578119E+00, 
    -.4744209684134443E+00,0.8925996493449884E-01, 
    0.1788338183287970E+00,-.5864723938922081E+00, 
    -.7750597333889675E+00,-.2335671290163587E+00, 
    0.1465570094281399E+00,0.4191985939494091E+00, 
    0.6100181397154107E+00,0.8429972982275405E+00, 
    -.6496814943785820E+00 };
  static double ws[] = {
    0.7673500038857188E-03,0.8868846753542908E-03, 
    0.2702370475569613E-02,0.5392375480513576E-02, 
    0.1761998115149276E-02,0.3252629117203137E-02, 
    0.4483878506691741E-02,0.4387445903456329E-02, 
    0.6293360138903417E-02,0.9110559046011608E-02, 
    0.4960172913990594E-02,0.4540132762921259E-02, 
    0.2380910242803054E-02,0.4079407684566781E-02, 
    0.7372795522491610E-02,0.5764065832069020E-02, 
    0.5996283579933687E-02,0.8231405547532079E-02, 
    0.7742564602285441E-02,0.9689437679198899E-02, 
    0.1992786357523802E-01,0.7578960870331433E-02, 
    0.7297515243700454E-02,0.5560565572125132E-02, 
    0.1431910744302743E-01,0.2083181247879573E-01, 
    0.1871125538452946E-01,0.1856480159069449E-01, 
    0.7926961712192319E-02,0.3149563402986980E-01, 
    0.1220979773721510E-01,0.3064077139403503E-01, 
    0.1246839879952511E-01,0.2727439670515271E-01, 
    0.2989955062354429E-01,0.2568927037462639E-01, 
    0.1731775418487197E-01,0.1651762549497732E-01, 
    0.2348412516559215E-01,0.2103194138017257E-01, 
    0.1218315588129536E-01,0.1587611356755889E-01, 
    0.2192984739792703E-01,0.3559892184688495E-01, 
    0.1598768389647512E-01,0.9998012603521544E-02, 
    0.3014828623151700E-01,0.4445807626369159E-01, 
    0.2303305529498344E-01,0.4177457820842124E-02, 
    0.2595692066159871E-01,0.4297075322654262E-01, 
    0.2734909126669048E-01,0.3576989995633175E-01, 
    0.8104038598451695E-02,0.5896633070443507E-01, 
    0.2530733931631170E-01,0.3397799558910337E-01, 
    0.3681905878262371E-01,0.3127418749936042E-01, 
    0.1043340150186936E-01,0.9877922117697804E-02, 
    0.2268684828187055E-01,0.5274200508727377E-01, 
    0.6668503169982437E-01,0.4053161003772505E-01, 
    0.1580741394786280E-01,0.4437844031312334E-01, 
    0.1465934495793213E-01,0.3316127050631540E-01, 
    0.1633165248572775E-01,0.6234782740688132E-01, 
    0.4701994260472551E-01,0.5761066620189915E-01, 
    0.4794170534889713E-01,0.5798541312289998E-02, 
    0.2505250386803961E-01,0.2055053414266298E-01, 
    0.8995847445576632E-02,0.3472177169167036E-01, 
    0.5444084929958415E-01,0.5971645153263631E-01, 
    0.4912587830256034E-01,0.2538968252075221E-01, 
    0.1542318408388366E-01,0.4013335870948910E-01, 
    0.3298688496369413E-01,0.2322125010689848E-01, 
    0.3671677541506277E-01,0.6139982772160751E-01, 
    0.6362592843510947E-01,0.9567538562556033E-02, 
    0.6459890039693177E-01,0.3562315409338049E-01, 
    0.5516254605923790E-01,0.1998312905915403E-01, 
    0.2052118732484228E-01,0.4458997988822495E-01, 
    0.5811519536921120E-01,0.2297223481430297E-01, 
    0.6766667516677112E-01,0.4532069980219111E-01, 
    0.4089559205492352E-01,0.5328231011100253E-01, 
    0.4369211916320918E-01,0.4895075796108013E-01, 
    0.3261352596485557E-01,0.5100190068919709E-02, 
    0.5182667283735970E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule25 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE25 returns the rule of degree 25.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.7103850079486462E+00,0.7103850079486462E+00, 
    0.2350136833954603E+00,-.2350136833954603E+00, 
    0.2319467644431913E+00,-.2319467644431913E+00, 
    0.1887551793322332E-01,-.1887551793322329E-01, 
    0.4762227841581165E+00,-.4762227841581164E+00, 
    0.9223578267424110E+00,-.9223578267424110E+00, 
    0.8612991449140567E+00,-.8612991449140567E+00, 
    -.2008507367686636E+00,0.2008507367686637E+00, 
    0.9910315414864325E+00,-.9910315414864325E+00, 
    0.3217353529776810E+00,-.3217353529776809E+00, 
    0.8933439162316077E+00,-.8933439162316077E+00, 
    -.5783268726892828E+00,0.5783268726892828E+00, 
    0.1876981715173633E+00,-.1876981715173632E+00, 
    -.4040937317707446E+00,0.4040937317707447E+00, 
    -.8086692861933220E+00,0.8086692861933220E+00, 
    0.7699071344451243E+00,-.7699071344451243E+00, 
    0.1349449671364110E+00,-.1349449671364109E+00, 
    0.7067483874126944E+00,-.7067483874126944E+00, 
    0.9941851173182269E+00,-.9941851173182271E+00, 
    0.9823116186120971E+00,-.9823116186120969E+00, 
    0.4238538814747814E-01,-.4238538814747803E-01, 
    -.1196214241956023E+00,0.1196214241956024E+00, 
    0.6534777415625678E+00,-.6534777415625675E+00, 
    0.9750429947957491E+00,-.9750429947957491E+00, 
    0.8117686312525225E+00,-.8117686312525225E+00, 
    0.6171167577694838E+00,-.6171167577694838E+00, 
    0.5721736374948944E+00,-.5721736374948944E+00, 
    -.3415910753286391E+00,0.3415910753286392E+00, 
    -.9491077472873352E+00,0.9491077472873354E+00, 
    0.9925898844290307E+00,-.9925898844290307E+00, 
    -.2361386817708837E+00,0.2361386817708838E+00, 
    0.4376492710177434E+00,-.4376492710177433E+00, 
    0.9907335874081320E+00,-.9907335874081320E+00, 
    0.3690348754186010E+00,-.3690348754186009E+00, 
    0.2356776284326877E+00,-.2356776284326877E+00, 
    0.9290327147765426E+00,-.9290327147765426E+00, 
    0.5315931005621322E+00,-.5315931005621322E+00, 
    -.7970927883912756E+00,0.7970927883912758E+00, 
    -.6382950957253196E+00,0.6382950957253198E+00, 
    0.9125245471252791E+00,-.9125245471252793E+00, 
    -.1417325966030156E+00,0.1417325966030157E+00, 
    0.4628221381215424E+00,-.4628221381215424E+00, 
    0.9500761399436652E+00,-.9500761399436652E+00, 
    0.6585082411431559E+00,-.6585082411431559E+00, 
    -.5874365186841102E-02,0.5874365186841170E-02, 
    0.9830805981225823E+00,-.9830805981225823E+00, 
    0.9607858073894220E+00,-.9607858073894220E+00, 
    0.8126469869288578E+00,-.8126469869288578E+00, 
    0.9197969485441327E+00,-.9197969485441327E+00, 
    -.5479961396095706E+00,0.5479961396095708E+00, 
    -.4470049599321014E+00,0.4470049599321015E+00, 
    0.7873844821702739E+00,-.7873844821702737E+00, 
    0.3012570662899303E-30,0.4368201344745598E+00, 
    -.4368201344745598E+00,0.8506134670179059E+00, 
    -.8506134670179059E+00,0.9113188530371741E+00, 
    -.9113188530371739E+00,0.9836230570941239E+00, 
    -.9836230570941239E+00,0.7402837708508578E+00, 
    -.7402837708508578E+00 };
  static double ys[] = {
    -.7506170245096466E+00,0.7506170245096466E+00, 
    -.6015670998122035E-01,0.6015670998122038E-01, 
    0.6313639161925437E-01,-.6313639161925434E-01, 
    -.2037217515148590E+00,0.2037217515148590E+00, 
    -.9891004866173655E+00,0.9891004866173655E+00, 
    -.3246552383950583E+00,0.3246552383950584E+00, 
    -.1540352101742460E+00,0.1540352101742461E+00, 
    -.3377511016941964E+00,0.3377511016941964E+00, 
    -.6383724442924701E+00,0.6383724442924703E+00, 
    -.8954389923107672E+00,0.8954389923107672E+00, 
    0.6633371340392643E+00,-.6633371340392643E+00, 
    -.6441606686218248E+00,0.6441606686218248E+00, 
    -.9825525591494771E+00,0.9825525591494771E+00, 
    -.4921744377774439E+00,0.4921744377774439E+00, 
    -.8285906504912744E+00,0.8285906504912744E+00, 
    0.5016460262480213E+00,-.5016460262480213E+00, 
    -.7648310056858745E+00,0.7648310056858745E+00, 
    -.7008712571932789E+00,0.7008712571932789E+00, 
    0.9429421827689296E+00,-.9429421827689294E+00, 
    -.9889464727951927E+00,0.9889464727951929E+00, 
    -.9395791578550448E+00,0.9395791578550448E+00, 
    -.8631057904957280E+00,0.8631057904957280E+00, 
    -.9411634392706579E+00,0.9411634392706579E+00, 
    0.7901512634507242E+00,-.7901512634507240E+00, 
    -.8620663504329279E+00,0.8620663504329279E+00, 
    0.3181861315943668E+00,-.3181861315943667E+00, 
    -.4998028335682857E+00,0.4998028335682858E+00, 
    -.9494727492703661E+00,0.9494727492703661E+00, 
    -.9900142583759359E+00,0.9900142583759357E+00, 
    -.2005292720209485E+00,0.2005292720209486E+00, 
    -.6936010697328517E+00,0.6936010697328517E+00, 
    -.9618522647938675E+00,0.9618522647938675E+00, 
    0.2336705239348026E+00,-.2336705239348024E+00, 
    -.6422765133729448E+00,0.6422765133729448E+00, 
    -.4058616631814986E+00,0.4058616631814986E+00, 
    -.7553494223525696E+00,0.7553494223525699E+00, 
    -.8131056508171524E+00,0.8131056508171524E+00, 
    -.9753835609823820E+00,0.9753835609823820E+00, 
    -.9152847871689571E+00,0.9152847871689571E+00, 
    0.9075396692621631E+00,-.9075396692621629E+00, 
    -.9916640269515081E+00,0.9916640269515081E+00, 
    -.2319599091440556E+00,0.2319599091440556E+00, 
    0.2744706063362164E-01,-.2744706063362152E-01, 
    -.4230837956698025E-01,0.4230837956698034E-01, 
    -.5553447391265243E+00,0.5553447391265243E+00, 
    0.5512083775806552E+00,-.5512083775806550E+00, 
    -.4486213975897800E+00,0.4486213975897801E+00, 
    0.1653383270128847E+00,-.1653383270128846E+00, 
    0.3695711955638312E+00,-.3695711955638311E+00, 
    -.9920191871388367E+00,0.9920191871388367E+00, 
    -.8211053506315364E+00,0.8211053506315364E+00, 
    -.9936569020855669E+00,0.9936569020855669E+00, 
    -.1899153506582438E-30,0.1367147397139465E+00, 
    -.1367147397139464E+00,-.5673411356087702E+00, 
    0.5673411356087702E+00,-.9545649683609251E+00, 
    0.9545649683609253E+00,-.8813580692740732E+00, 
    0.8813580692740735E+00,-.3448010185338413E+00, 
    0.3448010185338414E+00 };
  static double ws[] = {
    0.1890289035282082E-01,0.1890289035282082E-01, 
    0.4408519139496046E-01,0.4408519139496046E-01, 
    0.1760012528851595E-01,0.1760012528851595E-01, 
    0.5135792991631447E-01,0.5135792991631447E-01, 
    0.5903718396467129E-02,0.5903718396467129E-02, 
    0.1169110569868906E-01,0.1169110569868906E-01, 
    0.2797734321676553E-01,0.2797734321676553E-01, 
    0.5430655753317624E-01,0.5430655753317624E-01, 
    0.5579615408482969E-02,0.5579615408482969E-02, 
    0.2529071368793077E-01,0.2529071368793077E-01, 
    0.2144063823465160E-01,0.2144063823465160E-01, 
    0.3501744210151039E-01,0.3501744210151039E-01, 
    0.8405006765794200E-02,0.8405006765794200E-02, 
    0.4833682312227475E-01,0.4833682312227475E-01, 
    0.1594545160823835E-01,0.1594545160823835E-01, 
    0.3558939192260235E-01,0.3558939192260235E-01, 
    0.4226984758231494E-01,0.4226984758231494E-01, 
    0.3140501945419725E-01,0.3140501945419725E-01, 
    0.2024434452721722E-02,0.2024434452721722E-02, 
    0.1349622146875648E-02,0.1349622146875648E-02, 
    0.1534182110464278E-01,0.1534182110464278E-01, 
    0.2952379758820060E-01,0.2952379758820060E-01, 
    0.1647236890585741E-01,0.1647236890585741E-01, 
    0.8590118628187588E-02,0.8590118628187588E-02, 
    0.1949659427579714E-01,0.1949659427579714E-01, 
    0.4728589565124858E-01,0.4728589565124858E-01, 
    0.4410722929300222E-01,0.4410722929300222E-01, 
    0.1872860198784019E-01,0.1872860198784019E-01, 
    0.2411335871046539E-02,0.2411335871046539E-02, 
    0.6741232108638781E-02,0.6741232108638781E-02, 
    0.4411511228498587E-01,0.4411511228498587E-01, 
    0.4232733780262361E-02,0.4232733780262361E-02, 
    0.6812146781301679E-02,0.6812146781301679E-02, 
    0.4593608235278822E-01,0.4593608235278822E-01, 
    0.5792082516805559E-01,0.5792082516805559E-01, 
    0.1593195877796943E-01,0.1593195877796943E-01, 
    0.3056610289854859E-01,0.3056610289854859E-01, 
    0.8554835231990672E-02,0.8554835231990672E-02, 
    0.1956601012663267E-01,0.1956601012663267E-01, 
    0.1051117872480982E-01,0.1051117872480982E-01, 
    0.6521513110348300E-02,0.6521513110348300E-02, 
    0.5717074378851367E-01,0.5717074378851367E-01, 
    0.1929092010649333E-01,0.1929092010649333E-01, 
    0.5062498553433200E-01,0.5062498553433200E-01, 
    0.5402053691556560E-01,0.5402053691556560E-01, 
    0.8931127348310931E-02,0.8931127348310931E-02, 
    0.1271750255435475E-01,0.1271750255435475E-01, 
    0.3860453155876481E-01,0.3860453155876481E-01, 
    0.2372095848832025E-01,0.2372095848832025E-01, 
    0.5971489271210207E-02,0.5971489271210207E-02, 
    0.3192818361403054E-01,0.3192818361403054E-01, 
    0.4118255231945603E-02,0.4118255231945603E-02, 
    0.1843978737274201E-01,0.5393241195495463E-01, 
    0.5393241195495463E-01,0.2787139422709562E-01, 
    0.2787139422709562E-01,0.8002087241718343E-02, 
    0.8002087241718343E-02,0.5159834028120074E-02, 
    0.5159834028120074E-02,0.3908234388553434E-01, 
    0.3908234388553434E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule26 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE26 returns the rule of degree 26.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.4833333960601769E+00,0.8342454172040494E+00, 
    0.7576455307285274E+00,-.5897785388721970E+00, 
    -.1965092140379674E+00,-.7806410759767570E+00, 
    0.7735458124538029E+00,0.2194354253840281E+00, 
    -.6123258347540325E+00,-.3875195310855452E+00, 
    0.4210764575469916E+00,0.6717143218134665E+00, 
    0.4937567886903959E+00,0.4310467214948753E+00, 
    -.8221934052978648E-01,0.1474744561277340E+00, 
    -.7849971862231101E+00,0.2788743412255506E+00, 
    -.9174091785830387E+00,-.2680062461938560E+00, 
    -.6869484833848130E+00,-.5769800391396182E+00, 
    -.3175999245771031E+00,-.6900722826105213E+00, 
    -.4494495234589641E+00,0.2720327015362536E+00, 
    0.9844767989506050E+00,0.9539566427635388E+00, 
    -.6618689211018638E+00,0.9948047635608559E+00, 
    0.6791817577488417E+00,-.7660679944748241E+00, 
    -.3544317327177745E+00,0.8846710719994197E+00, 
    -.4769291205740912E+00,0.3374006730427574E+00, 
    -.5230046201031646E+00,0.9773104137186888E+00, 
    0.8460741659581372E+00,0.7811230642538761E+00, 
    -.7917954564682173E+00,0.6689548567706126E+00, 
    0.8995533103742481E+00,0.5803086536096742E+00, 
    0.9125294970192420E+00,0.9864324103956210E+00, 
    -.9634751680719612E+00,0.9171821405508755E+00, 
    -.8464504414046380E+00,0.1148682640972169E-01, 
    -.8893878433260284E+00,0.1985037063041746E+00, 
    0.6321992170210116E+00,-.8764667800594454E-01, 
    0.9612042270501360E+00,0.8449547678280973E+00, 
    -.3143374729428800E-01,-.7221561416927453E+00, 
    -.7408359655977623E+00,-.6573270993518512E+00, 
    0.7579600253683801E+00,0.1336947783582028E+00, 
    -.8801486144575378E+00,0.9863146730096674E+00, 
    0.5018023174101623E+00,0.5456331728950238E+00, 
    0.9966649790149599E+00,-.3664122751853194E-01, 
    0.3253875792212299E-01,0.9510103999741907E+00, 
    0.4259089402020555E+00,0.5239001877754490E+00, 
    0.6028877215019741E+00,0.3449303649336414E+00, 
    -.5679703819864291E+00,-.5650724066199918E+00, 
    -.8419709514596180E-01,0.8402671970302827E-01, 
    -.9108885214110795E+00,0.8192687649356658E+00, 
    -.9543907466564233E+00,-.7146063535728425E+00, 
    -.3837061453087470E+00,0.4746594317998059E-01, 
    -.1741432512718204E+00,0.2800861092084038E+00, 
    -.8207343144238776E+00,-.9935648502810311E+00, 
    -.8697847685325190E+00,-.2948353200874713E+00, 
    -.8584906175082293E+00,0.8306678364449214E+00, 
    -.9930049118920397E+00,-.4621304609915161E+00, 
    0.1460010936826870E+00,-.9839106590889559E+00, 
    -.4906721644588728E+00,-.9896275589251841E+00, 
    -.9930899945869894E+00,0.9038157695421652E-01, 
    0.9976519966165041E+00,0.9439357544987950E+00, 
    -.8915811088898837E+00,0.9306112013193892E+00, 
    -.9907255740995738E+00,-.9338876300481137E+00, 
    0.3077773696572407E+00,-.9878642619619029E+00, 
    0.6952834621750591E+00,0.8445020176382917E+00, 
    0.5386138090792357E+00,-.9518141860165342E+00, 
    -.1344815378332747E+00,-.3695179976651696E+00, 
    0.9873492554611523E+00,0.7062788639769675E+00, 
    -.6456167225977583E+00,-.3396836721033217E+00, 
    -.9362361756962438E+00,-.9626677628820013E+00, 
    0.3015160924588678E+00,0.7067532860197516E+00, 
    -.9571646423196823E+00,0.9839100297235785E+00, 
    -.2419685587001514E+00,0.9438884112343363E+00, 
    -.1742260222577426E+00 };
  static double ys[] = {
    -.4215090372415874E-01,-.3692934551145019E+00, 
    -.1931519618055097E+00,0.3988972577148074E+00, 
    -.1715260024620762E+00,0.9094199619403520E+00, 
    0.6714959797593258E+00,0.3368323178193660E+00, 
    0.6925643585799712E+00,-.3334911148152457E+00, 
    0.1579601812062493E+00,0.7212403727697902E+00, 
    0.7924318088830316E+00,0.4711520853394768E+00, 
    0.2329014884732836E+00,0.7688776471476298E-01, 
    -.2427208512979698E+00,-.9819844652166957E+00, 
    0.6275746962958220E+00,0.3958917640748951E+00, 
    0.8330411855203188E+00,-.4968412520921159E+00, 
    0.4887729302167374E-01,-.3604995467485130E+00, 
    0.5502876418517838E+00,-.1552754382216172E+00, 
    0.5038985691395298E+00,-.3255794611786079E+00, 
    0.3310081213549004E+00,-.6602065173414950E+00, 
    -.4894888916356721E+00,0.5399510111685960E+00, 
    0.7655039704316382E+00,-.5820392397031763E+00, 
    0.1995474060499386E+00,0.8781656258066206E+00, 
    -.1525213978780895E+00,-.5001983735790505E+00, 
    0.7986695769024170E+00,0.4824228603329770E+00, 
    0.1592378639694050E+00,0.1732234025449527E-02, 
    -.9985738973637792E-01,0.5955792594614818E+00, 
    0.3498450621988221E+00,0.2407129566690367E+00, 
    0.9947445171943684E+00,0.6288090977170689E+00, 
    0.7414927225667501E+00,-.9929600347714176E+00, 
    0.3684373427277031E+00,-.9535682594324043E+00, 
    0.3034752604286319E+00,0.8507770288950726E+00, 
    0.6305231708999399E-01,-.8678090617845265E+00, 
    -.4141884020649694E-01,0.9931812444078701E+00, 
    -.6492023240589931E+00,-.6795234046038703E-02, 
    -.7116015253442534E+00,-.5912100888422367E+00, 
    0.9669064714907090E+00,-.9896978591786247E+00, 
    -.6357772860044338E+00,-.2816935374305686E+00, 
    0.9465152907473502E+00,0.9835467134448914E+00, 
    0.5163644022959618E+00,0.9899894348223212E+00, 
    -.9119806917142169E+00,0.9594778630442550E+00, 
    -.8263115619511332E+00,-.4386676054811558E+00, 
    0.9540273279949707E+00,-.7748522899685367E+00, 
    -.7271449178503047E+00,-.8683831780306135E+00, 
    -.9979822139283684E+00,0.1525566715594689E+00, 
    -.3052870542039953E+00,-.8823794421118246E+00, 
    -.6417482148168809E+00,-.3240255976105304E+00, 
    0.6580511661237559E+00,0.6610461063468759E+00, 
    -.9615342517942168E+00,-.8967208208364033E-01, 
    -.4912350350325419E+00,-.8431987564696162E+00, 
    -.7917868223215425E+00,-.9934797396330936E+00, 
    0.4007666180051725E+00,0.8621553422717932E+00, 
    0.9342344412005841E+00,-.9702959974411298E+00, 
    -.9316660034779380E+00,0.7509009750438762E+00, 
    -.5058825420350775E+00,0.7766280231903698E+00, 
    -.1591105418860929E+00,-.7614693278871865E+00, 
    -.6171441191167790E-01,-.9548062395608765E+00, 
    -.8211827496355327E+00,-.9063342142005131E+00, 
    0.9909185703328497E+00,0.9397084665195838E+00, 
    0.8902356582373753E+00,0.9582397465258308E+00, 
    -.9888160497196922E+00,-.6701339599706453E+00, 
    -.9449104256261275E+00,0.9907698250163826E+00, 
    -.8840101170670746E+00,0.9954576420924693E+00, 
    -.9868058887952446E+00,-.9884037853120492E+00, 
    0.8605082508169065E+00,0.1586026433797343E+00, 
    -.7662425769553762E+00,-.9472631440912187E+00, 
    0.5782109944873408E+00,0.7561475499065624E+00, 
    0.9307512366430157E+00,0.8850736835159266E+00, 
    -.4904564371890563E+00 };
  static double ws[] = {
    0.2913207807844331E-01,0.2275827975953009E-01, 
    0.2791249243709578E-01,0.1938452906022129E-01, 
    0.4364917601570797E-01,0.1237342034594484E-01, 
    0.1138865003097454E-01,0.4923704893009620E-01, 
    0.3063932561386288E-01,0.4428512850707964E-01, 
    0.4629109684697653E-01,0.1691765251942423E-01, 
    0.2741532810734821E-01,0.4052530459111174E-01, 
    0.5261109228291892E-01,0.4417466900490126E-01, 
    0.3058410753870397E-01,0.6406355493167926E-02, 
    0.6769861004210370E-02,0.4884106698203454E-01, 
    0.1697345666295625E-01,0.3712422414211005E-01, 
    0.4783538888580537E-01,0.2351215352613192E-01, 
    0.3967485550125242E-01,0.5227057590404430E-01, 
    0.7706038595693665E-02,0.1393888240267354E-01, 
    0.2605232454276587E-01,0.3221602953245650E-02, 
    0.3639783162612337E-01,0.3177013439355753E-01, 
    0.2662783889362856E-01,0.2133361104980210E-01, 
    0.4324098371978138E-01,0.2186171413802145E-01, 
    0.3997793758416952E-01,0.7491762737168889E-02, 
    0.1718989677527449E-01,0.3150610588787735E-01, 
    0.3338285362580999E-01,0.3522688180840185E-01, 
    0.2250843481977659E-01,0.2767390941687784E-01, 
    0.2163076701762542E-01,0.6370202612529778E-02, 
    0.1292816837261490E-02,0.1772036481071653E-01, 
    0.2045516332866545E-01,0.5029942715423091E-02, 
    0.2555860092259501E-01,0.8908434653957875E-02, 
    0.4462964590144296E-01,0.2350361388650933E-01, 
    0.1300162042596963E-01,0.1619762342723978E-01, 
    0.3633406947058673E-01,0.4148575443002976E-02, 
    0.3043458935989431E-01,0.3780056077982088E-01, 
    0.2670761755118099E-01,0.4679248624556967E-01, 
    0.6655017468306593E-02,0.1016870301128343E-02, 
    0.3956208492384800E-01,0.4422000056024995E-01, 
    0.1600493333371462E-02,0.9542615574085509E-02, 
    0.5101463103747718E-01,0.2316486911009150E-02, 
    0.2082363944511977E-01,0.1450858408533323E-01, 
    0.2593506316945297E-01,0.4854771494052668E-01, 
    0.1406093874479688E-01,0.2980564799727076E-01, 
    0.4032704717924811E-01,0.2836242586524349E-01, 
    0.1662837317311870E-02,0.3331747639887452E-01, 
    0.1737244435785181E-01,0.1818621007558914E-01, 
    0.4099371056923931E-01,0.5481494042070843E-01, 
    0.4070179363666812E-01,0.4218182579123851E-01, 
    0.8533899793279886E-02,0.5949667491635394E-02, 
    0.2693305692079543E-01,0.3035811941473991E-01, 
    0.1815553496055310E-01,0.3365556943518273E-02, 
    0.5641437853076681E-02,0.1997142107835644E-01, 
    0.1803786451914521E-01,0.2269810849285899E-02, 
    0.1824056247133745E-01,0.4961745614234228E-02, 
    0.5229975396089418E-02,0.3090226019029179E-01, 
    0.4530557980884387E-02,0.1260165841456744E-01, 
    0.2688472868433625E-01,0.6305207130138189E-02, 
    0.4114195161424335E-02,0.8573396052945570E-02, 
    0.6467557193927746E-02,0.2703578258059107E-02, 
    0.1975997243239372E-01,0.9381005409995410E-02, 
    0.6472606174579578E-02,0.1382406712148664E-01, 
    0.1805324127586059E-01,0.6329200107849942E-02, 
    0.3895024567992489E-02,0.3809507629042853E-02, 
    0.6517542895566797E-02,0.7169776873936058E-02, 
    0.1060178841340027E-01,0.1641915380777178E-01, 
    0.3678722762764779E-01,0.1379818565014169E-01, 
    0.1123917596421086E-01,0.6410312514956618E-02, 
    0.1884563206279892E-01,0.8912626526178684E-02, 
    0.4962602708111043E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule27 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE27 returns the rule of degree 27.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.3179234344617641E+00,0.9945426115270097E+00, 
    0.3179234344617642E+00,-.9945426115270097E+00, 
    -.9658581524130557E+00,0.9842821312283567E+00, 
    0.9658581524130559E+00,-.9842821312283565E+00, 
    0.4546608217543911E+00,0.9794570195190775E+00, 
    -.4546608217543910E+00,-.9794570195190775E+00, 
    0.9203445363225692E+00,0.9966680243748195E+00, 
    -.9203445363225690E+00,-.9966680243748197E+00, 
    0.5921826557247571E+00,0.9227007102569557E+00, 
    -.5921826557247569E+00,-.9227007102569557E+00, 
    -.6265556238656111E+00,0.9943315661391944E+00, 
    0.6265556238656114E+00,-.9943315661391944E+00, 
    -.9106466233445932E+00,0.9127106311205702E+00, 
    0.9106466233445935E+00,-.9127106311205699E+00, 
    -.4796281374834562E+00,0.9651530210115425E+00, 
    0.4796281374834563E+00,-.9651530210115425E+00, 
    0.6298732844502039E+00,0.9921721529860120E+00, 
    -.6298732844502036E+00,-.9921721529860120E+00, 
    -.8335744040052592E+00,0.9824660135148449E+00, 
    0.8335744040052594E+00,-.9824660135148449E+00, 
    0.1959395767367455E+00,0.9899174626675955E+00, 
    -.1959395767367454E+00,-.9899174626675955E+00, 
    -.1386302484544642E+00,0.9830581303940688E+00, 
    0.1386302484544643E+00,-.9830581303940688E+00, 
    -.7022171903082602E+00,0.9292585099108952E+00, 
    0.7022171903082605E+00,-.9292585099108952E+00, 
    0.7970439287784883E+00,0.9675649927064782E+00, 
    -.7970439287784881E+00,-.9675649927064782E+00, 
    -.8020178388613470E+00,0.8160591260895965E+00, 
    0.8020178388613470E+00,-.8160591260895965E+00, 
    -.5398233377470387E+00,0.8364692379579785E+00, 
    0.5398233377470387E+00,-.8364692379579785E+00, 
    -.1520649596611693E+00,0.8369974777769484E+00, 
    0.1520649596611695E+00,-.8369974777769484E+00, 
    -.2998875490507344E+00,0.3409328206751328E+00, 
    0.2998875490507344E+00,-.3409328206751328E+00, 
    -.3209420000531037E+00,0.9175861963438381E+00, 
    0.3209420000531039E+00,-.9175861963438381E+00, 
    -.9020677393268002E-01,0.1670667038911821E+00, 
    0.9020677393268005E-01,-.1670667038911821E+00, 
    0.1926515070565072E+00,0.5293872532788142E+00, 
    -.1926515070565071E+00,-.5293872532788142E+00, 
    0.4976080851052046E+00,0.7953917795116896E+00, 
    -.4976080851052045E+00,-.7953917795116896E+00, 
    0.3002782637358375E-01,0.9362295643614138E+00, 
    -.3002782637358363E-01,-.9362295643614138E+00, 
    0.6881617137006012E+00,0.8842713205559528E+00, 
    -.6881617137006012E+00,-.8842713205559528E+00, 
    0.2200611683790181E+00,0.8481900686404982E+00, 
    -.2200611683790180E+00,-.8481900686404982E+00, 
    0.4531153802420229E-01,0.3769373245611116E+00, 
    -.4531153802420225E-01,-.3769373245611116E+00, 
    -.4941549461445193E+00,0.5167583593472795E+00, 
    0.4941549461445194E+00,-.5167583593472795E+00, 
    0.3405123948054840E+00,0.6649111597000639E+00, 
    -.3405123948054839E+00,-.6649111597000639E+00, 
    -.1521097283009047E+00,0.5602567958473701E+00, 
    0.1521097283009048E+00,-.5602567958473701E+00, 
    0.3710843473830662E+00,0.9336952377043988E+00, 
    -.3710843473830661E+00,-.9336952377043988E+00, 
    -.3536770731498163E+00,0.7123617422035845E+00, 
    0.3536770731498164E+00,-.7123617422035845E+00, 
    -.6668117984035637E+00,0.6750511711635266E+00, 
    0.6668117984035637E+00,-.6750511711635266E+00, 
    0.3880331998977081E-01,0.7201416766661869E+00, 
    -.3880331998977073E-01,-.7201416766661869E+00 };
  static double ys[] = {
    -.9945426115270097E+00,-.3179234344617641E+00, 
    0.9945426115270097E+00,0.3179234344617642E+00, 
    -.9842821312283568E+00,-.9658581524130558E+00, 
    0.9842821312283566E+00,0.9658581524130561E+00, 
    -.9794570195190775E+00,0.4546608217543910E+00, 
    0.9794570195190775E+00,-.4546608217543909E+00, 
    -.9966680243748194E+00,0.9203445363225691E+00, 
    0.9966680243748196E+00,-.9203445363225689E+00, 
    -.9227007102569557E+00,0.5921826557247570E+00, 
    0.9227007102569557E+00,-.5921826557247568E+00, 
    -.9943315661391944E+00,-.6265556238656113E+00, 
    0.9943315661391944E+00,0.6265556238656115E+00, 
    -.9127106311205703E+00,-.9106466233445933E+00, 
    0.9127106311205700E+00,0.9106466233445936E+00, 
    -.9651530210115425E+00,-.4796281374834563E+00, 
    0.9651530210115425E+00,0.4796281374834564E+00, 
    -.9921721529860120E+00,0.6298732844502037E+00, 
    0.9921721529860120E+00,-.6298732844502035E+00, 
    -.9824660135148449E+00,-.8335744040052593E+00, 
    0.9824660135148449E+00,0.8335744040052595E+00, 
    -.9899174626675955E+00,0.1959395767367454E+00, 
    0.9899174626675955E+00,-.1959395767367453E+00, 
    -.9830581303940688E+00,-.1386302484544642E+00, 
    0.9830581303940688E+00,0.1386302484544643E+00, 
    -.9292585099108952E+00,-.7022171903082604E+00, 
    0.9292585099108952E+00,0.7022171903082606E+00, 
    -.9675649927064782E+00,0.7970439287784882E+00, 
    0.9675649927064782E+00,-.7970439287784880E+00, 
    -.8160591260895965E+00,-.8020178388613470E+00, 
    0.8160591260895965E+00,0.8020178388613470E+00, 
    -.8364692379579785E+00,-.5398233377470387E+00, 
    0.8364692379579785E+00,0.5398233377470387E+00, 
    -.8369974777769484E+00,-.1520649596611694E+00, 
    0.8369974777769484E+00,0.1520649596611695E+00, 
    -.3409328206751328E+00,-.2998875490507344E+00, 
    0.3409328206751328E+00,0.2998875490507344E+00, 
    -.9175861963438381E+00,-.3209420000531038E+00, 
    0.9175861963438381E+00,0.3209420000531039E+00, 
    -.1670667038911821E+00,-.9020677393268003E-01, 
    0.1670667038911821E+00,0.9020677393268006E-01, 
    -.5293872532788142E+00,0.1926515070565072E+00, 
    0.5293872532788142E+00,-.1926515070565071E+00, 
    -.7953917795116896E+00,0.4976080851052045E+00, 
    0.7953917795116896E+00,-.4976080851052044E+00, 
    -.9362295643614138E+00,0.3002782637358369E-01, 
    0.9362295643614138E+00,-.3002782637358357E-01, 
    -.8842713205559528E+00,0.6881617137006012E+00, 
    0.8842713205559528E+00,-.6881617137006012E+00, 
    -.8481900686404982E+00,0.2200611683790180E+00, 
    0.8481900686404982E+00,-.2200611683790179E+00, 
    -.3769373245611116E+00,0.4531153802420227E-01, 
    0.3769373245611116E+00,-.4531153802420223E-01, 
    -.5167583593472795E+00,-.4941549461445193E+00, 
    0.5167583593472795E+00,0.4941549461445194E+00, 
    -.6649111597000639E+00,0.3405123948054840E+00, 
    0.6649111597000639E+00,-.3405123948054838E+00, 
    -.5602567958473701E+00,-.1521097283009047E+00, 
    0.5602567958473701E+00,0.1521097283009048E+00, 
    -.9336952377043988E+00,0.3710843473830662E+00, 
    0.9336952377043988E+00,-.3710843473830661E+00, 
    -.7123617422035845E+00,-.3536770731498164E+00, 
    0.7123617422035845E+00,0.3536770731498165E+00, 
    -.6750511711635266E+00,-.6668117984035637E+00, 
    0.6750511711635266E+00,0.6668117984035637E+00, 
    -.7201416766661869E+00,0.3880331998977077E-01, 
    0.7201416766661869E+00,-.3880331998977068E-01 };
  static double ws[] = {
    0.2672899331427035E-02,0.2672899331427035E-02, 
    0.2672899331427035E-02,0.2672899331427035E-02, 
    0.2445906596530064E-02,0.2445906596530064E-02, 
    0.2445906596530064E-02,0.2445906596530064E-02, 
    0.5784555125932030E-02,0.5784555125932030E-02, 
    0.5784555125932030E-02,0.5784555125932030E-02, 
    0.1774682379035372E-02,0.1774682379035372E-02, 
    0.1774682379035372E-02,0.1774682379035372E-02, 
    0.9666652283416993E-02,0.9666652283416993E-02, 
    0.9666652283416993E-02,0.9666652283416993E-02, 
    0.3475754833559290E-02,0.3475754833559290E-02, 
    0.3475754833559290E-02,0.3475754833559290E-02, 
    0.9738340699044858E-02,0.9738340699044858E-02, 
    0.9738340699044858E-02,0.9738340699044858E-02, 
    0.9913445575624551E-02,0.9913445575624551E-02, 
    0.9913445575624551E-02,0.9913445575624551E-02, 
    0.4024564097795350E-02,0.4024564097795350E-02, 
    0.4024564097795350E-02,0.4024564097795350E-02, 
    0.5467642658052182E-02,0.5467642658052182E-02, 
    0.5467642658052182E-02,0.5467642658052182E-02, 
    0.6791307131143990E-02,0.6791307131143990E-02, 
    0.6791307131143990E-02,0.6791307131143990E-02, 
    0.7849471971542231E-02,0.7849471971542231E-02, 
    0.7849471971542231E-02,0.7849471971542231E-02, 
    0.1481562369739479E-01,0.1481562369739479E-01, 
    0.1481562369739479E-01,0.1481562369739479E-01, 
    0.8652776980870136E-02,0.8652776980870136E-02, 
    0.8652776980870136E-02,0.8652776980870136E-02, 
    0.1938931217339045E-01,0.1938931217339045E-01, 
    0.1938931217339045E-01,0.1938931217339045E-01, 
    0.2728047644809397E-01,0.2728047644809397E-01, 
    0.2728047644809397E-01,0.2728047644809397E-01, 
    0.2833042504936177E-01,0.2833042504936177E-01, 
    0.2833042504936177E-01,0.2833042504936177E-01, 
    0.5141433190365619E-01,0.5141433190365619E-01, 
    0.5141433190365619E-01,0.5141433190365619E-01, 
    0.1775207295130232E-01,0.1775207295130232E-01, 
    0.1775207295130232E-01,0.1775207295130232E-01, 
    0.4944421261124536E-01,0.4944421261124536E-01, 
    0.4944421261124536E-01,0.4944421261124536E-01, 
    0.3632695571196309E-01,0.3632695571196309E-01, 
    0.3632695571196309E-01,0.3632695571196309E-01, 
    0.2905709329104403E-01,0.2905709329104403E-01, 
    0.2905709329104403E-01,0.2905709329104403E-01, 
    0.1861253661440120E-01,0.1861253661440120E-01, 
    0.1861253661440120E-01,0.1861253661440120E-01, 
    0.1362514538944166E-01,0.1362514538944166E-01, 
    0.1362514538944166E-01,0.1362514538944166E-01, 
    0.2778931208298291E-01,0.2778931208298291E-01, 
    0.2778931208298291E-01,0.2778931208298291E-01, 
    0.4496161760244370E-01,0.4496161760244370E-01, 
    0.4496161760244370E-01,0.4496161760244370E-01, 
    0.4391754850935400E-01,0.4391754850935400E-01, 
    0.4391754850935400E-01,0.4391754850935400E-01, 
    0.3461629438577501E-01,0.3461629438577501E-01, 
    0.3461629438577501E-01,0.3461629438577501E-01, 
    0.4733898929183223E-01,0.4733898929183223E-01, 
    0.4733898929183223E-01,0.4733898929183223E-01, 
    0.1448889464186536E-01,0.1448889464186536E-01, 
    0.1448889464186536E-01,0.1448889464186536E-01, 
    0.3863763742526743E-01,0.3863763742526743E-01, 
    0.3863763742526743E-01,0.3863763742526743E-01, 
    0.3295515782263461E-01,0.3295515782263461E-01, 
    0.3295515782263461E-01,0.3295515782263461E-01, 
    0.3809514391912308E-01,0.3809514391912308E-01, 
    0.3809514391912308E-01,0.3809514391912308E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule28 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE28 returns the rule of degree 28.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    0.2827621856762107E+00,-.8043038962721508E+00, 
    -.3395780944817125E+00,0.9664085101515743E+00, 
    -.1096329656209190E-01,0.8193786706976756E+00, 
    -.9286317198271657E+00,0.4198022792026188E+00, 
    -.8604918125918429E+00,0.8109783092388664E+00, 
    0.5614606479838853E+00,0.1250139096754087E-01, 
    0.6890874578865651E+00,-.8003785261499256E+00, 
    -.7186772992712774E+00,0.2097379198022608E-01, 
    -.9904175004752394E+00,0.5247683324321177E+00, 
    0.9879057318966782E+00,0.9872995514825285E+00, 
    -.2010899760026548E+00,-.7150754139309137E+00, 
    -.9744175518522036E+00,-.2743957404149748E+00, 
    -.8801155850915744E-01,-.9120488211036881E+00, 
    0.5997911603089924E-01,-.9822903323396509E+00, 
    0.9933058426309856E+00,-.1094425840746250E+00, 
    0.8016032845806299E+00,-.3785896040093871E+00, 
    0.9908494000362051E+00,-.9920099650417782E+00, 
    -.5334142633484058E+00,-.6081971936243052E+00, 
    -.9913559247973139E+00,-.4420258077139650E+00, 
    0.9904922021473850E+00,-.2701059581013934E+00, 
    -.6762870118531359E+00,0.4550203878192706E+00, 
    -.9933274734742550E+00,0.1025031015215653E+00, 
    -.9908021888968861E+00,0.4018136560549200E-02, 
    0.1830465948969627E+00,-.9190332253043175E+00, 
    -.7967293798867977E+00,0.8910586530553062E+00, 
    0.9738677077447094E+00,0.4922211332598602E+00, 
    -.9935974988322385E+00,-.7469015595759197E+00, 
    -.7474851457426346E+00,-.8461768034166064E+00, 
    0.7263311602661456E+00,0.1838655693690960E+00, 
    -.8924925604346934E+00,0.3473229569871696E+00, 
    0.7888084262939273E+00,-.1622547625661196E+00, 
    0.1708950639909230E+00,-.9858834062192522E+00, 
    0.8595474037303855E+00,-.2487721441258377E+00, 
    0.9897369384324803E+00,-.8592565754231390E+00, 
    0.3552890037180957E+00,0.9343486754942414E+00, 
    -.5982760380368914E+00,0.9012497934650560E+00, 
    -.8867898943082418E+00,-.2973823530897202E+00, 
    0.2929780512625694E+00,0.5109886162010564E+00, 
    -.6362138293509048E+00,0.9840555620580882E+00, 
    0.6242255274464619E+00,0.9500492760473710E+00, 
    0.6206677950945041E+00,0.3667200892151265E+00, 
    -.3354257519635134E+00,0.6532050712129780E+00, 
    0.8340712667685216E+00,0.9517120425530485E+00, 
    0.4939764190001616E+00,0.8519117390868273E+00, 
    -.4641713197529147E+00,0.3096864797777032E+00, 
    -.9951406550164114E+00,0.6849334791605266E+00, 
    0.5377487266764005E+00,0.6596506837798173E+00, 
    -.6190982139272722E+00,0.9676681958645598E+00, 
    -.8663734426507712E+00,0.9229408505602535E+00, 
    0.8925298503917507E+00,0.1376340175855351E+00, 
    -.4255923534036878E+00,-.5069417238664446E+00, 
    -.9581264401591332E+00,0.7956423904664977E+00, 
    -.9438198089736733E+00,-.9602193183322451E+00, 
    -.9546324794056378E+00,0.7865853669573317E+00, 
    0.7640178616435913E+00,-.9037187541181603E-01, 
    0.6583886775017005E+00,0.9394184495408773E+00, 
    0.1210974568433837E+00,0.3181846883529323E+00, 
    -.9476649688000099E+00,-.7447346202738520E+00, 
    0.7669083321105266E+00,0.3177180418892223E+00, 
    0.9005374034994568E+00,0.4905235558364557E+00, 
    -.2799102574207112E+00,-.7756668408683119E+00, 
    -.7847585523489610E+00,-.5947448678197160E+00, 
    -.6835954896209653E-01,-.4603444798372620E+00, 
    0.1267789889405305E+00,0.9957007235251634E+00, 
    -.9555617449567378E+00,0.4926175951876377E+00, 
    0.9580744421864001E+00,0.6544409855907934E+00, 
    -.2589485669037206E+00,-.4625383187120728E+00, 
    0.8779137327105354E+00,-.8753270142023138E+00, 
    -.7680237343653598E+00,-.4646141113930034E+00, 
    -.6237643000374615E+00,-.1004225708457810E+00, 
    0.2934408497776114E+00,-.8838031019815079E+00, 
    -.8252198025376763E-01,0.1013709068378068E+00, 
    -.6355163316600663E+00,-.2814473462899598E+00, 
    0.9991999747547105E+00 };
  static double ys[] = {
    -.2146915742023843E+00,-.3602259227928870E+00, 
    -.2366070789246893E+00,-.1067737835329393E+00, 
    -.7832467667692294E-01,0.5090533388398083E-02, 
    0.5655320173536176E-02,-.9091460408400718E+00, 
    -.3610436689009546E+00,0.1173540113566495E+00, 
    -.7580822207698504E+00,0.7636857906206649E+00, 
    -.6742700466367159E+00,0.3605821467647248E+00, 
    0.2565264593758684E+00,0.8682397249203467E+00, 
    -.9887795883901583E+00,-.8777099834380453E+00, 
    -.9937283708213308E+00,0.9906296802548810E+00, 
    0.8788092453036578E+00,0.9706849956765652E+00, 
    -.7829696251737753E-01,0.9937649402523892E+00, 
    -.1524877190213388E+00,-.2242135513466236E+00, 
    0.9919534212785273E+00,0.9960588898544035E+00, 
    -.5034894383669408E+00,0.9599951355029216E+00, 
    -.5645937597081716E+00,0.9447996831444674E+00, 
    0.2711989824047218E+00,0.1447359435070455E+00, 
    0.9871535260388329E+00,0.9140753757822478E+00, 
    -.2960350108035981E+00,0.8371402678277077E+00, 
    -.7650309750102933E+00,-.3188681363434392E+00, 
    -.9938368658433656E+00,-.6339545710762906E+00, 
    -.6212095493273496E+00,-.9934382054756534E+00, 
    0.4715389357436094E+00,-.7215460671816014E+00, 
    0.6773640866690520E+00,0.9779955899159012E+00, 
    0.9970513707078691E+00,-.1551145670991118E+00, 
    -.3024178806712221E+00,-.9926567100374114E+00, 
    -.8660993570031557E+00,-.4956056862948234E+00, 
    0.8382998942855415E+00,0.9259687805728221E+00, 
    -.8351882927550272E+00,0.9386358290968080E+00, 
    -.9900384261478909E+00,-.4464745325070345E+00, 
    -.9936162118016632E+00,-.8235872679034352E+00, 
    -.5983594352543780E+00,0.9398154839397104E+00, 
    -.7371773517980565E+00,0.7366500643247871E+00, 
    0.8982625204679442E+00,0.1211489804261566E+00, 
    0.9859353129167986E+00,-.8465334600907201E+00, 
    0.7342041846853576E+00,-.4349009915613299E+00, 
    0.5008348384608579E+00,-.9942925485328575E+00, 
    -.9626593056228334E+00,-.2709700923096490E+00, 
    0.5944472663027801E-01,-.9248481403512147E+00, 
    -.4757878778805171E+00,0.9921234070894396E-01, 
    0.9864359054124004E+00,0.5564466223050852E+00, 
    -.9075561035062022E+00,0.6130603473347437E+00, 
    -.9196402734279348E+00,-.6285797630023995E+00, 
    0.7360652264674404E+00,0.9946281004658378E+00, 
    -.1056502846403645E+00,0.8343858336015567E+00, 
    0.7795434558480316E+00,0.2437361149690468E+00, 
    0.4083791315159580E+00,-.8237952273107227E-01, 
    -.6450626080823243E+00,0.4764961946576270E+00, 
    0.7429354103100736E+00,-.9752458082257867E+00, 
    0.3137605530793695E+00,0.7602258278505522E-01, 
    0.6028197206610887E+00,-.9665994542595457E+00, 
    -.9448687562971092E+00,0.7450892979512898E+00, 
    0.8580812769827701E+00,0.6469685322672528E+00, 
    -.4730667721016463E+00,0.4691067254414554E+00, 
    0.9468223400568799E+00,0.5747271977509141E+00, 
    0.8487224420002684E+00,0.9604420521431418E+00, 
    0.4325756220233750E+00,-.1154193690840114E+00, 
    0.2996859346038543E+00,0.6079262704780987E+00, 
    -.3014753479396839E+00,0.2694764821728903E+00, 
    0.6210587156046703E+00,0.9250623306784818E+00, 
    0.4091950530277537E+00,-.1149400908487048E+00, 
    -.9509596164554358E+00,0.4489210023756193E+00, 
    0.2491367754635779E+00,0.2328542094943754E+00, 
    -.3147181280829211E+00,0.6575404381717325E+00, 
    -.7572545223612748E+00,0.8744761196977546E-01, 
    0.7809511457782177E+00,-.9573612252124380E+00, 
    0.5984333435621010E-01,-.4795365437914892E+00, 
    0.8758201084343584E+00,-.6288667697841542E+00, 
    -.7704960778314179E+00,-.7807547746986493E+00, 
    -.2948221971565224E+00,-.9595773298541288E+00, 
    -.7859464731741874E+00,-.8701241408788146E+00, 
    -.4848613986697829E+00,-.8903571922698699E+00, 
    -.8832911479917025E+00,-.6435327316304572E+00, 
    -.5075800418230667E-01 };
  static double ws[] = {
    0.1740144063407695E-02,0.6784563755290567E-02, 
    0.1209466540188899E-01,0.5276845120506419E-02, 
    0.2120637597472411E-01,0.1500950169216209E-01, 
    0.1016605179698056E-01,0.1052237804962158E-01, 
    0.1389531092331450E-01,0.1799231574573814E-01, 
    0.1681749672744976E-01,0.2190516172058484E-01, 
    0.1802050074136916E-01,0.1929827700816552E-01, 
    0.2322655354121691E-01,0.1715984331158577E-01, 
    0.7585152723116851E-03,0.1423586238543194E-01, 
    0.6209204642340611E-03,0.8055582944619972E-03, 
    0.1686431398199053E-01,0.6130199135277467E-02, 
    0.8266675752855135E-02,0.3947441053087870E-02, 
    0.3763432904885534E-01,0.1544537020393182E-01, 
    0.4895887118324937E-02,0.5766157576611153E-03, 
    0.3911408899973994E-02,0.1134055149899166E-01, 
    0.2045540821981484E-01,0.1254460489172730E-01, 
    0.5370936278809426E-02,0.5162980339438905E-02, 
    0.5510957955189620E-02,0.1348681146014646E-01, 
    0.5293461261550956E-02,0.2135158963906704E-01, 
    0.3872040956962797E-02,0.4041863730898066E-01, 
    0.3480075985572175E-02,0.3082881220603644E-01, 
    0.3934841673430547E-02,0.4955048643733777E-02, 
    0.5332713420324274E-02,0.3127092112918890E-01, 
    0.3263868066569824E-01,0.3643906126193626E-02, 
    0.1773410002914162E-02,0.2031653637703585E-01, 
    0.9778965531187402E-02,0.4720821658621592E-02, 
    0.2492177802034716E-02,0.2644508244904856E-01, 
    0.1656291578773864E-01,0.9176500057842589E-02, 
    0.1732347887472145E-01,0.1561854636309510E-01, 
    0.2892151505244423E-02,0.3887473612274506E-01, 
    0.3107151229520027E-02,0.2613845455837795E-01, 
    0.3698667210865119E-01,0.2566493481700192E-02, 
    0.1615348439627276E-01,0.3075142974941159E-01, 
    0.2904460597417912E-02,0.2398826862979595E-01, 
    0.7145392996633393E-02,0.8978366822810722E-02, 
    0.2612889620260116E-01,0.1874378199539706E-01, 
    0.1917559308928981E-01,0.4661196658021693E-02, 
    0.1235879177909568E-01,0.4010025698621002E-01, 
    0.3740253247098928E-01,0.3148001303878885E-02, 
    0.3330951974159902E-01,0.1500177718625442E-01, 
    0.6070126982071918E-02,0.3772728641321498E-01, 
    0.1934239152098538E-01,0.2980908782028305E-01, 
    0.1060149985804505E-01,0.1169140727671868E-01, 
    0.2913863637730209E-01,0.2512926081791887E-02, 
    0.4408214653268561E-01,0.2616523021560373E-01, 
    0.2824238355893381E-02,0.3534124084616964E-01, 
    0.3855284013804392E-01,0.3769601120920466E-01, 
    0.3068556194531641E-01,0.1107117178829623E-01, 
    0.1675502525295356E-01,0.4250234458185585E-02, 
    0.2160027234544727E-01,0.4964511475611218E-01, 
    0.3661675186805007E-01,0.1113946092997709E-01, 
    0.4594039996241282E-02,0.2058836038539646E-01, 
    0.8608231530163511E-02,0.1071439010120374E-01, 
    0.1325231032613585E-01,0.2793585524817994E-01, 
    0.1060087465526644E-01,0.4186289109506627E-01, 
    0.2038937560848930E-01,0.4779376144301057E-02, 
    0.4600352201496827E-01,0.4850236286164095E-01, 
    0.1563580434067219E-01,0.2738173661329469E-01, 
    0.3152680760913545E-01,0.4699759394230683E-01, 
    0.1758154363368765E-01,0.1705163560097511E-01, 
    0.4564749915571795E-01,0.3274246294590765E-01, 
    0.9932190596257133E-02,0.3764556185875566E-01, 
    0.5058095830353233E-01,0.4537296603002153E-01, 
    0.4951708484325965E-01,0.3279405659962938E-02, 
    0.9987676851635598E-02,0.4573491590126628E-01, 
    0.9304988468918474E-02,0.1136613167624352E-01, 
    0.5119306902579393E-01,0.4198377464395852E-01, 
    0.1182329164018337E-01,0.1992962132997380E-01, 
    0.2179036392120167E-01,0.2957452226482192E-01, 
    0.3999530391333647E-01,0.1478968751484580E-01, 
    0.3170728715126719E-01,0.1199170942901081E-01, 
    0.4711944799135396E-01,0.2445849079621999E-01, 
    0.1955624669135905E-01,0.4043330460830744E-01, 
    0.2815879734180029E-02 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule29 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE29 returns the rule of degree 29.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.9957609056803307E+00,0.9554778819247284E+00, 
    0.9957609056803309E+00,-.9554778819247282E+00, 
    -.8018419403764031E+00,0.9936244517402659E+00, 
    0.8018419403764033E+00,-.9936244517402659E+00, 
    -.9579902382554616E+00,0.9905651653257478E+00, 
    0.9579902382554618E+00,-.9905651653257476E+00, 
    0.2638223412555635E-01,0.9935202992997069E+00, 
    -.2638223412555623E-01,-.9935202992997069E+00, 
    -.5066682706262957E+00,0.9913300131756428E+00, 
    0.5066682706262959E+00,-.9913300131756428E+00, 
    -.2298993811511364E+00,0.9820954997507334E+00, 
    0.2298993811511365E+00,-.9820954997507334E+00, 
    0.4377228178889339E+00,0.9940035328985276E+00, 
    -.4377228178889338E+00,-.9940035328985276E+00, 
    -.6971536278750922E+00,0.9699469973610456E+00, 
    0.6971536278750924E+00,-.9699469973610456E+00, 
    -.1287796452538595E+00,0.9374257301963881E+00, 
    0.1287796452538596E+00,-.9374257301963881E+00, 
    0.2219765871034010E+00,0.9669818379796637E+00, 
    -.2219765871034008E+00,-.9669818379796637E+00, 
    0.8793284704221047E+00,0.9627913549967865E+00, 
    -.8793284704221045E+00,-.9627913549967865E+00, 
    -.4224857795905392E+00,0.8997395988035505E+00, 
    0.4224857795905393E+00,-.8997395988035505E+00, 
    0.7561015490560310E+00,0.9929163478894405E+00, 
    -.7561015490560308E+00,-.9929163478894405E+00, 
    0.2640620995843525E+00,0.8673709698403427E+00, 
    -.2640620995843524E+00,-.8673709698403427E+00, 
    -.8876763767236916E+00,0.9499938690143425E+00, 
    0.8876763767236918E+00,-.9499938690143425E+00, 
    0.6107736243515250E+00,0.9600495691709722E+00, 
    -.6107736243515248E+00,-.9600495691709722E+00, 
    -.7868750412424903E+00,0.8719483642340994E+00, 
    0.7868750412424903E+00,-.8719483642340994E+00, 
    -.1141990095980696E-01,0.9167780797645280E+00, 
    0.1141990095980707E-01,-.9167780797645280E+00, 
    0.6200952646798407E+00,0.8398067949907553E+00, 
    -.6200952646798407E+00,-.8398067949907553E+00, 
    -.4280702467551835E+00,0.3950643020883297E+00, 
    0.4280702467551835E+00,-.3950643020883297E+00, 
    0.3401810242198736E+00,0.6973120201514252E+00, 
    -.3401810242198735E+00,-.6973120201514252E+00, 
    -.2617868665129255E+00,0.2355937519408776E+00, 
    0.2617868665129255E+00,-.2355937519408775E+00, 
    0.1411576681950311E+00,0.8161158957508043E+00, 
    -.1411576681950309E+00,-.8161158957508043E+00, 
    -.3695797105378154E-01,0.6981295844594165E+00, 
    0.3695797105378162E-01,-.6981295844594165E+00, 
    -.6312274766498698E+00,0.9353634542809267E+00, 
    0.6312274766498700E+00,-.9353634542809267E+00, 
    -.1150463641736545E+00,0.9674061994537117E-01, 
    0.1150463641736545E+00,-.9674061994537114E-01, 
    -.3747271550737231E+00,0.9438572336357420E+00, 
    0.3747271550737232E+00,-.9438572336357420E+00, 
    0.5008869009053424E+00,0.7904714944143627E+00, 
    -.5008869009053424E+00,-.7904714944143627E+00, 
    0.1608699923234913E+00,0.5513908795637696E+00, 
    -.1608699923234913E+00,-.5513908795637696E+00, 
    -.5657668007893579E+00,0.8210097847827826E+00, 
    0.5657668007893579E+00,-.8210097847827826E+00, 
    0.7712590259177344E+00,0.8970730060299638E+00, 
    -.7712590259177344E+00,-.8970730060299638E+00, 
    0.4163780011554897E+00,0.9113544664316420E+00, 
    -.4163780011554896E+00,-.9113544664316420E+00, 
    -.5836842069804584E+00,0.5579917371493338E+00, 
    0.5836842069804584E+00,-.5579917371493338E+00, 
    -.2215359986069006E+00,0.8254893061838994E+00, 
    0.2215359986069007E+00,-.8254893061838994E+00, 
    -.4009765829120037E+00,0.6957123574117627E+00, 
    0.4009765829120038E+00,-.6957123574117627E+00, 
    -.2343918950839540E-01,0.3763993540852892E+00, 
    0.2343918950839545E-01,-.3763993540852892E+00, 
    -.7105079759089011E+00,0.7179044391714098E+00, 
    0.7105079759089011E+00,-.7179044391714098E+00, 
    -.2200751569931638E+00,0.5431094815533654E+00, 
    0.2200751569931638E+00,-.5431094815533654E+00 };
  static double ys[] = {
    -.9554778819247285E+00,-.9957609056803308E+00, 
    0.9554778819247283E+00,0.9957609056803310E+00, 
    -.9936244517402659E+00,-.8018419403764032E+00, 
    0.9936244517402659E+00,0.8018419403764034E+00, 
    -.9905651653257479E+00,-.9579902382554617E+00, 
    0.9905651653257477E+00,0.9579902382554619E+00, 
    -.9935202992997069E+00,0.2638223412555629E-01, 
    0.9935202992997069E+00,-.2638223412555617E-01, 
    -.9913300131756428E+00,-.5066682706262958E+00, 
    0.9913300131756428E+00,0.5066682706262960E+00, 
    -.9820954997507334E+00,-.2298993811511365E+00, 
    0.9820954997507334E+00,0.2298993811511366E+00, 
    -.9940035328985276E+00,0.4377228178889339E+00, 
    0.9940035328985276E+00,-.4377228178889337E+00, 
    -.9699469973610456E+00,-.6971536278750923E+00, 
    0.9699469973610456E+00,0.6971536278750925E+00, 
    -.9374257301963881E+00,-.1287796452538595E+00, 
    0.9374257301963881E+00,0.1287796452538596E+00, 
    -.9669818379796637E+00,0.2219765871034009E+00, 
    0.9669818379796637E+00,-.2219765871034008E+00, 
    -.9627913549967865E+00,0.8793284704221046E+00, 
    0.9627913549967865E+00,-.8793284704221044E+00, 
    -.8997395988035505E+00,-.4224857795905393E+00, 
    0.8997395988035505E+00,0.4224857795905394E+00, 
    -.9929163478894405E+00,0.7561015490560309E+00, 
    0.9929163478894405E+00,-.7561015490560307E+00, 
    -.8673709698403427E+00,0.2640620995843524E+00, 
    0.8673709698403427E+00,-.2640620995843523E+00, 
    -.9499938690143425E+00,-.8876763767236917E+00, 
    0.9499938690143425E+00,0.8876763767236919E+00, 
    -.9600495691709722E+00,0.6107736243515249E+00, 
    0.9600495691709722E+00,-.6107736243515247E+00, 
    -.8719483642340994E+00,-.7868750412424903E+00, 
    0.8719483642340994E+00,0.7868750412424903E+00, 
    -.9167780797645280E+00,-.1141990095980701E-01, 
    0.9167780797645280E+00,0.1141990095980712E-01, 
    -.8398067949907553E+00,0.6200952646798407E+00, 
    0.8398067949907553E+00,-.6200952646798407E+00, 
    -.3950643020883297E+00,-.4280702467551835E+00, 
    0.3950643020883297E+00,0.4280702467551835E+00, 
    -.6973120201514252E+00,0.3401810242198735E+00, 
    0.6973120201514252E+00,-.3401810242198734E+00, 
    -.2355937519408776E+00,-.2617868665129255E+00, 
    0.2355937519408775E+00,0.2617868665129255E+00, 
    -.8161158957508043E+00,0.1411576681950310E+00, 
    0.8161158957508043E+00,-.1411576681950309E+00, 
    -.6981295844594165E+00,-.3695797105378158E-01, 
    0.6981295844594165E+00,0.3695797105378167E-01, 
    -.9353634542809267E+00,-.6312274766498699E+00, 
    0.9353634542809267E+00,0.6312274766498701E+00, 
    -.9674061994537118E-01,-.1150463641736545E+00, 
    0.9674061994537116E-01,0.1150463641736545E+00, 
    -.9438572336357420E+00,-.3747271550737232E+00, 
    0.9438572336357420E+00,0.3747271550737233E+00, 
    -.7904714944143627E+00,0.5008869009053424E+00, 
    0.7904714944143627E+00,-.5008869009053424E+00, 
    -.5513908795637696E+00,0.1608699923234913E+00, 
    0.5513908795637696E+00,-.1608699923234912E+00, 
    -.8210097847827826E+00,-.5657668007893579E+00, 
    0.8210097847827826E+00,0.5657668007893579E+00, 
    -.8970730060299638E+00,0.7712590259177344E+00, 
    0.8970730060299638E+00,-.7712590259177344E+00, 
    -.9113544664316420E+00,0.4163780011554896E+00, 
    0.9113544664316420E+00,-.4163780011554895E+00, 
    -.5579917371493338E+00,-.5836842069804584E+00, 
    0.5579917371493338E+00,0.5836842069804584E+00, 
    -.8254893061838994E+00,-.2215359986069006E+00, 
    0.8254893061838994E+00,0.2215359986069007E+00, 
    -.6957123574117627E+00,-.4009765829120038E+00, 
    0.6957123574117627E+00,0.4009765829120039E+00, 
    -.3763993540852892E+00,-.2343918950839542E-01, 
    0.3763993540852892E+00,0.2343918950839547E-01, 
    -.7179044391714098E+00,-.7105079759089011E+00, 
    0.7179044391714098E+00,0.7105079759089011E+00, 
    -.5431094815533654E+00,-.2200751569931638E+00, 
    0.5431094815533654E+00,0.2200751569931639E+00 };
  static double ws[] = {
    0.1227297430951934E-02,0.1227297430951934E-02, 
    0.1227297430951934E-02,0.1227297430951934E-02, 
    0.2705540778783016E-02,0.2705540778783016E-02, 
    0.2705540778783016E-02,0.2705540778783016E-02, 
    0.1737911403412906E-02,0.1737911403412906E-02, 
    0.1737911403412906E-02,0.1737911403412906E-02, 
    0.4392077971620324E-02,0.4392077971620324E-02, 
    0.4392077971620324E-02,0.4392077971620324E-02, 
    0.5022826333532844E-02,0.5022826333532844E-02, 
    0.5022826333532844E-02,0.5022826333532844E-02, 
    0.7874070819611865E-02,0.7874070819611865E-02, 
    0.7874070819611865E-02,0.7874070819611865E-02, 
    0.4436852998533855E-02,0.4436852998533855E-02, 
    0.4436852998533855E-02,0.4436852998533855E-02, 
    0.5267024709906290E-02,0.5267024709906290E-02, 
    0.5267024709906290E-02,0.5267024709906290E-02, 
    0.5674802470601283E-02,0.5674802470601283E-02, 
    0.5674802470601283E-02,0.5674802470601283E-02, 
    0.1214093209674204E-01,0.1214093209674204E-01, 
    0.1214093209674204E-01,0.1214093209674204E-01, 
    0.6621349971001042E-02,0.6621349971001042E-02, 
    0.6621349971001042E-02,0.6621349971001042E-02, 
    0.1231662091800690E-01,0.1231662091800690E-01, 
    0.1231662091800690E-01,0.1231662091800690E-01, 
    0.3418869012187221E-02,0.3418869012187221E-02, 
    0.3418869012187221E-02,0.3418869012187221E-02, 
    0.1030951784492058E-01,0.1030951784492058E-01, 
    0.1030951784492058E-01,0.1030951784492058E-01, 
    0.7362609227691093E-02,0.7362609227691093E-02, 
    0.7362609227691093E-02,0.7362609227691093E-02, 
    0.1132084139272881E-01,0.1132084139272881E-01, 
    0.1132084139272881E-01,0.1132084139272881E-01, 
    0.1579946795763282E-01,0.1579946795763282E-01, 
    0.1579946795763282E-01,0.1579946795763282E-01, 
    0.1762701971434555E-01,0.1762701971434555E-01, 
    0.1762701971434555E-01,0.1762701971434555E-01, 
    0.1445134560321702E-01,0.1445134560321702E-01, 
    0.1445134560321702E-01,0.1445134560321702E-01, 
    0.4104428749850738E-01,0.4104428749850738E-01, 
    0.4104428749850738E-01,0.4104428749850738E-01, 
    0.3339495242559154E-01,0.3339495242559154E-01, 
    0.3339495242559154E-01,0.3339495242559154E-01, 
    0.4184373496097664E-01,0.4184373496097664E-01, 
    0.4184373496097664E-01,0.4184373496097664E-01, 
    0.2495764669470555E-01,0.2495764669470555E-01, 
    0.2495764669470555E-01,0.2495764669470555E-01, 
    0.3707130727927116E-01,0.3707130727927116E-01, 
    0.3707130727927116E-01,0.3707130727927116E-01, 
    0.9609787849027714E-02,0.9609787849027714E-02, 
    0.9609787849027714E-02,0.9609787849027714E-02, 
    0.3355691943165765E-01,0.3355691943165765E-01, 
    0.3355691943165765E-01,0.3355691943165765E-01, 
    0.8580982368089280E-02,0.8580982368089280E-02, 
    0.8580982368089280E-02,0.8580982368089280E-02, 
    0.1755951049274731E-01,0.1755951049274731E-01, 
    0.1755951049274731E-01,0.1755951049274731E-01, 
    0.4298498004935582E-01,0.4298498004935582E-01, 
    0.4298498004935582E-01,0.4298498004935582E-01, 
    0.2200937149220828E-01,0.2200937149220828E-01, 
    0.2200937149220828E-01,0.2200937149220828E-01, 
    0.1430848314758049E-01,0.1430848314758049E-01, 
    0.1430848314758049E-01,0.1430848314758049E-01, 
    0.1667506010717100E-01,0.1667506010717100E-01, 
    0.1667506010717100E-01,0.1667506010717100E-01, 
    0.3444495645764802E-01,0.3444495645764802E-01, 
    0.3444495645764802E-01,0.3444495645764802E-01, 
    0.2897820038333119E-01,0.2897820038333119E-01, 
    0.2897820038333119E-01,0.2897820038333119E-01, 
    0.3492622185053053E-01,0.3492622185053053E-01, 
    0.3492622185053053E-01,0.3492622185053053E-01, 
    0.4881429507850928E-01,0.4881429507850928E-01, 
    0.4881429507850928E-01,0.4881429507850928E-01, 
    0.2338133304387865E-01,0.2338133304387865E-01, 
    0.2338133304387865E-01,0.2338133304387865E-01, 
    0.4325777192033273E-01,0.4325777192033273E-01, 
    0.4325777192033273E-01,0.4325777192033273E-01 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void rule30 ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    RULE30 returns the rule of degree 30.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int N, the number of nodes.
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  static double xs[] = {
    -.2489884108477549E+00,-.8955668996881347E+00, 
    -.9323001501704753E+00,0.7405445449992548E+00, 
    -.9340642229281805E+00,-.9095962664223526E+00, 
    0.5608633663495270E+00,0.8616511147917938E+00, 
    0.1660321315312064E+00,0.6574874238191415E+00, 
    -.1433381379425067E+00,0.1329569336990351E+00, 
    -.8829748195890315E+00,0.6632673829575200E+00, 
    0.4936163200514849E+00,-.2027212042708638E+00, 
    -.9958809661720823E+00,-.3073287548883581E+00, 
    -.8217303027626981E+00,0.8261783389660409E+00, 
    -.4397988250817157E+00,0.5010528086431736E+00, 
    0.5909729749119785E+00,-.3540551267912417E+00, 
    -.1341504650225714E+00,0.8797782958287345E+00, 
    0.7856276574168712E+00,0.2511140307013671E-01, 
    0.6793669279007186E+00,0.5990617334262452E+00, 
    -.6571411220532614E-01,0.2486338305340984E+00, 
    0.2517441361617084E+00,0.4956264496869945E+00, 
    0.6407808309272243E+00,0.3783356039497430E+00, 
    0.4110110511714637E+00,-.8395840683102380E+00, 
    -.8086877648916926E+00,0.1348481347528898E+00, 
    0.7911211318084024E-01,0.1744159174802842E+00, 
    0.3212567387892802E-01,-.9336794010864112E-02, 
    0.7322797884447179E+00,0.3357167967967313E-01, 
    0.5935826720491633E+00,0.8968431432339444E+00, 
    0.1799598313156058E+00,0.7269821230572000E+00, 
    -.1856842034490502E+00,0.6913593560408486E+00, 
    -.4397815672742683E+00,-.7758135419229234E+00, 
    -.9616538724444206E+00,0.3193596929830660E+00, 
    0.3777165903330749E+00,-.7679662364686681E+00, 
    -.7298701218806570E+00,0.4065079920371972E+00, 
    -.9227465688784468E-01,-.5950684639841661E+00, 
    -.9013518967720980E+00,0.8836909793143869E+00, 
    -.9967771339953647E+00,-.9944392441466685E+00, 
    -.2195956599465894E+00,-.4351118739522196E+00, 
    -.4349934895265168E+00,0.8129226389660974E+00, 
    0.9133061964092235E+00,-.7096773902564911E+00, 
    -.7019587221008244E-01,0.3087832511180271E+00, 
    -.4541723274189980E+00,0.1380191907138354E+00, 
    0.9551957279557286E+00,-.7407772548705425E+00, 
    -.6078105587102869E+00,-.8435179784664049E+00, 
    0.5490559994565902E+00,0.2444037338918357E+00, 
    0.2400810439649010E+00,-.5106820540745040E-01, 
    0.5284630091276127E+00,0.8294520284860720E+00, 
    0.9644131173686892E+00,0.5513347774529967E+00, 
    -.7997415011464921E+00,0.9063983320825396E-01, 
    0.5014531965594803E+00,-.3847230717105717E+00, 
    -.2781131563704615E+00,-.9344492367354887E+00, 
    -.9030917852729532E+00,-.4172241922741589E-01, 
    -.2640772153341865E+00,0.9665539194330073E+00, 
    -.2367623099741220E+00,-.8626350293444497E+00, 
    -.5924624188363723E+00,-.2560781746707642E+00, 
    -.9697346473592514E+00,-.5129714977732606E+00, 
    0.9006954577534996E+00,0.4485570029784398E+00, 
    -.9929762013100525E+00,-.4556636018911340E+00, 
    -.6294980126555781E+00,-.3125079163601901E+00, 
    -.7343713176297997E+00,-.6080741370996756E+00, 
    0.7969083548584857E+00,0.9691866547839980E+00, 
    0.3151361183804022E+00,0.3344055950815164E+00, 
    0.9472309841785602E+00,-.1247403647674209E+00, 
    -.9779732163470612E+00,-.7256364179379601E+00, 
    0.9970921358526692E+00,-.9954560187756241E+00, 
    0.9951087682452733E+00,-.6883166208655062E+00, 
    0.4290000273640588E+00,-.9880739319478465E+00, 
    0.9072764372199305E+00,0.9928546893888752E+00, 
    -.9917918025218534E+00,-.4805464361132344E+00, 
    -.8469642292616028E+00,-.8369560183233818E+00, 
    -.9779640479084950E+00,0.9687313192733549E+00, 
    -.7547772532654357E+00,0.7252641746474808E+00, 
    0.6646022965583286E+00,-.5644091993909820E+00, 
    0.8026432404348592E+00,0.9543957566085060E+00, 
    -.3557861147500951E+00,0.9951864328005721E+00, 
    -.9294765518337259E+00,-.9726603739446762E+00, 
    0.9908893344007184E+00,0.9593076671548199E+00, 
    -.9203911103942219E+00,-.5020036235062807E+00, 
    0.8473616010099635E+00,0.7880730531161636E+00, 
    0.1506248287699496E+00,0.9937769310041399E+00, 
    0.9894968458978942E+00,0.5755438442317007E+00, 
    0.9058856705910572E+00,0.9189816123793266E+00, 
    0.6951883785315084E+00,-.6234868878005112E+00, 
    -.9575663408305919E+00,-.9618986262036308E+00, 
    0.8232277101227666E+00,0.9715176930569319E+00, 
    -.6414127133627360E+00,-.9547454138171275E+00, 
    -.9923199594879274E+00,-.8900231674458331E+00, 
    0.9908836025186869E+00 };
  static double ys[] = {
    -.2460485077113066E+00,-.4965476117162992E+00, 
    -.6487707118457230E+00,0.3697665250429530E+00, 
    -.2243181000931358E+00,-.1272532317878749E+00, 
    0.6675821699023886E-01,-.3572065887354983E+00, 
    -.5859511801028319E+00,0.2705477662356954E+00, 
    -.1384585540894461E+00,-.4984429846834076E+00, 
    -.6385136736297631E+00,-.9972347998858430E+00, 
    0.9914800246598282E+00,0.7647568592516673E+00, 
    0.5617037705569573E+00,0.7101121364821149E+00, 
    -.9361722473508081E+00,0.4495309095930795E+00, 
    0.9980930920696178E+00,0.1901925696164289E+00, 
    0.9914568702475882E+00,-.3600355152827119E+00, 
    0.5732016635520147E+00,-.4984298059079483E+00, 
    -.2335476149195165E+00,0.6688792604932079E+00, 
    -.5195288082893980E-01,-.2801670391820981E+00, 
    0.8352476010432909E+00,-.7081859480565704E+00, 
    -.8979072023199524E+00,-.4437687184357049E+00, 
    -.5837384948634179E+00,0.3614240616611519E+00, 
    -.8230684981188165E+00,0.3172654777326074E-01, 
    -.9623528941703298E+00,0.7321217620214555E+00, 
    -.8096610148202545E+00,-.1586165019412888E+00, 
    0.2323785604104418E-01,-.3309606689238450E+00, 
    -.4388893082509347E+00,0.4126373875717694E+00, 
    0.7615965342755699E+00,0.2500590258181676E+00, 
    0.2336722577585314E+00,0.8396217552527814E+00, 
    -.5038771822936119E+00,0.5974085228316085E+00, 
    0.9608827453877556E+00,0.9958546191278792E+00, 
    -.4499429575894099E+00,0.5431833943398334E-01, 
    -.5768346162092601E+00,0.9424496317962137E+00, 
    0.2149562696494443E+00,-.9626198942478345E+00, 
    -.9006091929700322E+00,0.3935553572852833E+00, 
    -.9903441091784717E+00,-.9730346543555619E-01, 
    0.2076974492755955E+00,-.6015458100234288E+00, 
    -.7914402846211577E+00,-.9950700583024972E+00, 
    0.5506915379662010E+00,0.6946004094913842E+00, 
    0.5611276703738952E+00,-.9902901583614324E+00, 
    -.9923523090719414E+00,-.3380059761467528E+00, 
    0.2166811227803833E+00,0.9084360793204729E+00, 
    -.6099079065118447E+00,-.1613429998348757E+00, 
    0.2829088743783339E-01,-.3475035436069812E+00, 
    0.4843880673732963E+00,-.9945285505783277E+00, 
    0.5463996157581324E+00,0.9658858240647616E+00, 
    -.7106923954489180E+00,0.9126926988703539E+00, 
    0.9995360204936274E+00,-.9899747511065149E+00, 
    -.7645030064887491E+00,-.9580271636184000E+00, 
    0.8948233603717312E+00,-.8974768652542269E+00, 
    0.3866240780739671E+00,0.7812537795763034E+00, 
    -.8567508357551102E+00,-.6608100052801699E+00, 
    0.9104312794921889E+00,0.3794498701460580E+00, 
    0.9917517335327534E+00,0.8711471531031385E+00, 
    0.6865345494324926E+00,-.9676402895588128E+00, 
    -.7614862881147397E+00,-.7974426670548660E+00, 
    0.8007048860105550E+00,-.1342520685119810E+00, 
    -.2726893919561137E+00,0.8259904245622020E+00, 
    -.3535762642754187E+00,0.2570082752527717E-01, 
    0.5408483080025682E+00,0.9034516527900712E+00, 
    0.1094396192245524E+00,-.9898880600690686E+00, 
    0.8171118504095591E+00,0.9624825870234074E+00, 
    -.3058427697709690E+00,0.2042047511467314E+00, 
    0.6881014856247000E+00,0.8038659949683351E+00, 
    -.9483358470837810E+00,-.8774593970142823E+00, 
    0.8159614789170281E+00,-.8758831247154076E+00, 
    0.6649096212014286E+00,-.9861910838287823E+00, 
    0.9701981255999362E+00,-.7285580220078556E+00, 
    0.9697520026151714E+00,-.1654759673552738E+00, 
    0.3730458721291450E+00,0.6880375541948958E+00, 
    -.1463794215119850E-01,0.6931515112743655E+00, 
    -.5310794389210227E+00,-.9586758167636964E+00, 
    0.9527988373282164E+00,-.9549818827769936E+00, 
    -.6636898088045536E+00,0.4709151069275523E-01, 
    -.6618488641031038E+00,0.5367182637758624E+00, 
    0.1898490002080080E+00,0.3931950072226694E+00, 
    0.9627353314786278E+00,0.8959541550605267E+00, 
    0.5524608659159009E+00,-.5236758593851664E+00, 
    -.9898656273130529E+00,0.9909061133703810E+00, 
    0.9951465511067116E+00,0.1849675400643637E+00, 
    -.1455654267961535E+00,-.9037901015984503E+00, 
    -.7708757339972588E+00,-.9444760769791583E+00, 
    -.7980610860639374E+00,0.9837596096190924E+00, 
    -.9414770920443709E+00,0.9976903799315803E+00, 
    -.8814229711931103E+00,-.8637004475570209E+00, 
    -.6773031038124658E+00,0.9224881482845090E+00, 
    0.8455113672684367E+00,0.9760969611534657E+00, 
    -.4540561623942587E+00 };
  static double ws[] = {
    0.1656982041847443E-01,0.6241957476859499E-02, 
    0.5216414751791839E-02,0.1560407651121931E-01, 
    0.9922352722380101E-02,0.1044885486074365E-01, 
    0.1879260280925844E-01,0.7828025997082013E-02, 
    0.6646815879562722E-02,0.2535391408286021E-01, 
    0.3596968831636211E-01,0.3381660954529117E-01, 
    0.1200214214883407E-01,0.1160844683478245E-02, 
    0.3455740863097736E-02,0.1153404181974989E-01, 
    0.2326794602666712E-02,0.2497552316469688E-01, 
    0.6829794315843321E-02,0.1834140231596609E-01, 
    0.1335225445575319E-02,0.2522798332446686E-01, 
    0.1595934032110033E-02,0.3325015027182070E-01, 
    0.3467011496498743E-01,0.1611839351364714E-01, 
    0.2222188812186740E-01,0.1902350040458492E-01, 
    0.2953427938744245E-01,0.2950256499869300E-01, 
    0.2375703035217899E-01,0.2641584105673240E-01, 
    0.1690793382324671E-01,0.2648926396425241E-01, 
    0.2312689012710705E-01,0.3659471133442508E-01, 
    0.2037361988036312E-01,0.2428599858933432E-01, 
    0.2744225331234096E-02,0.2314679166528731E-01, 
    0.2431857706237597E-01,0.3987956826883379E-01, 
    0.4134270745058791E-01,0.4048978066628536E-01, 
    0.2251078653857245E-01,0.4125065697406070E-01, 
    0.2175419992100145E-01,0.1864439716615348E-01, 
    0.3968876145073878E-01,0.1521710713143089E-01, 
    0.3736662827386808E-01,0.2360382881992710E-01, 
    0.9970157722257745E-02,0.1786030245435301E-02, 
    0.1077371157065800E-01,0.3789581847120879E-01, 
    0.2650233252873021E-01,0.8678377005773784E-02, 
    0.3088113765197352E-01,0.1042906365431807E-01, 
    0.1836452487305106E-01,0.3403335573613761E-01, 
    0.2495937103585967E-02,0.1987508148017507E-01, 
    0.3001373057974097E-02,0.3414083443427366E-02, 
    0.2598635525971543E-01,0.3548865418196895E-02, 
    0.3429487078025957E-01,0.1785504393755270E-01, 
    0.1456023929593576E-01,0.3984363472067331E-02, 
    0.4363711971919699E-02,0.3712059926790805E-01, 
    0.3986478887569449E-01,0.1913617153849744E-01, 
    0.1049511173281814E-01,0.3060958199220223E-01, 
    0.3647607908842043E-01,0.2174227606463140E-01, 
    0.3117137400117779E-01,0.3796109966676785E-02, 
    0.3772008807242438E-01,0.1190095135082278E-01, 
    0.2275866415131214E-01,0.9507317753716060E-02, 
    0.6602435439732473E-03,0.3562338138346624E-02, 
    0.1790794872114957E-01,0.1166720090620105E-01, 
    0.1835513327961018E-01,0.1782419810205321E-01, 
    0.4039879332467389E-01,0.8544698960003014E-02, 
    0.1019822565974266E-01,0.3318239678979050E-01, 
    0.1843639498443942E-01,0.1049436999441544E-01, 
    0.4401594721967702E-02,0.1035186584479928E-01, 
    0.2774663878371582E-01,0.1038891259079699E-01, 
    0.7200083574761216E-02,0.2284618397210542E-01, 
    0.1143761979472999E-01,0.3626034532564225E-01, 
    0.4552809166608951E-02,0.2292163519800657E-01, 
    0.3339012139070554E-01,0.4305115781657949E-01, 
    0.2689701371924873E-01,0.1397250038349306E-01, 
    0.2657597849484922E-01,0.1419885923832320E-02, 
    0.2482810620945660E-01,0.1189134969439302E-01, 
    0.1376786292971261E-01,0.4412810694947523E-01, 
    0.5396384182228113E-02,0.1824530377035380E-01, 
    0.1045375639904588E-02,0.1889471014408118E-02, 
    0.2161857112903824E-02,0.1668110451401195E-01, 
    0.3023692617299402E-01,0.9232726779228300E-03, 
    0.4418205190510442E-02,0.3242974435474479E-02, 
    0.1103547225797253E-02,0.3921202767387279E-01, 
    0.2351260224707112E-01,0.1749878242380525E-01, 
    0.9478556031843440E-02,0.7823816828409005E-02, 
    0.2512073023506692E-01,0.8748252477568070E-02, 
    0.1020629440191134E-01,0.1142386634483088E-01, 
    0.2036353481009834E-01,0.1287466988683176E-01, 
    0.3056687396551737E-01,0.3145655833062538E-02, 
    0.1700461827659660E-01,0.9197009770615088E-02, 
    0.1441218157540074E-02,0.5481672105012487E-02, 
    0.1455840676178983E-01,0.3167421637373095E-01, 
    0.3035466841344861E-02,0.3337345728148663E-02, 
    0.4142211764571325E-02,0.4237586226174837E-02, 
    0.5970530714448440E-02,0.1594194686963999E-01, 
    0.1264147552673667E-01,0.5857851156194439E-02, 
    0.1967714628271757E-01,0.5157040758112765E-02, 
    0.4236825269003221E-02,0.7120803451144430E-03, 
    0.1238091173770127E-01,0.5529561908151683E-02, 
    0.2467962146273505E-01,0.4839010235595776E-02, 
    0.2559138088025228E-02,0.4161994856862944E-02, 
    0.4970131407881912E-02 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void square_arbq ( int degree, int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_ARBQ returns a quadrature rule for the symmetric square.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the degree of the quadrature rule.
//    1 <= DEGREE <= 30.
//
//    Input, int N, the number of nodes.
//    This can be determined by a call to SQUARE_ARBQ_SIZE(DEGREE).
//
//    Output, double X[2*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  int i;
  double w_sum;

  if ( degree == 1 )
  {
    rule01 ( n, x, w );
  }
  else if ( degree == 2 )
  {
    rule02 ( n, x, w );
  }
  else if ( degree == 3 )
  {
    rule03 ( n, x, w );
  }
  else if ( degree == 4 )
  {
    rule04 ( n, x, w );
  }
  else if ( degree == 5 )
  {
    rule05 ( n, x, w );
  }
  else if ( degree == 6 )
  {
    rule06 ( n, x, w );
  }
  else if ( degree == 7 )
  {
    rule07 ( n, x, w );
  }
  else if ( degree == 8 )
  {
    rule08 ( n, x, w );
  }
  else if ( degree == 9 )
  {
    rule09 ( n, x, w );
  }
  else if ( degree == 10 )
  {
    rule10 ( n, x, w );
  }
  else if ( degree == 11 )
  {
    rule11 ( n, x, w );
  }
  else if ( degree == 12 )
  {
    rule12 ( n, x, w );
  }
  else if ( degree == 13 )
  {
    rule13 ( n, x, w );
  }
  else if ( degree == 14 )
  {
    rule14 ( n, x, w );
  }
  else if ( degree == 15 )
  {
    rule15 ( n, x, w );
  }
  else if ( degree == 16 )
  {
    rule16 ( n, x, w );
  }
  else if ( degree == 17 )
  {
    rule17 ( n, x, w );
  }
  else if ( degree == 18 )
  {
    rule18 ( n, x, w );
  }
  else if ( degree == 19 )
  {
    rule19 ( n, x, w );
  }
  else if ( degree == 20 )
  {
    rule20 ( n, x, w );
  }
  else if ( degree == 21 )
  {
    rule21 ( n, x, w );
  }
  else if ( degree == 22 )
  {
    rule22 ( n, x, w );
  }
  else if ( degree == 23 )
  {
    rule23 ( n, x, w );
  }
  else if ( degree == 24 )
  {
    rule24 ( n, x, w );
  }
  else if ( degree == 25 )
  {
    rule25 ( n, x, w );
  }
  else if ( degree == 26 )
  {
    rule26 ( n, x, w );
  }
  else if ( degree == 27 )
  {
    rule27 ( n, x, w );
  }
  else if ( degree == 28 )
  {
    rule28 ( n, x, w );
  }
  else if ( degree == 29 )
  {
    rule29 ( n, x, w );
  }
  else if ( degree == 30 )
  {
    rule30 ( n, x, w );
  }
  else
  {
    cerr << "\n";
    cerr << "SQUARE_ARBQ - Fatal error\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  w_sum = r8vec_sum ( n, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = 4.0 * w[i] / w_sum;
  }

  return;
}
//****************************************************************************80

void square_arbq_gnuplot ( int n, double x[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_ARBQ_GNUPLOT: plot of a quadrature rule for the symmetric square.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    11 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, double N, the number of nodes.
//
//    Input, double X[2*N], the coordinates of the nodes.
//
//    Input, string HEADER, a string to be used to identify
//    the files created.
//
{
  string command_filename;
  ofstream command_unit;
  int i;
  int j;
  int l;
  string node_filename;
  ofstream node_unit;
  string plot_filename;
  string vertex_filename;
  ofstream vertex_unit;
//
//  Create the vertex file.
//
  vertex_filename = header + "_vertices.txt";
  vertex_unit.open ( vertex_filename.c_str ( ) );
  vertex_unit << -1.0E+00 << "  "
              << -1.0E+00 << "\n";
  vertex_unit << +1.0E+00 << "  "
              << -1.0E+00 << "\n";
  vertex_unit << +1.0E+00 << "  "
              << +1.0E+00 << "\n";
  vertex_unit << -1.0E+00 << "  "
              << +1.0E+00 << "\n";
  vertex_unit << -1.0E+00 << "  "
              << -1.0E+00 << "\n";
  vertex_unit.close ( );
  cout << "\n";
  cout << "  Created vertex file '" << vertex_filename << "'\n";
//
//  Create node file.
//
  node_filename = header + "_nodes.txt";
  node_unit.open ( node_filename.c_str ( ) );
  for ( j = 0; j < n; j++ )
  {
    node_unit << x[0+j*2] << "  "
              << x[1+j*2] << "\n";
  }
  node_unit.close ( );
  cout << "  Created node file '" << node_filename << "'\n";
//
//  Create graphics command file.
//
  command_filename = header + "_commands.txt";
  command_unit.open ( command_filename.c_str ( ) );
  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  plot_filename = header + ".png";
  command_unit << "set output '" << plot_filename << "'\n";
  command_unit << "set xlabel '<--- X --->'\n";
  command_unit << "set ylabel '<--- Y --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set size ratio -1\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "plot '" << vertex_filename << "' with lines lw 3, \\\n";
  command_unit << "     '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

int square_arbq_size ( int degree )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_ARBQ_SIZE: size of quadrature rule for the symmetric square.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2014
//
//  Author:
//
//    John Burkardt.
//
//  Reference:
//
//    Hong Xiao, Zydrunas Gimbutas,
//    A numerical algorithm for the construction of efficient quadrature
//    rules in two and higher dimensions,
//    Computers and Mathematics with Applications,
//    Volume 59, 2010, pages 663-676.
//
//  Parameters:
//
//    Input, int DEGREE, the desired degree of exactness.
//    1 <= DEGREE <= 30.
//
//    Output, int SQUARE_ARBQ_SIZE, the number of points in the
//    corresponding rule.
//
{
  int n;
  const int n_save[30] = {
      1,   3,   4,   6,   7, 
     10,  12,  16,  17,  22, 
     24,  31,  33,  41,  44, 
     52,  55,  64,  68,  78, 
     82,  93,  98, 109, 115, 
    127, 132, 147, 152, 167 };

  if ( degree < 1 || 30 < degree )
  {
    cerr << "\n";
    cerr << "SQUARE_ARBQ_SIZE - Fatal error!\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  n = n_save[degree-1];

  return n;
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
//    08 July 2009
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
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
