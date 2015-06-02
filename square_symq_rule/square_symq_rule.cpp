# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "square_symq_rule.hpp"

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
//    02 July 2014
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
//    02 July 2014
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

int rule_full_size ( int degree )

//****************************************************************************80
//
//  Purpose:
//
//    RULE_FULL_SIZE returns the full size of the requested quadrature rule.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    02 July 2014
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
//    Input, int DEGREE, the degree of the quadrature (the 
//    maximum degree of the polynomials of two variables that are integrated
//    exactly.  1 <= DEGREE <= 20.
//
//    Output, int RULE_FULL_SIZE, the number of nodes in the full rule.
//
{
  int n;
  const int n_save[20] = {
      1,   4,   4,   7,   7,  12,  12,  17,  17,  24, 
     24,  33,  33,  44,  44,  55,  55,  68,  68,  81 };

  if ( 1 <= degree && degree <= 20 )
  {
    n = n_save[degree-1];
  }
  else
  {
    cerr << "\n";
    cerr << "RULE_FULL_SIZE - Fatal error!\n";
    cerr << "  Degree DEGREE must be between 1 and 20.\n";
    exit ( 1 );
  }

  return n;
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
//    02 July 2014
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
       0.28284271247461904E+01 };

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
//    02 July 2014
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
  double xs[4] = { 
   -0.5773502691896256E+00, 
    0.5773502691896260E+00, 
    0.5773502691896256E+00, 
   -0.5773502691896260E+00 };
  double ys[4] = { 
   -0.5773502691896260E+00, 
   -0.5773502691896256E+00, 
    0.5773502691896260E+00, 
    0.5773502691896256E+00 };
  double ws[4] = { 
    0.7071067811865476E+00, 
    0.7071067811865476E+00, 
    0.7071067811865476E+00, 
    0.7071067811865476E+00 };

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
//    02 July 2014
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
  double xs[4] = { 
   -0.5773502691896256E+00, 
    0.5773502691896260E+00, 
    0.5773502691896256E+00, 
   -0.5773502691896260E+00 };
  double ys[4] = { 
   -0.5773502691896260E+00, 
   -0.5773502691896256E+00, 
    0.5773502691896260E+00, 
    0.5773502691896256E+00 };
  double ws[4] = { 
    0.7071067811865476E+00, 
    0.7071067811865476E+00, 
    0.7071067811865476E+00, 
    0.7071067811865476E+00 };

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
//    02 July 2014
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
  double xs[7] = { 
    0.3683480503448356E+00, 
   -0.3683480503448355E+00, 
    0.8881837963234579E+00, 
   -0.8881837963234579E+00, 
   -0.6849278434806340E+00, 
    0.6849278434806340E+00, 
    0.1035042199756803E-32 };
  double ys[7] = { 
   -0.8931142408116063E+00, 
    0.8931142408116063E+00, 
   -0.3800827242611582E+00, 
    0.3800827242611583E+00, 
   -0.6813275148988932E+00, 
    0.6813275148988932E+00, 
   -0.4874534345070689E-33 };
  double ws[7] = { 
    0.2922561796990344E+00, 
    0.2922561796990344E+00, 
    0.2970097006317383E+00, 
    0.2970097006317383E+00, 
    0.4208866642214383E+00, 
    0.4208866642214383E+00, 
    0.8081220356417685E+00 };

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
//    02 July 2014
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
  double xs[7] = { 
    0.1775868202077551E-01, 
   -0.1775868202077539E-01, 
    0.7788710544649639E+00, 
   -0.7788710544649639E+00, 
   -0.7703781288541645E+00, 
    0.7703781288541645E+00, 
   -0.7490353914168658E-33 };
  double ys[7] = { 
   -0.9659285494001192E+00, 
    0.9659285494001192E+00, 
   -0.5715708301251639E+00, 
    0.5715708301251639E+00, 
   -0.5829672991828014E+00, 
    0.5829672991828014E+00, 
    0.1356144833394667E-33 };
  double ws[7] = { 
    0.2246199725165690E+00, 
    0.2246199725165690E+00, 
    0.3901817339168917E+00, 
    0.3901817339168917E+00, 
    0.3953508381187504E+00, 
    0.3953508381187504E+00, 
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
//    02 July 2014
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
  double xs[12] = { 
    0.4595981103653579E-16, 
    0.9258200997725515E+00, 
    0.6742045114073804e-16, 
   -0.9258200997725515E+00, 
   -0.3805544332083157E+00, 
    0.3805544332083157E+00, 
    0.3805544332083157E+00, 
   -0.3805544332083157E+00, 
   -0.8059797829185990E+00, 
    0.8059797829185988E+00, 
    0.8059797829185990E+00, 
   -0.8059797829185988E+00 };
  double ys[12] = { 
   -0.9258200997725515E+00, 
   -0.1073032005210112E-16, 
    0.9258200997725515E+00, 
    0.1241105822293750e-15, 
   -0.3805544332083157E+00, 
   -0.3805544332083157E+00, 
    0.3805544332083157E+00, 
    0.3805544332083157E+00, 
   -0.8059797829185988E+00, 
   -0.8059797829185990E+00, 
    0.8059797829185988E+00, 
    0.8059797829185990E+00 };
  double ws[12] = { 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00 };

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
//    02 July 2014
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
  double xs[12] = { 
    0.4595981103653579E-16, 
    0.9258200997725515E+00, 
    0.6742045114073804E-16, 
   -0.9258200997725515E+00, 
   -0.3805544332083157E+00, 
    0.3805544332083157E+00, 
    0.3805544332083157E+00, 
   -0.3805544332083157E+00, 
   -0.8059797829185990E+00, 
    0.8059797829185988E+00, 
    0.8059797829185990E+00, 
   -0.8059797829185988E+00 };
  double ys[12] = { 
   -0.9258200997725515E+00, 
   -0.1073032005210112E-16, 
    0.9258200997725515E+00, 
    0.1241105822293750E-15, 
   -0.3805544332083157E+00, 
   -0.3805544332083157E+00, 
    0.3805544332083157E+00, 
    0.3805544332083157E+00, 
   -0.8059797829185988E+00, 
   -0.8059797829185990E+00, 
    0.8059797829185988E+00, 
    0.8059797829185990E+00 };
  double ws[12] = { 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.1711023816204485E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.3681147816131979E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00, 
    0.1678896179529011E+00 };

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
//    02 July 2014
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
  double xs[17] = { 
    0.6306801197316689E+00, 
    0.9688499663619776E+00, 
   -0.6306801197316687E+00, 
   -0.9688499663619776E+00, 
   -0.7502770999789002E+00, 
    0.9279616459595696E+00, 
    0.7502770999789005E+00, 
   -0.9279616459595696E+00, 
   -0.7620832819261708E-01, 
    0.8526157293336623E+00, 
    0.7620832819261719E-01, 
   -0.8526157293336623E+00, 
   -0.5237358202144292E+00, 
    0.4533398211356472E+00, 
    0.5237358202144292E+00, 
   -0.4533398211356471E+00, 
    0.1018964154952896E-32 };
  double ys[17] = { 
   -0.9688499663619776E+00, 
    0.6306801197316688E+00, 
    0.9688499663619776E+00, 
   -0.6306801197316686E+00, 
   -0.9279616459595696E+00, 
   -0.7502770999789004E+00, 
    0.9279616459595696E+00, 
    0.7502770999789006E+00, 
   -0.8526157293336623E+00, 
   -0.7620832819261714E-01, 
    0.8526157293336623E+00, 
    0.7620832819261725E-01, 
   -0.4533398211356472E+00, 
   -0.5237358202144292E+00, 
    0.4533398211356471E+00, 
    0.5237358202144292E+00, 
   -0.7403196379681869E-32 };
  double ws[17] = { 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.3724677695139019E+00 };

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
//    02 July 2014
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
  double xs[17] = { 
    0.6306801197316689E+00, 
    0.9688499663619776E+00, 
   -0.6306801197316687E+00, 
   -0.9688499663619776E+00, 
   -0.7502770999789002E+00, 
    0.9279616459595696E+00, 
    0.7502770999789005E+00, 
   -0.9279616459595696E+00, 
   -0.7620832819261708E-01, 
    0.8526157293336623E+00, 
    0.7620832819261719E-01, 
   -0.8526157293336623E+00, 
   -0.5237358202144292E+00, 
    0.4533398211356472E+00, 
    0.5237358202144292E+00, 
   -0.4533398211356471E+00, 
    0.1018964154952896E-32 };
  double ys[17] = { 
   -0.9688499663619776E+00, 
    0.6306801197316688E+00, 
    0.9688499663619776E+00, 
   -0.6306801197316686E+00, 
   -0.9279616459595696E+00, 
   -0.7502770999789004E+00, 
    0.9279616459595696E+00, 
    0.7502770999789006E+00, 
   -0.8526157293336623E+00, 
   -0.7620832819261714E-01, 
    0.8526157293336623E+00, 
    0.7620832819261725E-01, 
   -0.4533398211356472E+00, 
   -0.5237358202144292E+00, 
    0.4533398211356471E+00, 
    0.5237358202144292E+00, 
   -0.7403196379681869E-32 };
  double ws[17] = { 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.6284721101179121E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.7926638883415160E-01, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.1902480253324004E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.2816282136297291E+00, 
    0.3724677695139019E+00 };

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
//    02 July 2014
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
  double xs[24] = { 
   -0.6980761045495689E+00, 
    0.9826392235408554E+00, 
    0.6980761045495691E+00, 
   -0.9826392235408554E+00, 
    0.8257758359029634E+00, 
    0.9394863828167371E+00, 
   -0.8257758359029632E+00, 
   -0.9394863828167371E+00, 
    0.1885861387186400E+00, 
    0.9535395282015321E+00, 
   -0.1885861387186399E+00, 
   -0.9535395282015321E+00, 
   -0.7120019130753369E+00, 
    0.5253202503645465E+00, 
    0.7120019130753369E+00, 
   -0.5253202503645465E+00, 
   -0.3156234329152560E+00, 
    0.8125205483048131E+00, 
    0.3156234329152561E+00, 
   -0.8125205483048131E+00, 
   -0.4248472488486695E+00, 
    0.4165807191202114E-01, 
    0.4248472488486695E+00, 
   -0.4165807191202109E-01 };
  double ys[24] = { 
   -0.9826392235408554E+00, 
   -0.6980761045495690E+00, 
    0.9826392235408554E+00, 
    0.6980761045495693E+00, 
   -0.9394863828167371E+00, 
    0.8257758359029633E+00, 
    0.9394863828167371E+00, 
   -0.8257758359029631E+00, 
   -0.9535395282015321E+00, 
    0.1885861387186400E+00, 
    0.9535395282015321E+00, 
   -0.1885861387186399E+00, 
   -0.5253202503645465E+00, 
   -0.7120019130753369E+00, 
    0.5253202503645465E+00, 
    0.7120019130753369E+00, 
   -0.8125205483048131E+00, 
   -0.3156234329152560E+00, 
    0.8125205483048131E+00, 
    0.3156234329152561E+00, 
   -0.4165807191202117E-01, 
   -0.4248472488486695E+00, 
    0.4165807191202112E-01, 
    0.4248472488486695E+00 };
  double ws[24] = { 
    0.3395580740305119E-01, 
    0.3395580740305119E-01, 
    0.3395580740305119E-01, 
    0.3395580740305119E-01, 
    0.4671948489426219E-01, 
    0.4671948489426219E-01, 
    0.4671948489426219E-01, 
    0.4671948489426219E-01, 
    0.6886285066821875E-01, 
    0.6886285066821875E-01, 
    0.6886285066821875E-01, 
    0.6886285066821875E-01, 
    0.1595417182608940E+00, 
    0.1595417182608940E+00, 
    0.1595417182608940E+00, 
    0.1595417182608940E+00, 
    0.1497202089079447E+00, 
    0.1497202089079447E+00, 
    0.1497202089079447E+00, 
    0.1497202089079447E+00, 
    0.2483067110521768E+00, 
    0.2483067110521768E+00, 
    0.2483067110521768E+00, 
    0.2483067110521768E+00 };

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
//    02 July 2014
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
  double xs[24] = { 
    0.1885861387186414E+00, 
    0.9535395282015320E+00, 
   -0.1885861387186413E+00, 
   -0.9535395282015320E+00, 
   -0.6980761045495679E+00, 
    0.9826392235408555E+00, 
    0.6980761045495681E+00, 
   -0.9826392235408555E+00, 
   -0.9394863828167370E+00, 
    0.8257758359029639E+00, 
    0.9394863828167370E+00, 
   -0.8257758359029637E+00, 
   -0.7120019130753364E+00, 
    0.5253202503645475E+00, 
    0.7120019130753364E+00, 
   -0.5253202503645475E+00, 
   -0.3156234329152547E+00, 
    0.8125205483048131E+00, 
    0.3156234329152548E+00, 
   -0.8125205483048131E+00, 
   -0.4248472488486694E+00, 
    0.4165807191202203E-01, 
    0.4248472488486694E+00, 
   -0.4165807191202197E-01 };
  double ys[24] = { 
   -0.9535395282015320E+00, 
    0.1885861387186414E+00, 
    0.9535395282015320E+00, 
   -0.1885861387186413E+00, 
   -0.9826392235408555E+00, 
   -0.6980761045495680E+00, 
    0.9826392235408555E+00, 
    0.6980761045495683E+00, 
   -0.8257758359029640E+00, 
   -0.9394863828167370E+00, 
    0.8257758359029638E+00, 
    0.9394863828167370E+00, 
   -0.5253202503645475E+00, 
   -0.7120019130753364E+00, 
    0.5253202503645475E+00, 
    0.7120019130753364E+00, 
   -0.8125205483048131E+00, 
   -0.3156234329152547E+00, 
    0.8125205483048131E+00, 
    0.3156234329152549E+00, 
   -0.4165807191202205E-01, 
   -0.4248472488486694E+00, 
    0.4165807191202200E-01, 
    0.4248472488486694E+00 };
  double ws[24] = { 
    0.6886285066821880E-01, 
    0.6886285066821880E-01, 
    0.6886285066821880E-01, 
    0.6886285066821880E-01, 
    0.3395580740305121E-01, 
    0.3395580740305121E-01, 
    0.3395580740305121E-01, 
    0.3395580740305121E-01, 
    0.4671948489426224E-01, 
    0.4671948489426224E-01, 
    0.4671948489426224E-01, 
    0.4671948489426224E-01, 
    0.1595417182608939E+00, 
    0.1595417182608939E+00, 
    0.1595417182608939E+00, 
    0.1595417182608939E+00, 
    0.1497202089079448E+00, 
    0.1497202089079448E+00, 
    0.1497202089079448E+00, 
    0.1497202089079448E+00, 
    0.2483067110521767E+00, 
    0.2483067110521767E+00, 
    0.2483067110521767E+00, 
    0.2483067110521767E+00 };

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
//    02 July 2014
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
  double xs[33] = { 
   -0.9572976997863073E+00, 
    0.8595560056416388E+00, 
    0.9572976997863073E+00, 
   -0.8595560056416386E+00, 
   -0.7788097115544194E+00, 
    0.9834866824398721E+00, 
    0.7788097115544196E+00, 
   -0.9834866824398721E+00, 
   -0.4758086252182758E+00, 
    0.8500766736997486E+00, 
    0.4758086252182759E+00, 
   -0.8500766736997486E+00, 
   -0.7558053565720815E+00, 
    0.6478216371870107E+00, 
    0.7558053565720815E+00, 
   -0.6478216371870107E+00, 
   -0.3427165560404068E+00, 
    0.4093045616940387E+00, 
    0.3427165560404068E+00, 
   -0.4093045616940387E+00, 
   -0.1381834598624653E+00, 
    0.9589251702875349E+00, 
    0.1381834598624654E+00, 
   -0.9589251702875349E+00, 
    0.7074150899644485E-01, 
    0.6962500784917494E+00, 
   -0.7074150899644477E-01, 
   -0.6962500784917494E+00, 
    0.3907362161294610E+00, 
    0.9413272258729252E+00, 
   -0.3907362161294609E+00, 
   -0.9413272258729252E+00, 
   -0.3126032252245169E-31 };
  double ys[33] = { 
   -0.8595560056416389E+00, 
   -0.9572976997863073E+00, 
    0.8595560056416387E+00, 
    0.9572976997863073E+00, 
   -0.9834866824398721E+00, 
   -0.7788097115544195E+00, 
    0.9834866824398721E+00, 
    0.7788097115544197E+00, 
   -0.8500766736997486E+00, 
   -0.4758086252182758E+00, 
    0.8500766736997486E+00, 
    0.4758086252182759E+00, 
   -0.6478216371870107E+00, 
   -0.7558053565720815E+00, 
    0.6478216371870107E+00, 
    0.7558053565720815E+00, 
   -0.4093045616940387E+00, 
   -0.3427165560404068E+00, 
    0.4093045616940387E+00, 
    0.3427165560404068E+00, 
   -0.9589251702875349E+00, 
   -0.1381834598624653E+00, 
    0.9589251702875349E+00, 
    0.1381834598624654E+00, 
   -0.6962500784917494E+00, 
    0.7074150899644481E-01, 
    0.6962500784917494E+00, 
   -0.7074150899644473E-01, 
   -0.9413272258729252E+00, 
    0.3907362161294610E+00, 
    0.9413272258729252E+00, 
   -0.3907362161294609E+00, 
   -0.1114446878059780E-31 };
  double ws[33] = { 
    0.2699339218118220E-01, 
    0.2699339218118220E-01, 
    0.2699339218118220E-01, 
    0.2699339218118220E-01, 
    0.2120743264134161E-01, 
    0.2120743264134161E-01, 
    0.2120743264134161E-01, 
    0.2120743264134161E-01, 
    0.8403587015611028E-01, 
    0.8403587015611028E-01, 
    0.8403587015611028E-01, 
    0.8403587015611028E-01, 
    0.9175668641747105E-01, 
    0.9175668641747105E-01, 
    0.9175668641747105E-01, 
    0.9175668641747105E-01, 
    0.1816350488471703E+00, 
    0.1816350488471703E+00, 
    0.1816350488471703E+00, 
    0.1816350488471703E+00, 
    0.4272687338421145E-01, 
    0.4272687338421145E-01, 
    0.4272687338421145E-01, 
    0.4272687338421145E-01, 
    0.1508552789574408E+00, 
    0.1508552789574408E+00, 
    0.1508552789574408E+00, 
    0.1508552789574408E+00, 
    0.5479564090947486E-01, 
    0.5479564090947486E-01, 
    0.5479564090947486E-01, 
    0.5479564090947486E-01, 
    0.2124022307685798E+00 };

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
//    02 July 2014
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
  double xs[33] = { 
   -0.9572976997863074E+00, 
    0.8595560056416388E+00, 
    0.9572976997863074E+00, 
   -0.8595560056416386E+00, 
   -0.7788097115544195E+00, 
    0.9834866824398722E+00, 
    0.7788097115544197E+00, 
   -0.9834866824398722E+00, 
   -0.4758086252182752E+00, 
    0.8500766736997490E+00, 
    0.4758086252182753E+00, 
   -0.8500766736997490E+00, 
    0.3907362161294613E+00, 
    0.9413272258729251E+00, 
   -0.3907362161294612E+00, 
   -0.9413272258729251E+00, 
   -0.1381834598624646E+00, 
    0.9589251702875351E+00, 
    0.1381834598624647E+00, 
   -0.9589251702875351E+00, 
    0.6478216371870111E+00, 
    0.7558053565720809E+00, 
   -0.6478216371870111E+00, 
   -0.7558053565720809E+00, 
    0.7074150899644462E-01, 
    0.6962500784917495E+00, 
   -0.7074150899644453E-01, 
   -0.6962500784917495E+00, 
   -0.3427165560404070E+00, 
    0.4093045616940387E+00, 
    0.3427165560404070E+00, 
   -0.4093045616940387E+00, 
   -0.7375869198366919E-30 };
  double ys[33] = { 
   -0.8595560056416389E+00, 
   -0.9572976997863074E+00, 
    0.8595560056416387E+00, 
    0.9572976997863074E+00, 
   -0.9834866824398722E+00, 
   -0.7788097115544196E+00, 
    0.9834866824398722E+00, 
    0.7788097115544198E+00, 
   -0.8500766736997490E+00, 
   -0.4758086252182752E+00, 
    0.8500766736997490E+00, 
    0.4758086252182753E+00, 
   -0.9413272258729251E+00, 
    0.3907362161294612E+00, 
    0.9413272258729251E+00, 
   -0.3907362161294611E+00, 
   -0.9589251702875351E+00, 
   -0.1381834598624647E+00, 
    0.9589251702875351E+00, 
    0.1381834598624648E+00, 
   -0.7558053565720809E+00, 
    0.6478216371870111E+00, 
    0.7558053565720809E+00, 
   -0.6478216371870111E+00, 
   -0.6962500784917495E+00, 
    0.7074150899644457E-01, 
    0.6962500784917495E+00, 
   -0.7074150899644449E-01, 
   -0.4093045616940387E+00, 
   -0.3427165560404070E+00, 
    0.4093045616940387E+00, 
    0.3427165560404070E+00, 
   -0.6522588594679827E-30 };
  double ws[33] = { 
    0.2699339218118215E-01, 
    0.2699339218118215E-01, 
    0.2699339218118215E-01, 
    0.2699339218118215E-01, 
    0.2120743264134157E-01, 
    0.2120743264134157E-01, 
    0.2120743264134157E-01, 
    0.2120743264134157E-01, 
    0.8403587015611026E-01, 
    0.8403587015611026E-01, 
    0.8403587015611026E-01, 
    0.8403587015611026E-01, 
    0.5479564090947502E-01, 
    0.5479564090947502E-01, 
    0.5479564090947502E-01, 
    0.5479564090947502E-01, 
    0.4272687338421139E-01, 
    0.4272687338421139E-01, 
    0.4272687338421139E-01, 
    0.4272687338421139E-01, 
    0.9175668641747110E-01, 
    0.9175668641747110E-01, 
    0.9175668641747110E-01, 
    0.9175668641747110E-01, 
    0.1508552789574409E+00, 
    0.1508552789574409E+00, 
    0.1508552789574409E+00, 
    0.1508552789574409E+00, 
    0.1816350488471704E+00, 
    0.1816350488471704E+00, 
    0.1816350488471704E+00, 
    0.1816350488471704E+00, 
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
//    02 July 2014
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
  double xs[44] = { 
   -0.6714783701550190E+00, 
    0.9859876542016408E+00, 
    0.6714783701550192E+00, 
   -0.9859876542016408E+00, 
   -0.9318844245957986E+00, 
    0.9382770335701854E+00, 
    0.9318844245957988E+00, 
   -0.9382770335701852E+00, 
    0.6776977793098985E+00, 
    0.9773357693271729E+00, 
   -0.6776977793098983E+00, 
   -0.9773357693271729E+00, 
    0.4073679548284153E+00, 
    0.8648066658739809E+00, 
   -0.4073679548284151E+00, 
   -0.8648066658739809E+00, 
    0.6518175069036650E-01, 
    0.9759935658724420E+00, 
   -0.6518175069036639E-01, 
   -0.9759935658724420E+00, 
   -0.7473119631960774E+00, 
    0.7834652444128232E+00, 
    0.7473119631960774E+00, 
   -0.7834652444128232E+00, 
    0.1328305205898269E+00, 
    0.6241210323620054E+00, 
   -0.1328305205898269E+00, 
   -0.6241210323620054E+00, 
   -0.4781379108769722E+00, 
    0.5501448214169192E+00, 
    0.4781379108769723E+00, 
   -0.5501448214169192E+00, 
   -0.1803286643164523E+00, 
    0.8053335984690123E+00, 
    0.1803286643164524E+00, 
   -0.8053335984690123E+00, 
   -0.4134760830488010E+00, 
    0.9261965849285028E+00, 
    0.4134760830488011E+00, 
   -0.9261965849285028E+00, 
   -0.1307639250027494E+00, 
    0.2910908755606336E+00, 
    0.1307639250027494E+00, 
   -0.2910908755606336E+00 };
  double ys[44] = { 
   -0.9859876542016408E+00, 
   -0.6714783701550191E+00, 
    0.9859876542016408E+00, 
    0.6714783701550193E+00, 
   -0.9382770335701855E+00, 
   -0.9318844245957987E+00, 
    0.9382770335701853E+00, 
    0.9318844245957989E+00, 
   -0.9773357693271729E+00, 
    0.6776977793098984E+00, 
    0.9773357693271729E+00, 
   -0.6776977793098982E+00, 
   -0.8648066658739809E+00, 
    0.4073679548284152E+00, 
    0.8648066658739809E+00, 
   -0.4073679548284151E+00, 
   -0.9759935658724420E+00, 
    0.6518175069036644E-01, 
    0.9759935658724420E+00, 
   -0.6518175069036633E-01, 
   -0.7834652444128232E+00, 
   -0.7473119631960774E+00, 
    0.7834652444128232E+00, 
    0.7473119631960774E+00, 
   -0.6241210323620054E+00, 
    0.1328305205898269E+00, 
    0.6241210323620054E+00, 
   -0.1328305205898269E+00, 
   -0.5501448214169192E+00, 
   -0.4781379108769723E+00, 
    0.5501448214169192E+00, 
    0.4781379108769724E+00, 
   -0.8053335984690123E+00, 
   -0.1803286643164524E+00, 
    0.8053335984690123E+00, 
    0.1803286643164525E+00, 
   -0.9261965849285028E+00, 
   -0.4134760830488011E+00, 
    0.9261965849285028E+00, 
    0.4134760830488012E+00, 
   -0.2910908755606336E+00, 
   -0.1307639250027494E+00, 
    0.2910908755606336E+00, 
    0.1307639250027494E+00 };
  double ws[44] = { 
    0.1410384661573933E-01, 
    0.1410384661573933E-01, 
    0.1410384661573933E-01, 
    0.1410384661573933E-01, 
    0.1896652423210582E-01, 
    0.1896652423210582E-01, 
    0.1896652423210582E-01, 
    0.1896652423210582E-01, 
    0.2088141025507279E-01, 
    0.2088141025507279E-01, 
    0.2088141025507279E-01, 
    0.2088141025507279E-01, 
    0.7331394692154988E-01, 
    0.7331394692154988E-01, 
    0.7331394692154988E-01, 
    0.7331394692154988E-01, 
    0.3078002143226069E-01, 
    0.3078002143226069E-01, 
    0.3078002143226069E-01, 
    0.3078002143226069E-01, 
    0.6693059666394105E-01, 
    0.6693059666394105E-01, 
    0.6693059666394105E-01, 
    0.6693059666394105E-01, 
    0.1122840307920054E+00, 
    0.1122840307920054E+00, 
    0.1122840307920054E+00, 
    0.1122840307920054E+00, 
    0.1159261595200915E+00, 
    0.1159261595200915E+00, 
    0.1159261595200915E+00, 
    0.1159261595200915E+00, 
    0.7346051498025349E-01, 
    0.7346051498025349E-01, 
    0.7346051498025349E-01, 
    0.7346051498025349E-01, 
    0.4099703937729331E-01, 
    0.4099703937729331E-01, 
    0.4099703937729331E-01, 
    0.4099703937729331E-01, 
    0.1394626903962344E+00, 
    0.1394626903962344E+00, 
    0.1394626903962344E+00, 
    0.1394626903962344E+00 };

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
//    02 July 2014
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
  double xs[44] = { 
    0.7749527857778351E+00, 
    0.9885448991378063E+00, 
   -0.7749527857778349E+00, 
   -0.9885448991378063E+00, 
   -0.9070374303651182E+00, 
    0.9571446613308432E+00, 
    0.9070374303651184E+00, 
   -0.9571446613308430E+00, 
   -0.4303978306869286E+00, 
    0.9769578054468787E+00, 
    0.4303978306869287E+00, 
   -0.9769578054468787E+00, 
   -0.9756646723906326E+00, 
    0.1107064048513496E+00, 
    0.9756646723906326E+00, 
   -0.1107064048513495E+00, 
   -0.7388921437312957E+00, 
    0.7868610204187212E+00, 
    0.7388921437312957E+00, 
   -0.7868610204187212E+00, 
    0.1995220876718269E+00, 
    0.6659287668239546E+00, 
   -0.1995220876718268E+00, 
   -0.6659287668239546E+00, 
   -0.1934983412061240E+00, 
    0.8412271039808018E+00, 
    0.1934983412061241E+00, 
   -0.8412271039808018E+00, 
    0.4882189227791580E+00, 
    0.8922368778153702E+00, 
   -0.4882189227791579E+00, 
   -0.8922368778153702E+00, 
   -0.5772265461040059E+00, 
    0.9526539504944950E+00, 
    0.5772265461040061E+00, 
   -0.9526539504944950E+00, 
   -0.4474426063114782E+00, 
    0.5675455860909890E+00, 
    0.4474426063114783E+00, 
   -0.5675455860909890E+00, 
   -0.7044956995149931E-01, 
    0.3256679896817100E+00, 
    0.7044956995149934E-01, 
   -0.3256679896817100E+00 };
  double ys[44] = { 
   -0.9885448991378063E+00, 
    0.7749527857778350E+00, 
    0.9885448991378063E+00, 
   -0.7749527857778348E+00, 
   -0.9571446613308433E+00, 
   -0.9070374303651183E+00, 
    0.9571446613308431E+00, 
    0.9070374303651185E+00, 
   -0.9769578054468787E+00, 
   -0.4303978306869286E+00, 
    0.9769578054468787E+00, 
    0.4303978306869287E+00, 
   -0.1107064048513496E+00, 
   -0.9756646723906326E+00, 
    0.1107064048513495E+00, 
    0.9756646723906326E+00, 
   -0.7868610204187212E+00, 
   -0.7388921437312957E+00, 
    0.7868610204187212E+00, 
    0.7388921437312957E+00, 
   -0.6659287668239546E+00, 
    0.1995220876718268E+00, 
    0.6659287668239546E+00, 
   -0.1995220876718268E+00, 
   -0.8412271039808018E+00, 
   -0.1934983412061240E+00, 
    0.8412271039808018E+00, 
    0.1934983412061241E+00, 
   -0.8922368778153702E+00, 
    0.4882189227791580E+00, 
    0.8922368778153702E+00, 
   -0.4882189227791578E+00, 
   -0.9526539504944950E+00, 
   -0.5772265461040060E+00, 
    0.9526539504944950E+00, 
    0.5772265461040063E+00, 
   -0.5675455860909890E+00, 
   -0.4474426063114783E+00, 
    0.5675455860909890E+00, 
    0.4474426063114784E+00, 
   -0.3256679896817100E+00, 
   -0.7044956995149933E-01, 
    0.3256679896817100E+00, 
    0.7044956995149936E-01 };
  double ws[44] = { 
    0.1443015463807196E-01, 
    0.1443015463807196E-01, 
    0.1443015463807196E-01, 
    0.1443015463807196E-01, 
    0.1816242330920956E-01, 
    0.1816242330920956E-01, 
    0.1816242330920956E-01, 
    0.1816242330920956E-01, 
    0.1290815898308381E-01, 
    0.1290815898308381E-01, 
    0.1290815898308381E-01, 
    0.1290815898308381E-01, 
    0.3010764365372140E-01, 
    0.3010764365372140E-01, 
    0.3010764365372140E-01, 
    0.3010764365372140E-01, 
    0.6540469907131932E-01, 
    0.6540469907131932E-01, 
    0.6540469907131932E-01, 
    0.6540469907131932E-01, 
    0.1197895531736646E+00, 
    0.1197895531736646E+00, 
    0.1197895531736646E+00, 
    0.1197895531736646E+00, 
    0.8473841548096289E-01, 
    0.8473841548096289E-01, 
    0.8473841548096289E-01, 
    0.8473841548096289E-01, 
    0.6453833756714425E-01, 
    0.6453833756714425E-01, 
    0.6453833756714425E-01, 
    0.6453833756714425E-01, 
    0.2403055376316494E-01, 
    0.2403055376316494E-01, 
    0.2403055376316494E-01, 
    0.2403055376316494E-01, 
    0.1196130510491228E+00, 
    0.1196130510491228E+00, 
    0.1196130510491228E+00, 
    0.1196130510491228E+00, 
    0.1533837904970821E+00, 
    0.1533837904970821E+00, 
    0.1533837904970821E+00, 
    0.1533837904970821E+00 };

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
//    02 July 2014
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
  double xs[55] = { 
    0.7331873192446229E+00, 
   -0.7331873192446227E+00, 
   -0.9811278880414770E+00, 
    0.9811278880414772E+00, 
   -0.8004995596996590E+00, 
    0.8004995596996592E+00, 
    0.2935594202060772E+00, 
   -0.2935594202060772E+00, 
    0.5019013651861420E+00, 
   -0.5019013651861418E+00, 
   -0.9240427888147712E+00, 
    0.9240427888147712E+00, 
   -0.7321159842417640E+00, 
    0.7321159842417640E+00, 
    0.9107218705094187E+00, 
   -0.9107218705094184E+00, 
    0.9799531606782582E+00, 
   -0.9799531606782582E+00, 
   -0.2536359436096021E+00, 
    0.2536359436096021E+00, 
    0.8800049697526030E+00, 
   -0.8800049697526030E+00, 
    0.7136219272623606E+00, 
   -0.7136219272623606E+00, 
    0.5185051092186185E+00, 
   -0.5185051092186185E+00, 
    0.9890262305049052E+00, 
   -0.9890262305049052E+00, 
    0.9865971248382277E+00, 
   -0.9865971248382277E+00, 
    0.4087785918187709E-01, 
   -0.4087785918187702E-01, 
    0.9650604144351506E+00, 
   -0.9650604144351506E+00, 
   -0.5228670170578392E+00, 
    0.5228670170578394E+00, 
   -0.2304316370092423E+00, 
    0.2304316370092424E+00, 
    0.7381821882135022E+00, 
   -0.7381821882135022E+00, 
   -0.4979206093242921E+00, 
    0.4979206093242922E+00, 
    0.8494669121845019E+00, 
   -0.8494669121845019E+00, 
    0.4390176422841324E+00, 
   -0.4390176422841323E+00, 
    0.1590601194183188E+00, 
   -0.1590601194183187E+00, 
    0.8973818517920210E+00, 
   -0.8973818517920210E+00, 
    0.6726312443333152E+00, 
   -0.6726312443333152E+00, 
   -0.1686064273871127E+00, 
    0.1686064273871128E+00, 
   -0.3548241530243386E-18 };
  double ys[55] = { 
   -0.9711078221435576E+00, 
    0.9711078221435576E+00, 
   -0.9668551959097115E+00, 
    0.9668551959097113E+00, 
   -0.9746926011666336E+00, 
    0.9746926011666336E+00, 
   -0.3231309208576288E+00, 
    0.3231309208576288E+00, 
   -0.9765444785368099E+00, 
    0.9765444785368099E+00, 
   -0.8490306235166675E+00, 
    0.8490306235166672E+00, 
   -0.7537198042004623E+00, 
    0.7537198042004623E+00, 
   -0.9737587969123404E+00, 
    0.9737587969123406E+00, 
   -0.3822148312292263E+00, 
    0.3822148312292264E+00, 
   -0.2988363050086515E+00, 
    0.2988363050086515E+00, 
    0.4849608774128832E+00, 
   -0.4849608774128831E+00, 
    0.2492237020321146E+00, 
   -0.2492237020321144E+00, 
   -0.3504141436316342E-01, 
    0.3504141436316349E-01, 
    0.6278936489285102E+00, 
   -0.6278936489285100E+00, 
   -0.8591476119499135E+00, 
    0.8591476119499137E+00, 
   -0.5892598635566724E+00, 
    0.5892598635566724E+00, 
    0.1438346146728415E+00, 
   -0.1438346146728414E+00, 
   -0.9289486752701194E+00, 
    0.9289486752701194E+00, 
   -0.8028060773786958E+00, 
    0.8028060773786958E+00, 
   -0.8651144139342870E+00, 
    0.8651144139342870E+00, 
   -0.5653829126627348E+00, 
    0.5653829126627348E+00, 
   -0.1574661586091270E+00, 
    0.1574661586091272E+00, 
   -0.7312745784466166E+00, 
    0.7312745784466166E+00, 
   -0.9115177107109407E+00, 
    0.9115177107109407E+00, 
   -0.6626783799774586E+00, 
    0.6626783799774586E+00, 
   -0.4696061222418765E+00, 
    0.4696061222418766E+00, 
   -0.9939228673343959E+00, 
    0.9939228673343959E+00, 
    0.3228625474392587E-19 };
  double ws[55] = { 
    0.3224472434909546E-02, 
    0.3224472434909546E-02, 
    0.4080157527921578E-02, 
    0.4080157527921578E-02, 
    0.1406321867924724E-01, 
    0.1406321867924724E-01, 
    0.1094478053582958E+00, 
    0.1094478053582958E+00, 
    0.2046021623250057E-01, 
    0.2046021623250057E-01, 
    0.2244481290183435E-01, 
    0.2244481290183435E-01, 
    0.5310357585578484E-01, 
    0.5310357585578484E-01, 
    0.1049750419840204E-01, 
    0.1049750419840204E-01, 
    0.2100735514277743E-01, 
    0.2100735514277743E-01, 
    0.1140510361065565E+00, 
    0.1140510361065565E+00, 
    0.4811709451294231E-01, 
    0.4811709451294231E-01, 
    0.7994419804097108E-01, 
    0.7994419804097108E-01, 
    0.1010005451633049E+00, 
    0.1010005451633049E+00, 
    0.1204195881877324E-01, 
    0.1204195881877324E-01, 
    0.9474459024829298E-02, 
    0.9474459024829298E-02, 
    0.1005514424347678E+00, 
    0.1005514424347678E+00, 
    0.3161642787539286E-01, 
    0.3161642787539286E-01, 
    0.3963833050663004E-01, 
    0.3963833050663004E-01, 
    0.7350586661049985E-01, 
    0.7350586661049985E-01, 
    0.4319417861510279E-01, 
    0.4319417861510279E-01, 
    0.8810251098693814E-01, 
    0.8810251098693814E-01, 
    0.6864316028539075E-01, 
    0.6864316028539075E-01, 
    0.8257746135731812E-01, 
    0.8257746135731812E-01, 
    0.5439632620644287E-01, 
    0.5439632620644287E-01, 
    0.4386704732153978E-01, 
    0.4386704732153978E-01, 
    0.8808225769982879E-01, 
    0.8808225769982879E-01, 
    0.1534893259270625E-01, 
    0.1534893259270625E-01, 
    0.1234624197629746E+00 };

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
//    02 July 2014
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
  double xs[55] = { 
   -0.7710386602263628E+00, 
    0.7710386602263630E+00, 
    0.9803457456469387E+00, 
   -0.9803457456469384E+00, 
   -0.2292639639675523E+00, 
    0.2292639639675524E+00, 
    0.4847176019505991E-03, 
   -0.4847176019504780E-03, 
   -0.6189416389750175E+00, 
    0.6189416389750177E+00, 
    0.9587315519802511E+00, 
   -0.9587315519802511E+00, 
    0.8409306922533750E+00, 
   -0.8409306922533748E+00, 
   -0.4308042054877432E+00, 
    0.4308042054877433E+00, 
    0.4761431266211590E+00, 
   -0.4761431266211589E+00, 
    0.8651144531733139E+00, 
   -0.8651144531733139E+00, 
    0.9846617345267017E+00, 
   -0.9846617345267017E+00, 
   -0.7981639404863030E+00, 
    0.7981639404863030E+00, 
    0.6877591943414725E+00, 
   -0.6877591943414725E+00, 
   -0.3038305486106544E+00, 
    0.3038305486106544E+00, 
    0.9852576255116258E+00, 
   -0.9852576255116258E+00, 
    0.9853756930046446E+00, 
   -0.9853756930046446E+00, 
    0.7024672194580522E+00, 
   -0.7024672194580522E+00, 
    0.4589513024499020E+00, 
   -0.4589513024499019E+00, 
   -0.5838938372432102E+00, 
    0.5838938372432102E+00, 
    0.4855363777625971E+00, 
   -0.4855363777625971E+00, 
    0.1909552287968119E+00, 
   -0.1909552287968118E+00, 
    0.1970910744873101E+00, 
   -0.1970910744873101E+00, 
    0.9070140000742543E+00, 
   -0.9070140000742543E+00, 
   -0.9370706813548184E+00, 
    0.9370706813548186E+00, 
   -0.1024098809482286E+00, 
    0.1024098809482287E+00, 
    0.9018657853563646E+00, 
   -0.9018657853563646E+00, 
    0.7422255782894629E+00, 
   -0.7422255782894629E+00, 
   -0.1975779250586182E-19 };
  double ys[55] = { 
   -0.9187170657318696E+00, 
    0.9187170657318696E+00, 
   -0.9679135253250817E+00, 
    0.9679135253250819E+00, 
   -0.9437800394025085E+00, 
    0.9437800394025085E+00, 
   -0.9886578344699537E+00, 
    0.9886578344699537E+00, 
   -0.9803491213417113E+00, 
    0.9803491213417113E+00, 
   -0.8226737868824753E+00, 
    0.8226737868824755E+00, 
   -0.9649601466712245E+00, 
    0.9649601466712245E+00, 
   -0.8370492275539414E+00, 
    0.8370492275539414E+00, 
   -0.9716943047473653E+00, 
    0.9716943047473653E+00, 
   -0.6326447362896030E+00, 
    0.6326447362896030E+00, 
    0.2029425559112923E+00, 
   -0.2029425559112922E+00, 
   -0.7906135688735062E+00, 
    0.7906135688735062E+00, 
   -0.8442560578129694E+00, 
    0.8442560578129694E+00, 
   -0.3117615836793495E+00, 
    0.3117615836793495E+00, 
    0.7701659795648228E+00, 
   -0.7701659795648226E+00, 
   -0.4379432170880169E+00, 
    0.4379432170880170E+00, 
   -0.3820619012323893E+00, 
    0.3820619012323894E+00, 
   -0.6514286057161101E+00, 
    0.6514286057161101E+00, 
   -0.5711068454496305E+00, 
    0.5711068454496305E+00, 
   -0.8072896746317025E-01, 
    0.8072896746317031E-01, 
   -0.8630149364726712E+00, 
    0.8630149364726712E+00, 
   -0.3872678175415290E+00, 
    0.3872678175415290E+00, 
    0.5103334842355030E+00, 
   -0.5103334842355027E+00, 
   -0.9584329986119476E+00, 
    0.9584329986119474E+00, 
   -0.6619201369182062E+00, 
    0.6619201369182062E+00, 
   -0.1238115372273944E+00, 
    0.1238115372273945E+00, 
    0.2071876599633523E+00, 
   -0.2071876599633522E+00, 
    0.5346688849930886E-20 };
  double ws[55] = { 
    0.1261638293838951E-01, 
    0.1261638293838951E-01, 
    0.3408339905429266E-02, 
    0.3408339905429266E-02, 
    0.2796862081921473E-01, 
    0.2796862081921473E-01, 
    0.1252812914329644E-01, 
    0.1252812914329644E-01, 
    0.1635296523785200E-01, 
    0.1635296523785200E-01, 
    0.1720881227455075E-01, 
    0.1720881227455075E-01, 
    0.1523407270818440E-01, 
    0.1523407270818440E-01, 
    0.5600796522816800E-01, 
    0.5600796522816800E-01, 
    0.2382823797668716E-01, 
    0.2382823797668716E-01, 
    0.4513279974663867E-01, 
    0.4513279974663867E-01, 
    0.1931215256841166E-01, 
    0.1931215256841166E-01, 
    0.4158804216001467E-01, 
    0.4158804216001467E-01, 
    0.4685849665862760E-01, 
    0.4685849665862760E-01, 
    0.1200522449400290E+00, 
    0.1200522449400290E+00, 
    0.1238565802221201E-01, 
    0.1238565802221201E-01, 
    0.1760077392303538E-01, 
    0.1760077392303538E-01, 
    0.8264937698824523E-01, 
    0.8264937698824523E-01, 
    0.8629133710270168E-01, 
    0.8629133710270168E-01, 
    0.8660536182880048E-01, 
    0.8660536182880048E-01, 
    0.1134857467272575E+00, 
    0.1134857467272575E+00, 
    0.6518861145910534E-01, 
    0.6518861145910534E-01, 
    0.1184802238173896E+00, 
    0.1184802238173896E+00, 
    0.4767526390300979E-01, 
    0.4767526390300979E-01, 
    0.1203076112968188E-01, 
    0.1203076112968188E-01, 
    0.1010849820160845E+00, 
    0.1010849820160845E+00, 
    0.5753445241741756E-01, 
    0.5753445241741756E-01, 
    0.8946701652955226E-01, 
    0.8946701652955226E-01, 
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
//    02 July 2014
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
  double xs[68] = { 
   -0.9669786385710223E+00, 
    0.9737001842077581E+00, 
    0.9669786385710225E+00, 
   -0.9737001842077578E+00, 
   -0.2156318842512505E+00, 
    0.9910931195695962E+00, 
    0.2156318842512506E+00, 
   -0.9910931195695962E+00, 
   -0.7389660590011030E+00, 
    0.9797385272966153E+00, 
    0.7389660590011032E+00, 
   -0.9797385272966153E+00, 
    0.7689094060317012E+00, 
    0.9882749272572955E+00, 
   -0.7689094060317010E+00, 
   -0.9882749272572955E+00, 
   -0.8922402234413791E+00, 
    0.8925564983087213E+00, 
    0.8922402234413791E+00, 
   -0.8925564983087213E+00, 
    0.2617471442719549E+00, 
    0.9844702542794935E+00, 
   -0.2617471442719548E+00, 
   -0.9844702542794935E+00, 
   -0.7742833119206508E+00, 
    0.7411227454690407E+00, 
    0.7742833119206508E+00, 
   -0.7411227454690407E+00, 
   -0.5506736485553229E+00, 
    0.8796491853095826E+00, 
    0.5506736485553229E+00, 
   -0.8796491853095826E+00, 
   -0.5792562772184127E+00, 
    0.5652337954199163E+00, 
    0.5792562772184127E+00, 
   -0.5652337954199163E+00, 
   -0.1014796206724937E-01, 
    0.9024857168797702E+00, 
    0.1014796206724948E-01, 
   -0.9024857168797702E+00, 
    0.5420066475220151E+00, 
    0.9210890053684702E+00, 
   -0.5420066475220149E+00, 
   -0.9210890053684702E+00, 
    0.2943587054075071E+00, 
    0.7683262972049428E+00, 
   -0.2943587054075070E+00, 
   -0.7683262972049428E+00, 
   -0.3513695172888806E+00, 
    0.3692629613410464E+00, 
    0.3513695172888806E+00, 
   -0.3692629613410464E+00, 
   -0.3707443881794703E+00, 
    0.9667097045148131E+00, 
    0.3707443881794704E+00, 
   -0.9667097045148131E+00, 
   -0.2686897439986438E+00, 
    0.7370294813846769E+00, 
    0.2686897439986439E+00, 
   -0.7370294813846769E+00, 
   -0.1140106895094741E+00, 
    0.1969733705383891E+00, 
    0.1140106895094742E+00, 
   -0.1969733705383891E+00, 
    0.3612358695381902E-01, 
    0.5430113079937613E+00, 
   -0.3612358695381895E-01, 
   -0.5430113079937613E+00 };
  double ys[68] = { 
   -0.9737001842077582E+00, 
   -0.9669786385710224E+00, 
    0.9737001842077579E+00, 
    0.9669786385710226E+00, 
   -0.9910931195695962E+00, 
   -0.2156318842512506E+00, 
    0.9910931195695962E+00, 
    0.2156318842512507E+00, 
   -0.9797385272966153E+00, 
   -0.7389660590011031E+00, 
    0.9797385272966153E+00, 
    0.7389660590011033E+00, 
   -0.9882749272572955E+00, 
    0.7689094060317011E+00, 
    0.9882749272572955E+00, 
   -0.7689094060317009E+00, 
   -0.8925564983087213E+00, 
   -0.8922402234413791E+00, 
    0.8925564983087213E+00, 
    0.8922402234413791E+00, 
   -0.9844702542794935E+00, 
    0.2617471442719548E+00, 
    0.9844702542794935E+00, 
   -0.2617471442719547E+00, 
   -0.7411227454690407E+00, 
   -0.7742833119206508E+00, 
    0.7411227454690407E+00, 
    0.7742833119206508E+00, 
   -0.8796491853095826E+00, 
   -0.5506736485553229E+00, 
    0.8796491853095826E+00, 
    0.5506736485553229E+00, 
   -0.5652337954199163E+00, 
   -0.5792562772184127E+00, 
    0.5652337954199163E+00, 
    0.5792562772184127E+00, 
   -0.9024857168797702E+00, 
   -0.1014796206724942E-01, 
    0.9024857168797702E+00, 
    0.1014796206724953E-01, 
   -0.9210890053684702E+00, 
    0.5420066475220150E+00, 
    0.9210890053684702E+00, 
   -0.5420066475220148E+00, 
   -0.7683262972049428E+00, 
    0.2943587054075071E+00, 
    0.7683262972049428E+00, 
   -0.2943587054075070E+00, 
   -0.3692629613410464E+00, 
   -0.3513695172888806E+00, 
    0.3692629613410464E+00, 
    0.3513695172888806E+00, 
   -0.9667097045148131E+00, 
   -0.3707443881794704E+00, 
    0.9667097045148131E+00, 
    0.3707443881794705E+00, 
   -0.7370294813846769E+00, 
   -0.2686897439986438E+00, 
    0.7370294813846769E+00, 
    0.2686897439986439E+00, 
   -0.1969733705383891E+00, 
   -0.1140106895094741E+00, 
    0.1969733705383891E+00, 
    0.1140106895094742E+00, 
   -0.5430113079937613E+00, 
    0.3612358695381898E-01, 
    0.5430113079937613E+00, 
   -0.3612358695381891E-01 };
  double ws[68] = { 
    0.4697922862445027E-02, 
    0.4697922862445027E-02, 
    0.4697922862445027E-02, 
    0.4697922862445027E-02, 
    0.7136263254607511E-02, 
    0.7136263254607511E-02, 
    0.7136263254607511E-02, 
    0.7136263254607511E-02, 
    0.1293239065568239E-01, 
    0.1293239065568239E-01, 
    0.1293239065568239E-01, 
    0.1293239065568239E-01, 
    0.9398347568166604E-02, 
    0.9398347568166604E-02, 
    0.9398347568166604E-02, 
    0.9398347568166604E-02, 
    0.1884626577476044E-01, 
    0.1884626577476044E-01, 
    0.1884626577476044E-01, 
    0.1884626577476044E-01, 
    0.1572887987347023E-01, 
    0.1572887987347023E-01, 
    0.1572887987347023E-01, 
    0.1572887987347023E-01, 
    0.4107161379502558E-01, 
    0.4107161379502558E-01, 
    0.4107161379502558E-01, 
    0.4107161379502558E-01, 
    0.4035221395931435E-01, 
    0.4035221395931435E-01, 
    0.4035221395931435E-01, 
    0.4035221395931435E-01, 
    0.6647952625116643E-01, 
    0.6647952625116643E-01, 
    0.6647952625116643E-01, 
    0.6647952625116643E-01, 
    0.4719480581715914E-01, 
    0.4719480581715914E-01, 
    0.4719480581715914E-01, 
    0.4719480581715914E-01, 
    0.3594938959356454E-01, 
    0.3594938959356454E-01, 
    0.3594938959356454E-01, 
    0.3594938959356454E-01, 
    0.6892712069447091E-01, 
    0.6892712069447091E-01, 
    0.6892712069447091E-01, 
    0.6892712069447091E-01, 
    0.8060688072749707E-01, 
    0.8060688072749707E-01, 
    0.8060688072749707E-01, 
    0.8060688072749707E-01, 
    0.1530603725863855E-01, 
    0.1530603725863855E-01, 
    0.1530603725863855E-01, 
    0.1530603725863855E-01, 
    0.7313001882369689E-01, 
    0.7313001882369689E-01, 
    0.7313001882369689E-01, 
    0.7313001882369689E-01, 
    0.7447739831288605E-01, 
    0.7447739831288605E-01, 
    0.7447739831288605E-01, 
    0.7447739831288605E-01, 
    0.9487170596399580E-01, 
    0.9487170596399580E-01, 
    0.9487170596399580E-01, 
    0.9487170596399580E-01 };

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
//    02 July 2014
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
  double xs[68] = { 
   -0.9734386316165470E+00, 
    0.9744990929036832E+00, 
    0.9734386316165472E+00, 
   -0.9744990929036830E+00, 
   -0.3841574585766744E+00, 
    0.9670641778942685E+00, 
    0.3841574585766745E+00, 
   -0.9670641778942685E+00, 
    0.2986734938364671E+00, 
    0.9905525689050123E+00, 
   -0.2986734938364670E+00, 
   -0.9905525689050123E+00, 
   -0.7396581737067777E+00, 
    0.9869464369033261E+00, 
    0.7396581737067779E+00, 
   -0.9869464369033261E+00, 
   -0.1425244970455050E+00, 
    0.9733021904515969E+00, 
    0.1425244970455051E+00, 
   -0.9733021904515969E+00, 
    0.7650240374639232E+00, 
    0.9804863471920530E+00, 
   -0.7650240374639230E+00, 
   -0.9804863471920530E+00, 
   -0.7599006633708002E+00, 
    0.7279453517455540E+00, 
    0.7599006633708002E+00, 
   -0.7279453517455540E+00, 
   -0.1192987760526789E+00, 
   -0.2637912058730560E-02, 
    0.1192987760526789E+00, 
    0.2637912058730575E-02, 
   -0.8850504442537889E+00, 
    0.9022234232868145E+00, 
    0.8850504442537889E+00, 
   -0.9022234232868145E+00, 
    0.5304174652462883E+00, 
    0.9125489607085608E+00, 
   -0.5304174652462881E+00, 
   -0.9125489607085608E+00, 
   -0.2858528945041368E+00, 
    0.2941600854694212E+00, 
    0.2858528945041368E+00, 
   -0.2941600854694212E+00, 
   -0.5671850101113227E+00, 
    0.8836081660895880E+00, 
    0.5671850101113227E+00, 
   -0.8836081660895880E+00, 
    0.3174295148500719E+00, 
    0.7293427112089215E+00, 
   -0.3174295148500718E+00, 
   -0.7293427112089215E+00, 
   -0.2492430513869149E+00, 
    0.7672563284436533E+00, 
    0.2492430513869150E+00, 
   -0.7672563284436533E+00, 
   -0.5087793568494521E+00, 
    0.5623738439118215E+00, 
    0.5087793568494521E+00, 
   -0.5623738439118215E+00, 
    0.7335719396414396E-01, 
    0.8930925855397183E+00, 
   -0.7335719396414385E-01, 
   -0.8930925855397183E+00, 
    0.8350775723842838E-02, 
    0.5392457387102469E+00, 
   -0.8350775723842772E-02, 
   -0.5392457387102469E+00 };
  double ys[68] = { 
   -0.9744990929036833E+00, 
   -0.9734386316165471E+00, 
    0.9744990929036831E+00, 
    0.9734386316165473E+00, 
   -0.9670641778942685E+00, 
   -0.3841574585766744E+00, 
    0.9670641778942685E+00, 
    0.3841574585766745E+00, 
   -0.9905525689050123E+00, 
    0.2986734938364670E+00, 
    0.9905525689050123E+00, 
   -0.2986734938364669E+00, 
   -0.9869464369033261E+00, 
   -0.7396581737067778E+00, 
    0.9869464369033261E+00, 
    0.7396581737067780E+00, 
   -0.9733021904515969E+00, 
   -0.1425244970455050E+00, 
    0.9733021904515969E+00, 
    0.1425244970455051E+00, 
   -0.9804863471920530E+00, 
    0.7650240374639231E+00, 
    0.9804863471920530E+00, 
   -0.7650240374639229E+00, 
   -0.7279453517455540E+00, 
   -0.7599006633708002E+00, 
    0.7279453517455540E+00, 
    0.7599006633708002E+00, 
    0.2637912058730553E-02, 
   -0.1192987760526789E+00, 
   -0.2637912058730568E-02, 
    0.1192987760526789E+00, 
   -0.9022234232868145E+00, 
   -0.8850504442537889E+00, 
    0.9022234232868145E+00, 
    0.8850504442537889E+00, 
   -0.9125489607085608E+00, 
    0.5304174652462882E+00, 
    0.9125489607085608E+00, 
   -0.5304174652462880E+00, 
   -0.2941600854694212E+00, 
   -0.2858528945041368E+00, 
    0.2941600854694212E+00, 
    0.2858528945041368E+00, 
   -0.8836081660895880E+00, 
   -0.5671850101113227E+00, 
    0.8836081660895880E+00, 
    0.5671850101113227E+00, 
   -0.7293427112089215E+00, 
    0.3174295148500719E+00, 
    0.7293427112089215E+00, 
   -0.3174295148500718E+00, 
   -0.7672563284436533E+00, 
   -0.2492430513869149E+00, 
    0.7672563284436533E+00, 
    0.2492430513869150E+00, 
   -0.5623738439118215E+00, 
   -0.5087793568494521E+00, 
    0.5623738439118215E+00, 
    0.5087793568494521E+00, 
   -0.8930925855397183E+00, 
    0.7335719396414390E-01, 
    0.8930925855397183E+00, 
   -0.7335719396414379E-01, 
   -0.5392457387102469E+00, 
    0.8350775723842805E-02, 
    0.5392457387102469E+00, 
   -0.8350775723842739E-02 };
  double ws[68] = { 
    0.4076118519980060E-02, 
    0.4076118519980060E-02, 
    0.4076118519980060E-02, 
    0.4076118519980060E-02, 
    0.1627326938099484E-01, 
    0.1627326938099484E-01, 
    0.1627326938099484E-01, 
    0.1627326938099484E-01, 
    0.1254157952509427E-01, 
    0.1254157952509427E-01, 
    0.1254157952509427E-01, 
    0.1254157952509427E-01, 
    0.1028929333936017E-01, 
    0.1028929333936017E-01, 
    0.1028929333936017E-01, 
    0.1028929333936017E-01, 
    0.1475928282295525E-01, 
    0.1475928282295525E-01, 
    0.1475928282295525E-01, 
    0.1475928282295525E-01, 
    0.1207323692393111E-01, 
    0.1207323692393111E-01, 
    0.1207323692393111E-01, 
    0.1207323692393111E-01, 
    0.4619184040692218E-01, 
    0.4619184040692218E-01, 
    0.4619184040692218E-01, 
    0.4619184040692218E-01, 
    0.3696173437828049E-01, 
    0.3696173437828049E-01, 
    0.3696173437828049E-01, 
    0.3696173437828049E-01, 
    0.2018069481193246E-01, 
    0.2018069481193246E-01, 
    0.2018069481193246E-01, 
    0.2018069481193246E-01, 
    0.3738944032940469E-01, 
    0.3738944032940469E-01, 
    0.3738944032940469E-01, 
    0.3738944032940469E-01, 
    0.9821701539315209E-01, 
    0.9821701539315209E-01, 
    0.9821701539315209E-01, 
    0.9821701539315209E-01, 
    0.3844110871724747E-01, 
    0.3844110871724747E-01, 
    0.3844110871724747E-01, 
    0.3844110871724747E-01, 
    0.7127049386881731E-01, 
    0.7127049386881731E-01, 
    0.7127049386881731E-01, 
    0.7127049386881731E-01, 
    0.6966749913838975E-01, 
    0.6966749913838975E-01, 
    0.6966749913838975E-01, 
    0.6966749913838975E-01, 
    0.7715964130310782E-01, 
    0.7715964130310782E-01, 
    0.7715964130310782E-01, 
    0.7715964130310782E-01, 
    0.4598470092336809E-01, 
    0.4598470092336809E-01, 
    0.4598470092336809E-01, 
    0.4598470092336809E-01, 
    0.9562983140360957E-01, 
    0.9562983140360957E-01, 
    0.9562983140360957E-01, 
    0.9562983140360957E-01 };

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
//    02 July 2014
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
  double xs[81] = { 
   -0.9795110740034025E+00, 
    0.9831906073122737E+00, 
    0.9795110740034028E+00, 
   -0.9831906073122735E+00, 
   -0.7431761069248197E+00, 
    0.9923743096061538E+00, 
    0.7431761069248199E+00, 
   -0.9923743096061538E+00, 
   -0.4283144128355606E+00, 
    0.9641460474769801E+00, 
    0.4283144128355607E+00, 
   -0.9641460474769801E+00, 
    0.2195391124808899E+00, 
    0.9631697483532271E+00, 
   -0.2195391124808898E+00, 
   -0.9631697483532271E+00, 
    0.6056140907858303E+00, 
    0.9331619907848750E+00, 
   -0.6056140907858301E+00, 
   -0.9331619907848750E+00, 
    0.4538625783394974E+00, 
    0.9980174969022684E+00, 
   -0.4538625783394973E+00, 
   -0.9980174969022684E+00, 
   -0.8095537467004988E+00, 
    0.7623591488515665E+00, 
    0.8095537467004988E+00, 
   -0.7623591488515665E+00, 
   -0.1187579119827596E+00, 
    0.9879801664420653E+00, 
    0.1187579119827597E+00, 
   -0.9879801664420653E+00, 
   -0.8923235157505165E+00, 
    0.9333621871500086E+00, 
    0.8923235157505167E+00, 
   -0.9333621871500086E+00, 
    0.8231553038658227E+00, 
    0.9792360167943942E+00, 
   -0.8231553038658225E+00, 
   -0.9792360167943942E+00, 
   -0.2288711050959638E+00, 
    0.8448136056975591E+00, 
    0.2288711050959640E+00, 
   -0.8448136056975591E+00, 
   -0.6414644180013116E+00, 
    0.8887383480333905E+00, 
    0.6414644180013116E+00, 
   -0.8887383480333905E+00, 
    0.2100118285690190E-01, 
    0.9154636292013463E+00, 
   -0.2100118285690179E-01, 
   -0.9154636292013463E+00, 
    0.2939039049089534E+00, 
    0.4700673563865673E+00, 
   -0.2939039049089532E+00, 
   -0.4700673563865673E+00, 
   -0.4701209495753256E+00, 
    0.7110849452816542E+00, 
    0.4701209495753257E+00, 
   -0.7110849452816542E+00, 
   -0.2561845423520469E+00, 
    0.1372468757285573E-01, 
    0.2561845423520469E+00, 
   -0.1372468757285570E-01, 
    0.5331634078426070E+00, 
    0.6746722584255035E+00, 
   -0.5331634078426070E+00, 
   -0.6746722584255035E+00, 
    0.3458330575650539E+00, 
    0.8408056203362516E+00, 
   -0.3458330575650538E+00, 
   -0.8408056203362516E+00, 
    0.6630732857737233E-01, 
    0.6973527543224615E+00, 
   -0.6630732857737225E-01, 
   -0.6973527543224615E+00, 
   -0.2157929992274237E+00, 
    0.5168584327986239E+00, 
    0.2157929992274237E+00, 
   -0.5168584327986239E+00, 
   -0.1195405968452537E-31 };
  double ys[81] = { 
   -0.9831906073122738E+00, 
   -0.9795110740034026E+00, 
    0.9831906073122736E+00, 
    0.9795110740034029E+00, 
   -0.9923743096061538E+00, 
   -0.7431761069248198E+00, 
    0.9923743096061538E+00, 
    0.7431761069248201E+00, 
   -0.9641460474769801E+00, 
   -0.4283144128355607E+00, 
    0.9641460474769801E+00, 
    0.4283144128355608E+00, 
   -0.9631697483532271E+00, 
    0.2195391124808899E+00, 
    0.9631697483532271E+00, 
   -0.2195391124808898E+00, 
   -0.9331619907848750E+00, 
    0.6056140907858302E+00, 
    0.9331619907848750E+00, 
   -0.6056140907858300E+00, 
   -0.9980174969022684E+00, 
    0.4538625783394974E+00, 
    0.9980174969022684E+00, 
   -0.4538625783394973E+00, 
   -0.7623591488515665E+00, 
   -0.8095537467004988E+00, 
    0.7623591488515665E+00, 
    0.8095537467004988E+00, 
   -0.9879801664420653E+00, 
   -0.1187579119827596E+00, 
    0.9879801664420653E+00, 
    0.1187579119827597E+00, 
   -0.9333621871500086E+00, 
   -0.8923235157505166E+00, 
    0.9333621871500086E+00, 
    0.8923235157505168E+00, 
   -0.9792360167943942E+00, 
    0.8231553038658226E+00, 
    0.9792360167943942E+00, 
   -0.8231553038658224E+00, 
   -0.8448136056975591E+00, 
   -0.2288711050959639E+00, 
    0.8448136056975591E+00, 
    0.2288711050959640E+00, 
   -0.8887383480333905E+00, 
   -0.6414644180013116E+00, 
    0.8887383480333905E+00, 
    0.6414644180013116E+00, 
   -0.9154636292013463E+00, 
    0.2100118285690184E-01, 
    0.9154636292013463E+00, 
   -0.2100118285690173E-01, 
   -0.4700673563865673E+00, 
    0.2939039049089533E+00, 
    0.4700673563865673E+00, 
   -0.2939039049089532E+00, 
   -0.7110849452816542E+00, 
   -0.4701209495753256E+00, 
    0.7110849452816542E+00, 
    0.4701209495753257E+00, 
   -0.1372468757285574E-01, 
   -0.2561845423520469E+00, 
    0.1372468757285571E-01, 
    0.2561845423520469E+00, 
   -0.6746722584255035E+00, 
    0.5331634078426070E+00, 
    0.6746722584255035E+00, 
   -0.5331634078426070E+00, 
   -0.8408056203362516E+00, 
    0.3458330575650538E+00, 
    0.8408056203362516E+00, 
   -0.3458330575650537E+00, 
   -0.6973527543224615E+00, 
    0.6630732857737229E-01, 
    0.6973527543224615E+00, 
   -0.6630732857737220E-01, 
   -0.5168584327986239E+00, 
   -0.2157929992274237E+00, 
    0.5168584327986239E+00, 
    0.2157929992274238E+00, 
    0.3240416764471269E-32 };
  double ws[81] = { 
    0.2403280128020245E-02, 
    0.2403280128020245E-02, 
    0.2403280128020245E-02, 
    0.2403280128020245E-02, 
    0.6918304937846545E-02, 
    0.6918304937846545E-02, 
    0.6918304937846545E-02, 
    0.6918304937846545E-02, 
    0.1998132824455828E-01, 
    0.1998132824455828E-01, 
    0.1998132824455828E-01, 
    0.1998132824455828E-01, 
    0.1612406542082527E-01, 
    0.1612406542082527E-01, 
    0.1612406542082527E-01, 
    0.1612406542082527E-01, 
    0.2451719014395468E-01, 
    0.2451719014395468E-01, 
    0.2451719014395468E-01, 
    0.2451719014395468E-01, 
    0.5618083393401648E-02, 
    0.5618083393401648E-02, 
    0.5618083393401648E-02, 
    0.5618083393401648E-02, 
    0.3267989661107104E-01, 
    0.3267989661107104E-01, 
    0.3267989661107104E-01, 
    0.3267989661107104E-01, 
    0.9643554633385169E-02, 
    0.9643554633385169E-02, 
    0.9643554633385169E-02, 
    0.9643554633385169E-02, 
    0.1438022637487432E-01, 
    0.1438022637487432E-01, 
    0.1438022637487432E-01, 
    0.1438022637487432E-01, 
    0.9462403050575163E-02, 
    0.9462403050575163E-02, 
    0.9462403050575163E-02, 
    0.9462403050575163E-02, 
    0.4414700234043260E-01, 
    0.4414700234043260E-01, 
    0.4414700234043260E-01, 
    0.4414700234043260E-01, 
    0.2997776103642255E-01, 
    0.2997776103642255E-01, 
    0.2997776103642255E-01, 
    0.2997776103642255E-01, 
    0.2217921802120890E-01, 
    0.2217921802120890E-01, 
    0.2217921802120890E-01, 
    0.2217921802120890E-01, 
    0.7979169324002153E-01, 
    0.7979169324002153E-01, 
    0.7979169324002153E-01, 
    0.7979169324002153E-01, 
    0.5450753416951606E-01, 
    0.5450753416951606E-01, 
    0.5450753416951606E-01, 
    0.5450753416951606E-01, 
    0.9164051342923195E-01, 
    0.9164051342923195E-01, 
    0.9164051342923195E-01, 
    0.9164051342923195E-01, 
    0.5417703706712328E-01, 
    0.5417703706712328E-01, 
    0.5417703706712328E-01, 
    0.5417703706712328E-01, 
    0.4265496337854927E-01, 
    0.4265496337854927E-01, 
    0.4265496337854927E-01, 
    0.4265496337854927E-01, 
    0.6713307669025259E-01, 
    0.6713307669025259E-01, 
    0.6713307669025259E-01, 
    0.6713307669025259E-01, 
    0.7903551107191877E-01, 
    0.7903551107191877E-01, 
    0.7903551107191877E-01, 
    0.7903551107191877E-01, 
    0.5365512134302086E-03 };

  r8mat_row_copy ( 2, n, 0, xs, x );
  r8mat_row_copy ( 2, n, 1, ys, x );
  r8vec_copy ( n, ws, w );

  return;
}
//****************************************************************************80

void square_symq ( int degree, int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_SYMQ returns a symmetric quadrature rule for the square.
//
//  Discussion:
//
//    This procedure returns a quadrature rule for smooth functions
//    on the unit square [-1,1]^2.
//
//    All quadratures are accurate to 15 digits
//    All weights are positive and inside the square
//
//    The nodes are symmetrically arranged.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2014
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
//    1 <= DEGREE <= 20.
//
//    Input, int N, the number of nodes.
//    This can be determined by a call to RULE_FULL_SIZE(DEGREE).
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
  else
  {
    cerr << "\n";
    cerr << "SQUARE_SYMQ - Fatal error\n";
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

void square_symq_gnuplot ( int n, double x[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    SQUARE_SYMQ_GNUPLOT: GNUPLOT plot of the symmetric square quadrature rule.
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
