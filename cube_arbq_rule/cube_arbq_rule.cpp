# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

# include "cube_arbq_rule.hpp"

//****************************************************************************80

void cube_arbq ( int degree, int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_ARBQ returns a quadrature rule for the symmetric cube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2014
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
//    1 <= DEGREE <= 15.
//
//    Input, int N, the number of nodes.
//    This can be determined by a call to CUBE_ARBQ_SIZE(DEGREE).
//
//    Output, double X[3*N], the coordinates of the nodes.
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
  else
  {
    cerr << "\n";
    cerr << "CUBE_ARBQ - Fatal error\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  w_sum = r8vec_sum ( n, w );

  for ( i = 0; i < n; i++ )
  {
    w[i] = 8.0 * w[i] / w_sum;
  }

  return;
}
//****************************************************************************80

void cube_arbq_gnuplot ( int n, double x[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_ARBQ_GNUPLOT: plot of a quadrature rule for the symmetric cube.
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
//    Input, double X[3*N], the coordinates of the nodes.
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
  vertex_unit << "-1.0, -1.0, -1.0\n";
  vertex_unit << "+1.0  -1.0  -1.0\n";;
  vertex_unit << "+1.0  +1.0 -1.0\n";
  vertex_unit << "-1.0  +1.0  -1.0\n";
  vertex_unit << "-1.0  -1.0  -1.0\n";
  vertex_unit << "\n";
  vertex_unit << "-1.0  -1.0  +1.0\n";
  vertex_unit << "+1.0  -1.0  +1.0\n";
  vertex_unit << "+1.0  +1.0  +1.0\n";
  vertex_unit << "-1.0  +1.0  +1.0\n";
  vertex_unit << "-1.0  -1.0  +1.0\n";
  vertex_unit << "\n";
  vertex_unit << "-1.0  -1.0  -1.0\n";
  vertex_unit << "-1.0  -1.0  +1.0\n";
  vertex_unit << "\n";
  vertex_unit << "+1.0  -1.0  -1.0\n";
  vertex_unit << "+1.0  -1.0  +1.0\n";
  vertex_unit << "\n";
  vertex_unit << "+1.0  +1.0  -1.0\n";
  vertex_unit << "+1.0  +1.0  +1.0\n";
  vertex_unit << "\n";
  vertex_unit << "-1.0  +1.0  -1.0\n";
  vertex_unit << "-1.0  +1.0  +1.0\n";
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
    for ( i = 0; i < 3; i++ )
    {
      node_unit << x[i+j*3] << "  ";
    }
    node_unit << "\n";
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
  command_unit << "set zlabel '<--- Z --->'\n";
  command_unit << "set title '" << header << "'\n";
  command_unit << "set grid\n";
  command_unit << "set key off\n";
  command_unit << "set view equal xyz\n";
  command_unit << "set style data lines\n";
  command_unit << "set timestamp\n";
  command_unit << "splot '" << vertex_filename << "' with lines lw 3, \\\n";
  command_unit << "      '" << node_filename << "' with points pt 7 lt 0\n";
  command_unit.close ( );

  cout << "  Created command file '" << command_filename << "'\n";

  return;
}
//****************************************************************************80

int cube_arbq_size ( int degree )

//****************************************************************************80
//
//  Purpose:
//
//    CUBE_ARBQ_SIZE returns the size of quadrature rule for a cube.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2014
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
//    1 <= DEGREE <= 15.
//
//    Output, int CUBE_ARBQ_SIZE, the number of points in the
//    corresponding rule.
//
{
  int n;
  const int n_save[15] = {
      1,   4,   6,  10,  13, 
     22,  26,  42,  50,  73, 
     84, 116, 130, 172, 190 };

  if ( degree < 1 || 15 < degree ) 
  {
    cerr << "\n";
    cerr << "CUBE_ARBQ_SIZE - Fatal error!\n";
    cerr << "  Illegal value of DEGREE.\n";
    exit ( 1 );
  }

  n = n_save[degree-1];

  return n;
}
//****************************************************************************80

double *lege3eva ( int degree, double z[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGE3EVA evaluates orthonormal polynomials in the cube.
//
//  Discussion:
//
//    The number of polynomials is
//      NPOLS = ( ( DEGREE + 1 ) * ( DEGREE + 2 ) * ( DEGREE + 3 ) ) / 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 July 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DEGREE, the maximum degree.
//
//    Input, double Z[3], the evaluation point.
//
//    Output, double LEGE3EVA[NPOLS], the polynomial values.
//
{
  double *f1;
  double *f2;
  double *f3;
  int kk;
  int m;
  int n1;
  int n2;
  int npols;
  double *pols;
  double scale;
  double t;

  f1 = llegepols1 ( degree, z[0] );
  f2 = llegepols1 ( degree, z[1] );
  f3 = llegepols1 ( degree, z[2] );

  npols = ( ( degree + 1 ) * ( degree + 2 ) * ( degree + 3 ) ) / 6;
  pols = new double[npols];

  kk = 0;
  for ( m = 0; m <= degree; m++ )
  {
    for ( n2 = 0; n2 <= m; n2++ )
    {
      for ( n1 = 0; n1 <= n2; n1++ )
      {
        pols[kk] = f1[m-n2] * f2[n2-n1] * f3[n1];
        scale = 1.0;
        t = 0.5 * ( double ) ( 1 + 2 * n1 );
        scale = scale * sqrt ( t );
        t = 0.5 * ( double ) ( 1 + 2 * n2 - 2 * n1 );
        scale = scale * sqrt ( t );
        t = 0.5 * ( double ) ( 1 + 2 * m - 2 * n2 );
        scale = scale * sqrt ( t );
        pols[kk] = pols[kk] * scale;
        kk = kk + 1;
      }
    }
  }

  delete [] f1;
  delete [] f2;
  delete [] f3;

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
//    09 July 2014
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[1] = {
       0.00000000000000000E+00 };
  double ys[1] = {
       0.00000000000000000E+00 };
  double zs[1] = {
       0.00000000000000000E+00 };
  double ws[1] = {
       0.2828427124746189E+01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[4] = {
   -0.4380038262474289E+00, -0.4908226851347194E+00, 
    0.6313544088573617E+00,  0.8010097448936463E+00 };
  double ys[4] = {
    0.6281382942978499E-01, -0.1242878410373149E+00, 
    0.8674258021608136E+00, -0.9533664350988082E+00 };
  double zs[4] = {
   -0.8444012341886235E+00,  0.6401086714464984E+00, 
    0.1317904550701903E+00, -0.9165855436522309E-01 };
  double ws[4] = {
    0.7590201299956376E+00, 0.9433497505402911E+00, 
    0.6278119491594441E+00, 0.4982452950508173E+00 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[6] = {
    -0.5588399103385127E+00, 0.5588399103385127E+00, 
    -0.6245748884711381E+00, 0.6245748884711382E+00, 
    0.5311975164163759E+00, -0.5311975164163759E+00 };
  double ys[6] = {
    -0.3347842945931215E+00, 0.3347842945931215E+00, 
    -0.3608970655525763E+00, 0.3608970655525763E+00, 
    -0.9062032945290301E+00, 0.9062032945290301E+00 };
  double zs[6] = {
    -0.8055865032240838E+00, 0.8055865032240838E+00, 
    0.5832521303475051E+00, -0.5832521303475051E+00, 
    0.8103733422256782E-02, -0.8103733422256782E-02 };
  double ws[6] = {
    0.4391890453578504E+00, 0.4391890453578504E+00, 
    0.5478113077968971E+00, 0.5478113077968971E+00, 
    0.4272132092183473E+00, 0.4272132092183473E+00 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[10] = {
    0.5378961355313995E+00, 0.9230515460764596E+00, 
    -0.3017729559458644E+00, -0.4286245373611857E-01, 
    -0.4919281030834148E+00, -0.6331783925165083E+00, 
    0.3433859697316836E+00, 0.7159140285867756E+00, 
    -0.9634379969028102E+00, 0.9533570645781855E+00 };
  double ys[10] = {
    -0.2921294864376587E+00, 0.1338187840340383E-01, 
    -0.6336712971880923E+00, 0.8637293692260690E+00, 
    0.6650083378678769E+00, -0.6579587923349518E+00, 
    0.1314461211661866E+00, -0.9334671020954506E+00, 
    0.2821274186555609E+00, 0.9271218349088852E+00 };
  double zs[10] = {
    -0.9904747768071651E+00, 0.9034506944137641E+00, 
    0.7478682467863593E+00, 0.6980492706389707E+00, 
    -0.6792172628059848E+00, -0.5142053190660802E+00, 
    -0.1491834467493042E-01, -0.1124003569988050E-01, 
    0.4370157998712708E+00, -0.4433572115706369E+00 };
  double ws[10] = {
    0.2015871863001034E+00, 0.1475782644463766E+00, 
    0.3427728700669111E+00, 0.2477511243780946E+00, 
    0.3268948471688580E+00, 0.3510375486490554E+00, 
    0.6548775481546729E+00, 0.2014557914699008E+00, 
    0.2371445377295266E+00, 0.1173274063826890E+00 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[13] = {
    -0.2648598061927048E-17, 0.8182612016952856E+00, 
    -0.8182612016952856E+00, 0.3320846059655003E+00, 
    -0.3320846059655003E+00, -0.7192063268453662E+00, 
    0.7192063268453662E+00, 0.8748312857031093E+00, 
    -0.8748312857031093E+00, 0.6703438980154780E+00, 
    -0.6703438980154780E+00, -0.1288932781094893E-01, 
    0.1288932781094893E-01 };
  double ys[13] = {
    0.2334819700920049E-17, 0.6996915596126003E+00, 
    -0.6996915596126003E+00, 0.8164106580812136E+00, 
    -0.8164106580812136E+00, 0.3823687870227997E+00, 
    -0.3823687870227997E+00, -0.4498159475180819E-01, 
    0.4498159475180820E-01, -0.8936284816367571E+00, 
    0.8936284816367571E+00, -0.6644842917167314E+00, 
    0.6644842917167314E+00 };
  double zs[13] = {
    -0.2514317405785843E-17, -0.3279435833702815E+00, 
    0.3279435833702815E+00, 0.6999000775245028E+00, 
    -0.6999000775245028E+00, 0.7766614685968243E+00, 
    -0.7766614685968243E+00, 0.7066212170288233E+00, 
    -0.7066212170288233E+00, 0.1368716985635272E+00, 
    -0.1368716985635272E+00, 0.9082737241365967E+00, 
    -0.9082737241365967E+00 };
  double ws[13] = {
    0.5954583420518295E+00, 0.1966394408155049E+00, 
    0.1966394408155048E+00, 0.1972244637798873E+00, 
    0.1972244637798874E+00, 0.2056648105428643E+00, 
    0.2056648105428642E+00, 0.1744527123656727E+00, 
    0.1744527123656728E+00, 0.1737372846337885E+00, 
    0.1737372846337885E+00, 0.1687656792094623E+00, 
    0.1687656792094623E+00 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[22] = {
    -0.6520106181245787E+00, -0.8839112733197223E+00, 
    -0.4110873304651753E+00, -0.5483810510367449E+00, 
    0.2892216600948648E+00, 0.6342733723572165E+00, 
    0.9277349401977213E+00, 0.1959145493384946E+00, 
    0.9698859394138759E+00, 0.2417507986598153E+00, 
    0.8736828548839980E+00, -0.2473250700893141E+00, 
    -0.8359900749051526E+00, -0.7010811197269031E+00, 
    0.6850729043876912E+00, 0.8706999495161931E+00, 
    0.4433962826908943E+00, -0.2932096565623538E+00, 
    -0.2633910887960743E+00, -0.9524343530691857E+00, 
    0.5584942649156998E+00, -0.9550226927399118E+00 };
  double ys[22] = {
    -0.9205233816190030E+00, 0.6254092057729300E+00, 
    -0.6525192767687262E+00, 0.6675277466837514E+00, 
    0.4977300744517685E+00, -0.8473122703810076E+00, 
    -0.2928745039351578E+00, -0.1497907581315417E+00, 
    0.3885854764755980E+00, -0.3158285918947235E+00, 
    0.8526532103620746E+00, 0.4591872919393147E+00, 
    -0.2440514778109439E+00, 0.4373390182884227E-01, 
    0.3426959807532622E+00, -0.5660924679198789E+00, 
    0.9081219217522495E+00, 0.9293150291789362E+00, 
    -0.7810388725033286E+00, 0.8260274928505457E+00, 
    -0.9767688598333689E+00, -0.8519217238424458E+00 };
  double zs[22] = {
    -0.9890506976739800E+00, 0.9406229500121331E+00, 
    0.9356142377961858E+00, -0.9192764154799976E+00, 
    0.8878969382797651E+00, -0.7870867837034067E+00, 
    0.8660830979298149E+00, -0.8024709071632448E+00, 
    -0.8880711289349531E+00, 0.4529703911468836E+00, 
    0.5290490417189219E+00, -0.1854203758534000E+00, 
    -0.5963095626351139E+00, 0.4902173455917232E+00, 
    -0.2880115314573995E-01, -0.2057637646023516E+00, 
    -0.6042328075390052E+00, 0.4966484117970202E+00, 
    -0.2095900180453476E+00, -0.2268797125337812E+00, 
    0.6071604819572268E+00, 0.3850513071757548E+00 };
  double ws[22] = {
    0.3178879561185747E-01, 0.4164076216957761E-01, 
    0.9948981399775451E-01, 0.9210852733177771E-01, 
    0.1401943913356885E+00, 0.7264966274445117E-01, 
    0.6525340756744888E-01, 0.2204173451372487E+00, 
    0.5094524937590311E-01, 0.2682685881175325E+00, 
    0.6962627523546265E-01, 0.2686201905068337E+00, 
    0.1619957777392285E+00, 0.2267482460842759E+00, 
    0.2223611442058410E+00, 0.1285943421530201E+00, 
    0.1171163684835690E+00, 0.1148693655898215E+00, 
    0.2337745730572889E+00, 0.6561672504107097E-01, 
    0.7668684384852936E-01, 0.5966072941200686E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[26] = {
    0.9181562358887896E+00, 0.4489608068872516E+00, 
    -0.3858744893115341E+00, -0.8375581084754601E+00, 
    0.8759137241702486E+00, 0.3582970649375837E+00, 
    -0.3606806413292835E-02, 0.8240785780348802E+00, 
    0.3960539405501686E+00, 0.6764493821739725E+00, 
    -0.7078657252089905E+00, -0.9910280819070909E+00, 
    -0.1644247691119722E-01, -0.9181562358887897E+00, 
    -0.4489608068872515E+00, 0.3858744893115342E+00, 
    0.8375581084754601E+00, -0.8759137241702486E+00, 
    -0.3582970649375837E+00, 0.3606806413292889E-02, 
    -0.8240785780348803E+00, -0.3960539405501686E+00, 
    -0.6764493821739724E+00, 0.7078657252089905E+00, 
    0.9910280819070909E+00, 0.1644247691119727E-01 };
  double ys[26] = {
    -0.8730646210868583E+00, 0.3003747003093036E+00, 
    -0.8537024640925601E+00, -0.5290274183292351E+00, 
    -0.2186272204366992E+00, -0.8524776744046263E+00, 
    0.4528069580583293E+00, -0.7958792991210972E+00, 
    -0.3231996926348866E+00, 0.9204855453579330E+00, 
    0.1300957145548008E+00, -0.6507820069674347E+00, 
    0.9452977797065001E+00, 0.8730646210868583E+00, 
    -0.3003747003093037E+00, 0.8537024640925601E+00, 
    0.5290274183292349E+00, 0.2186272204366992E+00, 
    0.8524776744046263E+00, -0.4528069580583293E+00, 
    0.7958792991210972E+00, 0.3231996926348866E+00, 
    -0.9204855453579330E+00, -0.1300957145548009E+00, 
    0.6507820069674347E+00, -0.9452977797065000E+00 };
  double zs[26] = {
    -0.8219732697110896E+00, 0.2702488394739114E+00, 
    -0.8405117117885644E+00, 0.8255893620180769E+00, 
    -0.2154761758370911E+00, -0.2372349012690425E+00, 
    -0.5304006751850494E+00, 0.5660607048365347E+00, 
    -0.8528418042921908E+00, -0.2586944089083810E+00, 
    -0.9317986614089844E+00, -0.5844907560422828E+00, 
    -0.9724859293584232E+00, 0.8219732697110896E+00, 
    -0.2702488394739114E+00, 0.8405117117885644E+00, 
    -0.8255893620180769E+00, 0.2154761758370910E+00, 
    0.2372349012690426E+00, 0.5304006751850494E+00, 
    -0.5660607048365347E+00, 0.8528418042921908E+00, 
    0.2586944089083810E+00, 0.9317986614089844E+00, 
    0.5844907560422828E+00, 0.9724859293584232E+00 };
  double ws[26] = {
    0.2818669485625658E-01, 0.2457778545108605E+00, 
    0.7322107549897901E-01, 0.7474403008707498E-01, 
    0.1424052690706277E+00, 0.1464998318560304E+00, 
    0.2315603483757853E+00, 0.8404353770943794E-01, 
    0.1402956703321254E+00, 0.8519832938411190E-01, 
    0.7845038486431180E-01, 0.5021162675106791E-01, 
    0.3361890907642531E-01, 0.2818669485625655E-01, 
    0.2457778545108605E+00, 0.7322107549897898E-01, 
    0.7474403008707500E-01, 0.1424052690706277E+00, 
    0.1464998318560304E+00, 0.2315603483757854E+00, 
    0.8404353770943793E-01, 0.1402956703321254E+00, 
    0.8519832938411191E-01, 0.7845038486431183E-01, 
    0.5021162675106793E-01, 0.3361890907642534E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[42] = {
    -0.4475988003404910E+00, 0.1190967237337572E+00, 
    -0.4490161785347289E+00, 0.8651765155230445E+00, 
    0.6505429233922998E+00, 0.2589210616796688E+00, 
    0.9656989388220407E+00, 0.2260702905987988E+00, 
    -0.9447161047888400E+00, -0.3175860090502710E+00, 
    -0.5229105079658517E+00, 0.7843159489865733E+00, 
    -0.8185789055564340E+00, 0.2825863387132066E+00, 
    -0.7607858417694694E+00, 0.5866226772335372E+00, 
    0.9548626078195578E+00, -0.8132498468766962E+00, 
    -0.1585020553751056E+00, 0.4773375480094558E+00, 
    0.1284049655553692E+00, 0.1291886447672361E+00, 
    -0.5532190735732102E+00, -0.3228367990078276E+00, 
    0.6814247401185878E+00, 0.7557777797079770E+00, 
    -0.9000019573243127E+00, -0.2575843081671350E+00, 
    0.9631832063187018E+00, 0.6908303030190355E+00, 
    0.1821056025197366E+00, -0.8221919687508521E+00, 
    0.8588518473343788E+00, -0.8972339906902843E+00, 
    -0.5358934078144110E+00, -0.5529992683228151E+00, 
    -0.2001426466255732E-02, 0.9378537204512910E+00, 
    -0.9715016981164607E+00, -0.9829083600782837E+00, 
    0.9958977682871862E+00, 0.5952915553901705E+00 };
  double ys[42] = {
    -0.4685908229573946E+00, -0.1796519363018103E+00, 
    0.8874428218089706E+00, 0.7669380410790798E+00, 
    -0.3322370805619765E+00, 0.4117713778755101E+00, 
    0.2677042158057942E+00, -0.6809206038797824E+00, 
    0.3428366182117646E+00, 0.9241061411826208E+00, 
    -0.9900623346092050E+00, -0.8016110436407171E+00, 
    0.4049960292973258E+00, 0.7556974668347299E+00, 
    -0.4707571060333716E+00, 0.1360876107297970E+00, 
    -0.9570347291094887E+00, -0.9999285474000510E+00, 
    -0.8417098248767442E+00, 0.7312551571808507E+00, 
    -0.8410385448275977E+00, 0.4262088590943245E+00, 
    -0.2958999490543837E+00, 0.3165179774473898E+00, 
    0.9500245471112998E+00, -0.6655395544878605E+00, 
    -0.5042046483103894E+00, -0.4401404838267679E-01, 
    0.1346519238430116E+00, -0.8848406668667183E-01, 
    -0.4368801480983787E+00, 0.8270836779598612E+00, 
    0.6421742345932883E+00, 0.1349979722216108E+00, 
    0.6936454950162161E+00, -0.8510571525270378E+00, 
    0.9631601523763739E+00, 0.9984857323994113E+00, 
    -0.8588323161863770E+00, 0.9585966990212117E+00, 
    -0.5890142881200693E+00, -0.9993743349557130E+00 };
  double zs[42] = {
    -0.9772314603244738E+00, 0.9815070306506343E+00, 
    0.9720526068033335E+00, -0.9045817066028675E+00, 
    -0.9182437310533620E+00, -0.8841947772730918E+00, 
    0.9557023227491765E+00, 0.7227916414981562E+00, 
    0.8941205664955980E+00, -0.8999986131977584E+00, 
    0.9755454910061744E+00, 0.8657944704400204E+00, 
    -0.8703196425126886E+00, -0.5456656708443869E+00, 
    0.8251281474414323E+00, -0.3076953975629710E+00, 
    -0.8869272628530901E+00, -0.9230410589688908E+00, 
    0.5313897417255299E+00, 0.8095112746222346E+00, 
    -0.7431662230334861E+00, 0.1742031738076680E+00, 
    0.2874814736221014E+00, 0.6894604248684468E+00, 
    -0.4353446530984859E+00, -0.4153364490406374E+00, 
    -0.6226993610484064E+00, -0.5868262229626008E+00, 
    -0.6014685458288908E+00, 0.5571508078000768E+00, 
    -0.2048268508291324E-01, 0.4917939745803034E+00, 
    0.1369232794659417E+00, 0.7756109558580779E-02, 
    -0.2857541181158879E+00, -0.2330139351444999E+00, 
    0.2669563650277833E+00, 0.6741006875162401E+00, 
    0.4066025204718843E+00, -0.5565937248525152E+00, 
    0.2510714338774305E+00, 0.1349294493439425E+00 };
  double ws[42] = {
    0.4086705877945525E-01, 0.5411955154283058E-01, 
    0.2383675183022706E-01, 0.2393935244990715E-01, 
    0.5813502964307203E-01, 0.7437498461386277E-01, 
    0.1984960102939174E-01, 0.6544081092072129E-01, 
    0.2659353157651851E-01, 0.3231832797003169E-01, 
    0.1304035568169426E-01, 0.3740726636480241E-01, 
    0.5743227080083509E-01, 0.7044771222211986E-01, 
    0.6713400141335041E-01, 0.1088370624350729E+00, 
    0.9430601398037914E-02, 0.9861052589481853E-02, 
    0.7173538891438543E-01, 0.8420277392308859E-01, 
    0.8298939034694322E-01, 0.1440733281281141E+00, 
    0.1295131217734937E+00, 0.1455803502566465E+00, 
    0.3601936507491939E-01, 0.8420800025743772E-01, 
    0.6488371524277974E-01, 0.1744028592454575E+00, 
    0.4595502978127373E-01, 0.1382478089544316E+00, 
    0.1725289949653584E+00, 0.6019513413223115E-01, 
    0.8271749698617868E-01, 0.9333260474256072E-01, 
    0.1346393673413504E+00, 0.9612888710317599E-01, 
    0.6100738206538034E-01, 0.1085509425473334E-01, 
    0.2871508442754113E-01, 0.1504859125563503E-01, 
    0.3740891876232848E-01, 0.4097311354933509E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[50] = {
    0.7677593774006314E+00, 
    -0.7677593774006314E+00, -0.3207821329562814E+00, 
    0.3207821329562814E+00, -0.4136209586513314E+00, 
    0.4136209586513314E+00, -0.6852061275740796E+00, 
    0.6852061275740795E+00, 0.9244405033583697E+00, 
    -0.9244405033583697E+00, -0.6411525460395265E+00, 
    0.6411525460395265E+00, -0.6524611766154725E+00, 
    0.6524611766154725E+00, -0.9637289125097334E+00, 
    0.9637289125097335E+00, -0.4042812911589966E-01, 
    0.4042812911589966E-01, 0.7010911265568170E+00, 
    -0.7010911265568170E+00, -0.4018375553483048E+00, 
    0.4018375553483048E+00, 0.2427999428831116E+00, 
    -0.2427999428831116E+00, -0.8888214064543165E+00, 
    0.8888214064543165E+00, 0.5686958127555934E+00, 
    -0.5686958127555934E+00, -0.1007305440999530E+00, 
    0.1007305440999530E+00, 0.9627142008988805E+00, 
    -0.9627142008988805E+00, 0.5575105029618763E+00, 
    -0.5575105029618763E+00, -0.2006401852932052E+00, 
    0.2006401852932052E+00, 0.1276245748755967E+00, 
    -0.1276245748755967E+00, 0.5324626645500558E+00, 
    -0.5324626645500558E+00, -0.8230657430079429E+00, 
    0.8230657430079429E+00, -0.9171428680981173E+00, 
    0.9171428680981173E+00, 0.9753289529764423E+00, 
    -0.9753289529764423E+00, 0.7278991004245323E+00, 
    -0.7278991004245322E+00, 0.9084671271535661E+00, 
    -0.9084671271535661E+00 };
  double ys[50] = {
    -0.5215705733515856E-02, 
    0.5215705733515811E-02, -0.9255213178288733E+00, 
    0.9255213178288733E+00, 0.3254858593442050E+00, 
    -0.3254858593442050E+00, -0.8673453037549068E+00, 
    0.8673453037549068E+00, 0.1651688473834196E+00, 
    -0.1651688473834196E+00, 0.9161123256468909E+00, 
    -0.9161123256468909E+00, 0.8789047081433894E+00, 
    -0.8789047081433894E+00, 0.7092728961533591E+00, 
    -0.7092728961533591E+00, 0.5224015774226531E+00, 
    -0.5224015774226531E+00, 0.4986074979684522E+00, 
    -0.4986074979684522E+00, -0.5887106445666494E-01, 
    0.5887106445666495E-01, 0.4841822578088601E+00, 
    -0.4841822578088601E+00, -0.7052476161004777E+00, 
    0.7052476161004776E+00, -0.7991952932799359E-01, 
    0.7991952932799361E-01, 0.2884264730944422E+00, 
    -0.2884264730944421E+00, -0.5964266662509132E+00, 
    0.5964266662509132E+00, -0.7120073930331048E+00, 
    0.7120073930331048E+00, 0.8586009498349154E+00, 
    -0.8586009498349154E+00, 0.9369688457657286E+00, 
    -0.9369688457657286E+00, 0.6794094006908223E+00, 
    -0.6794094006908225E+00, 0.4202573751253162E+00, 
    -0.4202573751253164E+00, -0.1970879922320003E+00, 
    0.1970879922320003E+00, 0.8523907907745764E+00, 
    -0.8523907907745764E+00, 0.9938423815326598E+00, 
    -0.9938423815326598E+00, -0.9848158730090135E+00, 
    0.9848158730090135E+00 };
  double zs[50] = {
    -0.9944420260442561E+00, 
    0.9944420260442561E+00, 0.9672106847608946E+00, 
    -0.9672106847608946E+00, -0.9450653904792801E+00, 
    0.9450653904792801E+00, -0.9208144764213119E+00, 
    0.9208144764213120E+00, 0.8857214694746194E+00, 
    -0.8857214694746194E+00, 0.8845112554423256E+00, 
    -0.8845112554423256E+00, -0.8251225389279271E+00, 
    0.8251225389279271E+00, -0.8150208317048079E+00, 
    0.8150208317048079E+00, 0.8824401782407778E+00, 
    -0.8824401782407778E+00, 0.5971438826916258E+00, 
    -0.5971438826916258E+00, -0.5274658210290475E+00, 
    0.5274658210290475E+00, -0.7239973446191893E+00, 
    0.7239973446191893E+00, 0.7834569443713458E+00, 
    -0.7834569443713458E+00, -0.6010889720637044E+00, 
    0.6010889720637044E+00, 0.5589543903569408E-01, 
    -0.5589543903569408E-01, -0.7052432618789399E+00, 
    0.7052432618789399E+00, -0.2951319795385468E+00, 
    0.2951319795385468E+00, -0.4212601390804162E+00, 
    0.4212601390804162E+00, 0.5022536274483731E+00, 
    -0.5022536274483730E+00, -0.3981899691586369E-01, 
    0.3981899691586369E-01, -0.3419107086504177E+00, 
    0.3419107086504178E+00, 0.1910346742820620E+00, 
    -0.1910346742820620E+00, 0.4032363262946108E+00, 
    -0.4032363262946109E+00, -0.3785656001274115E+00, 
    0.3785656001274115E+00, 0.8086398500198494E-02, 
    -0.8086398500198474E-02 };
  double ws[50] = {
    0.2538811854621882E-01, 
    0.2538811854621883E-01, 0.1762240779978733E-01, 
    0.1762240779978734E-01, 0.4976713220360957E-01, 
    0.4976713220360960E-01, 0.2355615731022458E-01, 
    0.2355615731022460E-01, 0.2724595042792551E-01, 
    0.2724595042792553E-01, 0.2403649336578707E-01, 
    0.2403649336578707E-01, 0.2946300694577264E-01, 
    0.2946300694577266E-01, 0.1432705784656010E-01, 
    0.1432705784656010E-01, 0.7802654597787825E-01, 
    0.7802654597787825E-01, 0.6022424126025965E-01, 
    0.6022424126025969E-01, 0.1068682771411129E+00, 
    0.1068682771411130E+00, 0.9993782502170498E-01, 
    0.9993782502170505E-01, 0.3727854515440352E-01, 
    0.3727854515440351E-01, 0.9828312978563304E-01, 
    0.9828312978563307E-01, 0.1296245426718050E+00, 
    0.1296245426718050E+00, 0.3099443142302383E-01, 
    0.3099443142302385E-01, 0.9399083827155912E-01, 
    0.9399083827155913E-01, 0.8169878897651819E-01, 
    0.8169878897651825E-01, 0.5623588432426985E-01, 
    0.5623588432426987E-01, 0.1056612200118129E+00, 
    0.1056612200118130E+00, 0.9074609788069554E-01, 
    0.9074609788069565E-01, 0.6975340566577869E-01, 
    0.6975340566577869E-01, 0.2198686280707488E-01, 
    0.2198686280707488E-01, 0.2302735019253151E-01, 
    0.2302735019253151E-01, 0.1846925136114678E-01, 
    0.1846925136114678E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[73] = {
    -0.9290008391594218E+00, 0.4559324934613976E+00, 
    -0.8503323372224584E+00, -0.6258718000763578E+00, 
    -0.7657521094877171E+00, 0.5194124679448002E+00, 
    -0.1876523482951094E-01, 0.7781249769033519E+00, 
    0.7907391104275550E+00, 0.1931186493871086E-01, 
    0.7973052188120806E-01, 0.8940978313791650E+00, 
    -0.3969913823357744E+00, -0.5840731119485252E+00, 
    -0.3054403489443687E+00, -0.6250814570086429E+00, 
    -0.4273170491243774E+00, 0.8887351719726477E+00, 
    -0.1425902095867877E+00, 0.1007385849466660E+00, 
    0.9914273856893784E+00, -0.2110597230830005E+00, 
    -0.9344242484094948E+00, 0.9254796648599809E+00, 
    0.5871962886124195E+00, -0.5644486965116591E+00, 
    -0.3450006491762153E+00, -0.1925706173126046E+00, 
    0.6150095319350374E+00, -0.9410296227058319E+00, 
    0.3990765859738513E+00, 0.3675511374298872E+00, 
    0.9468259750893930E+00, -0.8323347102578033E+00, 
    0.9487564158646586E+00, 0.6918585765278729E-01, 
    0.7257633299128992E+00, -0.3310667096656268E+00, 
    0.3501561292254360E+00, -0.2984756141933171E+00, 
    0.7448197310529086E+00, 0.8904795041635960E-01, 
    -0.8524054154380974E+00, -0.9797468974234023E+00, 
    0.2069917638147098E+00, -0.9222673238089781E+00, 
    -0.7404783544230075E+00, 0.4462428682423911E+00, 
    0.2499725669820603E+00, -0.5978592395842010E+00, 
    0.3880147299726153E+00, 0.9540597896735835E+00, 
    -0.1204215373750881E+00, 0.5635716786856889E+00, 
    -0.6198158466921444E+00, 0.5423781905790268E+00, 
    0.2041052298129538E+00, 0.8615414935964518E+00, 
    0.9625611083095826E+00, -0.5668586728612423E+00, 
    -0.9565835991972781E+00, 0.9290333570079982E+00, 
    0.9603761783535766E+00, -0.9001174228174761E+00, 
    0.7409829822994354E+00, -0.9186328686279674E+00, 
    0.7988074968670636E+00, -0.8055374206467150E+00, 
    -0.2303785439930381E+00, -0.7165954822802608E+00, 
    0.7098003466268242E+00, -0.2370105920082338E+00, 
    -0.9973208321714092E+00 };
  double ys[73] = {
    0.5674482687326375E+00, -0.2911223108663519E+00, 
    0.1975635105828826E+00, -0.6882273878784136E+00, 
    -0.9313977287710445E+00, -0.5339817326366849E+00, 
    0.5593088502110057E+00, 0.6938671668589624E+00, 
    0.7457324218656869E+00, -0.4591331058869646E-01, 
    0.5440689882793791E+00, -0.8616026786491190E+00, 
    0.7739289299329076E+00, -0.2114642239345504E+00, 
    -0.3131293885912573E-01, 0.9186402247572539E+00, 
    0.9877228633852757E+00, -0.8848784715526166E+00, 
    -0.8430339428373445E+00, -0.9121169131818918E+00, 
    0.8352108151249428E-01, -0.5108468353410851E+00, 
    0.7934672273226390E+00, 0.5843397847611630E-02, 
    0.1982224584099894E+00, 0.6529491925552273E+00, 
    0.3453783893258016E+00, -0.6965527925846071E+00, 
    -0.6030971553224019E+00, -0.9773971815452341E+00, 
    0.9306828192128160E+00, 0.1729189273247773E+00, 
    0.9423175295395478E+00, 0.3117471500716507E+00, 
    0.8484099142032112E+00, -0.4549656193618034E+00, 
    0.4364950089820075E+00, -0.2119737649763384E-01, 
    -0.7501696781614142E-01, 0.3377085041383336E+00, 
    0.5119386721405460E+00, -0.2861289818382340E+00, 
    -0.9435758438429453E+00, -0.5501358135716665E+00, 
    -0.8151761180652929E+00, -0.4382794463039796E+00, 
    0.1037779030310166E+00, -0.9327537847317822E+00, 
    0.7781924335471115E+00, -0.6696667136365475E+00, 
    0.9993628918113781E+00, 0.6532557117383904E+00, 
    0.8321683233238897E+00, 0.8340145881278882E+00, 
    -0.3990216364491641E+00, -0.7026952932032947E+00, 
    0.4592299742407868E+00, -0.4282015920212867E+00, 
    0.2004721513667550E+00, -0.8934779235340296E+00, 
    0.8942502209608590E+00, -0.8738897670233189E+00, 
    -0.5705526001256189E+00, -0.1379917314855049E+00, 
    -0.1665061725494578E+00, -0.7946830985787409E+00, 
    0.9767806321337382E+00, 0.9597550703525793E+00, 
    0.9603148869900205E+00, 0.6723297616898366E+00, 
    -0.9804400598708566E+00, -0.9949665182949334E+00, 
    0.4191106171156644E+00 };
  double zs[73] = {
    -0.9933978864296309E+00, 0.9901550709677349E+00, 
    0.9670485250396404E+00, 0.9686386295527810E+00, 
    -0.9482984392355279E+00, -0.9781218402197481E+00, 
    -0.9668923135113735E+00, 0.9569807092335304E+00, 
    -0.9255203116327040E+00, -0.8549048193170219E+00, 
    0.9206982447618198E+00, 0.9211244686937629E+00, 
    -0.7317954338646777E+00, -0.9289922028324820E+00, 
    0.8753848332253249E+00, -0.8882119979065144E+00, 
    0.9878899913703325E+00, -0.8880356823987890E+00, 
    -0.8903407082095307E+00, 0.8876467178896359E+00, 
    -0.9463172167846347E+00, 0.7267804723276519E+00, 
    0.8338672039848504E+00, 0.8514657617058597E+00, 
    -0.8105375103476049E+00, 0.7553648186259019E+00, 
    -0.5850428048791508E+00, -0.1054105129348611E+00, 
    0.7438442367847414E+00, 0.8658386707444876E+00, 
    0.7640316811327564E+00, 0.6612763241373982E+00, 
    0.7191463486728814E+00, -0.7362541429366616E+00, 
    -0.6578813706843174E+00, -0.6173927832617907E+00, 
    0.5688448669683742E+00, -0.2424692431476031E+00, 
    -0.3220397059827600E+00, 0.2741662573733258E+00, 
    -0.4245710530722387E+00, 0.3779936020945108E+00, 
    -0.4890743188433689E+00, -0.8427007062812649E+00, 
    0.3110361053541443E+00, 0.7265683900043480E+00, 
    0.5114587028650807E+00, -0.6590556333897338E+00, 
    -0.5656407036165850E+00, -0.5801968526465481E+00, 
    -0.8814901582749934E+00, 0.3049921985976026E+00, 
    0.4604030935118299E+00, 0.1644483100006628E+00, 
    0.1392355684867057E+00, -0.1804588400155518E+00, 
    -0.1659257531038997E-01, -0.6408450861337968E+00, 
    -0.2283594679036603E+00, 0.5423020975283729E+00, 
    -0.5488666370681280E+00, -0.1935317969193018E+00, 
    0.4390831472700056E+00, -0.3197962461756297E+00, 
    0.1476513989767390E+00, 0.1104669057496399E+00, 
    -0.1823430458926981E+00, 0.3668414192631035E+00, 
    -0.2174815141367195E+00, -0.1332596398840131E+00, 
    0.4442289032147039E+00, -0.1887850946760386E+00, 
    0.2545601348113754E+00 };
  double ws[73] = {
    0.7213361346970690E-02, 0.2333141303565080E-01, 
    0.1679006811176461E-01, 0.1914958033825363E-01, 
    0.7581033822063699E-02, 0.2434386571403711E-01, 
    0.2909344157431100E-01, 0.1700037166331166E-01, 
    0.1673233204518952E-01, 0.3947746719767090E-01, 
    0.3908321158531827E-01, 0.1015197562376277E-01, 
    0.2887452492906294E-01, 0.4051167157357858E-01, 
    0.5023577704859834E-01, 0.1585721420277247E-01, 
    0.7185787768043105E-02, 0.1081642505043930E-01, 
    0.3035331237698869E-01, 0.2526839495335092E-01, 
    0.8807297819671642E-02, 0.5593986671876056E-01, 
    0.1373773962571757E-01, 0.2661615086380097E-01, 
    0.5999969968705832E-01, 0.4972097477068091E-01, 
    0.6999983849508175E-01, 0.5598146788784706E-01, 
    0.5112821739971896E-01, 0.4853200937844245E-02, 
    0.2643720405354587E-01, 0.7488081622291798E-01, 
    0.8050832752047557E-02, 0.4291659968333356E-01, 
    0.1294625305397554E-01, 0.7403978691481330E-01, 
    0.5117832177948896E-01, 0.7869167027903835E-01, 
    0.7987922983101066E-01, 0.7935975973606037E-01, 
    0.5190634933748844E-01, 0.8887073925421256E-01, 
    0.1520654567592403E-01, 0.1437138357007375E-01, 
    0.5539120615704855E-01, 0.3308981979630965E-01, 
    0.6553371660845619E-01, 0.3169478052156868E-01, 
    0.6276849569730045E-01, 0.6211610103098184E-01, 
    0.1053947144115475E-01, 0.2331767126681452E-01, 
    0.6027945648014443E-01, 0.4978984285341077E-01, 
    0.7974451931076965E-01, 0.6688248399941439E-01, 
    0.9777913635381968E-01, 0.5098104982036934E-01, 
    0.2951246816860325E-01, 0.4353390845153333E-01, 
    0.1466253257039164E-01, 0.2081467883348036E-01, 
    0.2633254465708832E-01, 0.5423792444480278E-01, 
    0.8745492482441419E-01, 0.3119590040598747E-01, 
    0.1563786213616101E-01, 0.2105985290616827E-01, 
    0.3523362981701141E-01, 0.7539771307788777E-01, 
    0.1972567979376486E-01, 0.2513903851381411E-01, 
    0.2400953849626596E-01 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[84] = {
    0.2499788650260321E+00, 0.3921769824433600E+00, 
    -0.7499317903601719E+00, 0.1970279946664021E+00, 
    -0.4574811189815491E+00, 0.7920994458748165E+00, 
    0.9496878040642897E+00, 0.5494991726029024E+00, 
    -0.9415698074379264E+00, -0.1549790418658323E+00, 
    0.9014007983667046E+00, -0.4493885995928424E+00, 
    -0.7095579485335214E+00, -0.8736487328036152E+00, 
    -0.3264862525811954E-01, -0.1103424065156252E+00, 
    -0.5637539669697259E+00, 0.2945721905451988E+00, 
    -0.8603062445393469E+00, -0.7507058748066501E+00, 
    -0.2780639680677569E+00, 0.7047555601238907E+00, 
    -0.6295944388807833E+00, -0.9095323899548270E+00, 
    0.2964917498074751E+00, 0.8851424078101294E+00, 
    -0.7219563370216030E+00, -0.9230518045006708E+00, 
    0.5631064783902135E+00, -0.9323802418890973E+00, 
    0.9285744095618478E+00, -0.3446313474069602E-01, 
    0.5285124795902665E+00, -0.6564153211279398E+00, 
    0.9670556758395562E+00, 0.9746827810441796E+00, 
    -0.1497145509400533E+00, 0.5126157047506474E-01, 
    0.5771162069337331E+00, -0.9882152664587444E+00, 
    0.8664710684798440E+00, 0.1554480104048762E+00, 
    -0.2499788650260319E+00, -0.3921769824433600E+00, 
    0.7499317903601719E+00, -0.1970279946664021E+00, 
    0.4574811189815492E+00, -0.7920994458748165E+00, 
    -0.9496878040642897E+00, -0.5494991726029023E+00, 
    0.9415698074379265E+00, 0.1549790418658323E+00, 
    -0.9014007983667046E+00, 0.4493885995928425E+00, 
    0.7095579485335214E+00, 0.8736487328036152E+00, 
    0.3264862525811960E-01, 0.1103424065156252E+00, 
    0.5637539669697259E+00, -0.2945721905451988E+00, 
    0.8603062445393470E+00, 0.7507058748066501E+00, 
    0.2780639680677569E+00, -0.7047555601238906E+00, 
    0.6295944388807833E+00, 0.9095323899548271E+00, 
    -0.2964917498074751E+00, -0.8851424078101294E+00, 
    0.7219563370216030E+00, 0.9230518045006707E+00, 
    -0.5631064783902135E+00, 0.9323802418890973E+00, 
    -0.9285744095618478E+00, 0.3446313474069607E-01, 
    -0.5285124795902665E+00, 0.6564153211279398E+00, 
    -0.9670556758395561E+00, -0.9746827810441796E+00, 
    0.1497145509400534E+00, -0.5126157047506470E-01, 
    -0.5771162069337331E+00, 0.9882152664587445E+00, 
    -0.8664710684798441E+00, -0.1554480104048762E+00 };
  double ys[84] = {
    0.9380477277286925E+00, -0.3903208068416588E+00, 
    0.5539681565347837E+00, -0.9087002846713521E+00, 
    -0.6321385589135142E+00, 0.3182961009431599E+00, 
    -0.8236851175776945E+00, -0.8104327207737747E+00, 
    0.9326434060163864E+00, -0.4047646882920218E+00, 
    0.9551344910189569E+00, -0.8100975319542825E+00, 
    -0.2881465340470265E+00, 0.9393091372522897E+00, 
    0.7616291888714266E+00, -0.3792679791993436E+00, 
    -0.8198273601663221E+00, -0.6277875734879305E+00, 
    0.3601156105251353E+00, -0.6500193633731873E+00, 
    -0.1345149716413504E+00, 0.5453436457811436E+00, 
    0.9440747796649167E+00, -0.9592627915208141E+00, 
    -0.7303420755929871E-01, 0.8946123590730062E+00, 
    0.6073216328066294E-01, 0.2144997262319605E+00, 
    -0.1311566771466323E+00, -0.7186229429207273E+00, 
    -0.6749210564844365E+00, 0.9691840758487141E+00, 
    -0.6464321973320226E+00, 0.9661037973313142E+00, 
    0.5022924249551910E+00, 0.2011928667597957E+00, 
    -0.3711890466829925E+00, -0.7861623327708182E+00, 
    0.9831486697529307E+00, 0.1008010491222713E+00, 
    -0.6108831589158809E+00, 0.9546709186773804E+00, 
    -0.9380477277286926E+00, 0.3903208068416588E+00, 
    -0.5539681565347837E+00, 0.9087002846713521E+00, 
    0.6321385589135142E+00, -0.3182961009431599E+00, 
    0.8236851175776946E+00, 0.8104327207737747E+00, 
    -0.9326434060163864E+00, 0.4047646882920218E+00, 
    -0.9551344910189569E+00, 0.8100975319542825E+00, 
    0.2881465340470266E+00, -0.9393091372522897E+00, 
    -0.7616291888714266E+00, 0.3792679791993436E+00, 
    0.8198273601663221E+00, 0.6277875734879305E+00, 
    -0.3601156105251352E+00, 0.6500193633731873E+00, 
    0.1345149716413503E+00, -0.5453436457811436E+00, 
    -0.9440747796649167E+00, 0.9592627915208141E+00, 
    0.7303420755929881E-01, -0.8946123590730062E+00, 
    -0.6073216328066286E-01, -0.2144997262319605E+00, 
    0.1311566771466323E+00, 0.7186229429207273E+00, 
    0.6749210564844365E+00, -0.9691840758487141E+00, 
    0.6464321973320226E+00, -0.9661037973313142E+00, 
    -0.5022924249551909E+00, -0.2011928667597957E+00, 
    0.3711890466829925E+00, 0.7861623327708182E+00, 
    -0.9831486697529307E+00, -0.1008010491222714E+00, 
    0.6108831589158809E+00, -0.9546709186773804E+00 };
  double zs[84] = {
    0.8462611120890736E+00, 0.3547611988896818E+00, 
    0.5262314749145942E-02, 0.7237196107659919E+00, 
    -0.6425814432975854E-01, 0.1803373124222549E+00, 
    -0.7073699162717987E+00, -0.3664557776165722E+00, 
    0.1016191298800923E+00, -0.4999025408268907E+00, 
    0.9145661354846110E+00, 0.6962396504425696E+00, 
    0.4859552045824018E+00, -0.9060323834426769E+00, 
    -0.1822470452558573E+00, 0.7016603041053766E+00, 
    -0.6424102177565745E+00, -0.7635438798952285E+00, 
    0.5612530704657699E+00, 0.9279296695895975E+00, 
    -0.8932000387343114E+00, 0.9433964716928596E+00, 
    -0.3855860228559991E+00, 0.8267083885708849E+00, 
    -0.2085191046423283E+00, 0.2867455231288707E+00, 
    -0.6303696136436402E+00, -0.9480808897684662E+00, 
    -0.8851760571455803E+00, 0.3677243413181205E+00, 
    0.5707827021209619E+00, 0.3308168458013439E+00, 
    0.8942058612668059E+00, 0.8491966568619232E+00, 
    0.6997888056053100E+00, -0.8582587623738921E+00, 
    0.9980904532339944E+00, -0.9929983508078530E+00, 
    -0.2607372382203097E+00, -0.2706521005603409E-01, 
    -0.9886712504226799E+00, -0.9981234658396176E+00, 
    -0.8462611120890736E+00, -0.3547611988896818E+00, 
    -0.5262314749145976E-02, -0.7237196107659919E+00, 
    0.6425814432975864E-01, -0.1803373124222550E+00, 
    0.7073699162717988E+00, 0.3664557776165722E+00, 
    -0.1016191298800922E+00, 0.4999025408268908E+00, 
    -0.9145661354846109E+00, -0.6962396504425696E+00, 
    -0.4859552045824018E+00, 0.9060323834426768E+00, 
    0.1822470452558573E+00, -0.7016603041053766E+00, 
    0.6424102177565745E+00, 0.7635438798952286E+00, 
    -0.5612530704657699E+00, -0.9279296695895975E+00, 
    0.8932000387343114E+00, -0.9433964716928596E+00, 
    0.3855860228559991E+00, -0.8267083885708848E+00, 
    0.2085191046423283E+00, -0.2867455231288708E+00, 
    0.6303696136436402E+00, 0.9480808897684662E+00, 
    0.8851760571455803E+00, -0.3677243413181204E+00, 
    -0.5707827021209619E+00, -0.3308168458013439E+00, 
    -0.8942058612668059E+00, -0.8491966568619232E+00, 
    -0.6997888056053101E+00, 0.8582587623738921E+00, 
    -0.9980904532339944E+00, 0.9929983508078530E+00, 
    0.2607372382203098E+00, 0.2706521005603400E-01, 
    0.9886712504226799E+00, 0.9981234658396176E+00 };
  double ws[84] = {
    0.1322946412898770E-01, 0.7207978700091233E-01, 
    0.5018350497921869E-01, 0.2510540976305218E-01, 
    0.6381617344667434E-01, 0.5297826014608787E-01, 
    0.1107387819925981E-01, 0.4482499482382119E-01, 
    0.1156633883765873E-01, 0.8060095080318432E-01, 
    0.4983228901803450E-02, 0.3768915302119238E-01, 
    0.6105445139039273E-01, 0.6907527307462437E-02, 
    0.6642161794579643E-01, 0.6883290145149486E-01, 
    0.3765112561164340E-01, 0.5084567296013873E-01, 
    0.4165679702303787E-01, 0.1924373869332567E-01, 
    0.4590245342578355E-01, 0.2073261975977018E-01, 
    0.2506136162764897E-01, 0.6618443166679632E-02, 
    0.1133796667663198E+00, 0.2117476937280500E-01, 
    0.6105702586211274E-01, 0.1302086249789518E-01, 
    0.4287461768786696E-01, 0.2726768060699983E-01, 
    0.2693097474826901E-01, 0.2634041203981413E-01, 
    0.3577055224953832E-01, 0.1180436625397797E-01, 
    0.1895051399977877E-01, 0.1378259645552230E-01, 
    0.1414623939342679E-01, 0.1101406610125121E-01, 
    0.1973564062736466E-01, 0.2303337346288448E-01, 
    0.8809606694498163E-02, 0.6060743137742029E-02, 
    0.1322946412898769E-01, 0.7207978700091235E-01, 
    0.5018350497921865E-01, 0.2510540976305218E-01, 
    0.6381617344667437E-01, 0.5297826014608786E-01, 
    0.1107387819925980E-01, 0.4482499482382119E-01, 
    0.1156633883765873E-01, 0.8060095080318433E-01, 
    0.4983228901803448E-02, 0.3768915302119238E-01, 
    0.6105445139039275E-01, 0.6907527307462438E-02, 
    0.6642161794579643E-01, 0.6883290145149484E-01, 
    0.3765112561164340E-01, 0.5084567296013873E-01, 
    0.4165679702303786E-01, 0.1924373869332567E-01, 
    0.4590245342578356E-01, 0.2073261975977018E-01, 
    0.2506136162764896E-01, 0.6618443166679627E-02, 
    0.1133796667663198E+00, 0.2117476937280501E-01, 
    0.6105702586211271E-01, 0.1302086249789518E-01, 
    0.4287461768786696E-01, 0.2726768060699984E-01, 
    0.2693097474826901E-01, 0.2634041203981412E-01, 
    0.3577055224953834E-01, 0.1180436625397797E-01, 
    0.1895051399977878E-01, 0.1378259645552230E-01, 
    0.1414623939342678E-01, 0.1101406610125121E-01, 
    0.1973564062736465E-01, 0.2303337346288442E-01, 
    0.8809606694498165E-02, 0.6060743137742026E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[116] = {
    -0.7715535323463264E+00, 0.3118087350366243E-01, 
    0.7202332131024749E+00, 0.4477007231810433E+00, 
    0.4990977036666501E+00, 0.7276854447788592E+00, 
    0.5641552377783545E-01, -0.6804449077772355E+00, 
    -0.3147755705468945E-01, 0.9656886599696646E+00, 
    -0.5809498997793892E+00, -0.2289254671435773E+00, 
    0.7876771526454306E+00, -0.4068670890700770E+00, 
    -0.8972365649862598E+00, -0.3532895192280955E+00, 
    -0.9567571227498393E+00, -0.1216832508787129E+00, 
    0.5700762717211539E+00, 0.5051003810086915E+00, 
    0.9701289189066667E+00, -0.8467694257635688E+00, 
    -0.9792766700328013E+00, 0.9239050051895750E+00, 
    0.3182661293383706E+00, 0.2128643425667921E+00, 
    -0.9491244031998681E+00, 0.7992362723313403E+00, 
    -0.2460884124375043E+00, 0.1342871260315162E+00, 
    -0.2373280028547328E+00, -0.5865009282285342E+00, 
    0.4626440026822248E+00, 0.3297957744728636E+00, 
    -0.9420531528791128E+00, -0.4341804961950992E+00, 
    -0.9225452596590470E+00, 0.6871688089857060E+00, 
    0.2448760905505635E+00, 0.5741688284173365E+00, 
    -0.8594595285871465E+00, -0.5554936117139078E+00, 
    -0.6494964119774425E+00, 0.5732782989920666E+00, 
    0.3752071332901741E-02, 0.6435232718547543E+00, 
    0.8238997691962470E+00, 0.6248144351825363E-01, 
    -0.5415679713524282E+00, -0.6912095128076139E-01, 
    0.4168395206084785E+00, 0.7671340542346751E+00, 
    0.9499134274201164E+00, 0.6612862602874549E+00, 
    0.4988292488539028E+00, -0.3980774278492163E+00, 
    0.2143015046710522E+00, -0.8245870486521025E+00, 
    0.9547288142651947E+00, -0.8256353544128002E+00, 
    0.8043198685414620E+00, -0.6534154451891698E+00, 
    -0.9452902258113793E+00, -0.6217330929720845E+00, 
    -0.9665133496506871E+00, 0.1521701688630693E+00, 
    -0.9371839994014188E+00, 0.9495343847338622E+00, 
    -0.3463029817783037E+00, -0.7895101820283607E+00, 
    -0.8159632548392504E-01, 0.2855115686118832E+00, 
    0.9535458359476985E+00, -0.2998645589788054E+00, 
    0.8771573781780831E+00, 0.5455017995441601E+00, 
    0.9166715004534769E-01, 0.9320968927813329E+00, 
    -0.4371224793112961E+00, -0.8347365350568826E+00, 
    -0.6455342940659587E+00, 0.8555205523715417E+00, 
    -0.2072617599943643E+00, -0.2113416898671757E+00, 
    -0.9486945106212309E+00, -0.9718805025138264E+00, 
    -0.1940500299791751E-02, 0.8641147432354730E+00, 
    -0.7571701781562347E+00, 0.4947987121201377E+00, 
    0.2223337905054378E+00, 0.9496122896148295E+00, 
    0.8640108234726359E+00, 0.7680642862414261E+00, 
    0.4625223529706101E+00, 0.9726683521171305E+00, 
    -0.7964752179472936E+00, -0.7061109658827748E+00, 
    -0.8721608252584128E+00, -0.6991324850192889E+00, 
    0.9831710538535281E+00, 0.3063724124769707E+00, 
    0.9886078377543190E+00, -0.6436734499882067E+00, 
    -0.1696864839354078E+00, -0.9581984720461905E+00, 
    0.6404741270364127E+00, -0.9954098177933381E+00, 
    -0.9998048311370020E+00, -0.3608374768729438E+00, 
    -0.4193491067426698E+00, -0.2591511667594117E+00, 
    0.5950366411277217E+00, 0.8362890578433639E+00, 
    -0.9989746911020813E+00, 0.9819256903295290E+00 };
  double ys[116] = {
    0.3846975186938493E+00, -0.6224498777026040E+00, 
    0.1248067846342114E+00, 0.8759502561317533E+00, 
    0.3319258914157911E+00, -0.5633971521527585E+00, 
    0.8534622243379885E+00, 0.9174064836961019E+00, 
    0.6613043682070445E-01, 0.8826351471721087E+00, 
    0.3382028556583220E+00, 0.7699804117571577E+00, 
    -0.3555309699183173E+00, 0.5163363528866298E+00, 
    0.9578029706408058E+00, 0.2101811715995720E+00, 
    -0.6596351003887070E+00, 0.2923658777723023E-01, 
    -0.2160515820663783E+00, -0.7325716420674306E+00, 
    -0.9320850402391798E+00, 0.2535810500632735E+00, 
    -0.3368010225765684E+00, 0.9401901116430753E+00, 
    0.7005289547119918E-01, -0.4318126034481278E+00, 
    -0.7540912429552113E+00, 0.7738338149718788E+00, 
    -0.2709207207715034E+00, -0.4558191555703927E+00, 
    -0.9010499919081281E+00, 0.5077988769546107E+00, 
    -0.3368014224099352E+00, -0.1248384730156639E+00, 
    -0.9290918735449377E+00, -0.7020388511680736E+00, 
    0.8388875287519405E+00, 0.3950015262650662E+00, 
    0.4698716842353510E+00, 0.7927540313537204E-01, 
    -0.4958172798073285E+00, 0.9158035761509341E-01, 
    0.6706864246998623E+00, 0.7918494529907293E+00, 
    -0.6319475420779336E+00, 0.6495695968654632E+00, 
    -0.9456860423718022E+00, -0.8989638999178421E+00, 
    -0.3735186503913843E+00, 0.4704635101750327E+00, 
    -0.8716224029261982E+00, 0.5966448338821769E+00, 
    0.2261581260452755E-01, -0.7330157165466431E+00, 
    -0.7471944764154522E+00, 0.9683805827306203E+00, 
    -0.5519164542917125E+00, -0.8971331427558276E+00, 
    -0.8173334319953630E+00, 0.8059812120592472E+00, 
    -0.6969456308437438E+00, -0.2340344455553026E+00, 
    -0.6410900230088896E+00, -0.7876018852914574E+00, 
    -0.2676406878547541E+00, 0.3656010440914417E+00, 
    0.7517669252844423E+00, 0.9299482185727768E+00, 
    0.8095208519693138E+00, -0.9105465384851248E-01, 
    0.5754559865383174E-01, 0.7987833771715905E+00, 
    0.4062222029583985E+00, -0.2917475248364298E+00, 
    0.8341598642457054E+00, 0.9565675836400581E+00, 
    0.9495501483452489E+00, -0.5027178324567889E+00, 
    -0.6094294091525906E+00, -0.2041586923215633E+00, 
    0.8735321107965007E+00, 0.8176411928945507E-01, 
    0.5798716521493317E+00, 0.9747533229160877E+00, 
    0.4095645350870887E+00, 0.2069198533007474E+00, 
    -0.8800808072488089E+00, -0.2586645109021137E+00, 
    -0.6284140921870087E+00, -0.9837391845052924E+00, 
    0.6897557755650331E+00, 0.6214627043620098E+00, 
    -0.9147932185865665E+00, 0.9683589014439020E+00, 
    0.9722720243469519E+00, 0.3545312716884971E+00, 
    -0.9265634816103003E+00, 0.9823681379288890E+00, 
    0.6027221797612410E+00, 0.1550378502770828E+00, 
    -0.7188423853605198E+00, -0.9880539732209580E+00, 
    -0.3538081678676123E+00, -0.9848881983673942E+00, 
    0.9503771401811550E+00, 0.9895536743716908E+00, 
    -0.9877214291859133E+00, 0.6282240608226360E-01, 
    0.8687973020715331E+00, -0.9993768714765972E+00, 
    -0.8816819132287129E+00, -0.9024706682879141E+00, 
    0.3779679396631332E+00, -0.1517734130495604E+00, 
    -0.9719850622545042E+00, -0.9258604067050853E+00 };
  double zs[116] = {
    -0.9149704001000553E+00, -0.5996714764210891E+00, 
    0.7907026496157833E+00, 0.9860392703559183E+00, 
    -0.4660715186534219E-01, -0.1393752230506733E+00, 
    0.8674767286359658E+00, 0.5467250815629719E+00, 
    -0.9763780125306059E+00, 0.5121364231430013E+00, 
    -0.6798617661241456E+00, 0.6371639769399203E+00, 
    0.2198362375814198E+00, -0.9280616749626068E+00, 
    0.7168598119569591E+00, 0.7368659134196429E+00, 
    -0.6061233370696343E+00, -0.7956884470646045E+00, 
    0.6572737265408953E+00, -0.9403846902844000E+00, 
    -0.6276502147372000E+00, 0.1269919055587806E+00, 
    0.4130533850992282E+00, 0.9393352681302030E+00, 
    0.8990311398428759E+00, -0.9212501242203712E+00, 
    0.5134369934005911E+00, 0.6109131534473928E+00, 
    -0.4181668086036411E+00, 0.5711408979823358E+00, 
    0.1117558004736949E+00, 0.3763350218603200E+00, 
    -0.5509292330117498E+00, -0.1436673317026735E-01, 
    0.8947995961383289E+00, 0.5293478643731013E+00, 
    0.9681968180532952E+00, 0.3548727821114203E+00, 
    0.6370668273720668E+00, -0.7744664564128696E+00, 
    -0.1174975410572717E+00, -0.2455457407049606E+00, 
    0.8684746243106054E+00, -0.3023424151630013E+00, 
    -0.8384309459289455E-01, 0.9285129719767410E+00, 
    -0.1479858042465700E+00, -0.7406305201907566E+00, 
    0.2042255456501517E+00, 0.9724059801463540E+00, 
    -0.3537644656876312E+00, -0.6281373274195593E+00, 
    0.4150304386610159E+00, 0.8808795896002272E+00, 
    0.3686704712426682E+00, 0.9460942659633094E+00, 
    0.9681464254739912E+00, 0.1328695747867607E+00, 
    -0.9698622847731962E+00, -0.9424723737502007E+00, 
    -0.6724135642153274E+00, -0.9637049413487057E+00, 
    -0.9654826579439407E+00, -0.3959451854631504E+00, 
    0.9553150406508956E+00, -0.4818859928369797E+00, 
    0.2896997475970485E+00, -0.5366726979456555E+00, 
    -0.6510037274971031E+00, 0.6185934968970651E+00, 
    0.3210185330943278E+00, 0.2400288119511415E+00, 
    0.8510838518999654E+00, 0.8754236562515145E+00, 
    -0.9470610777853996E+00, -0.8205588369782968E+00, 
    -0.4018976232642799E+00, 0.7751784036421725E+00, 
    -0.8147980972514376E+00, -0.6601638637102712E+00, 
    -0.1260885622828708E+00, -0.3084734983651212E+00, 
    -0.8179658955267972E-01, 0.2700716241482205E+00, 
    0.7694412333217397E+00, -0.8836683466383671E+00, 
    0.7669737715528676E+00, -0.9137690852575695E+00, 
    0.8605953113236016E+00, 0.8953797328026272E+00, 
    -0.8573304152865556E+00, -0.1765576765676193E-01, 
    0.5663399409423517E+00, 0.8216340568717669E-01, 
    0.7178509062543948E+00, -0.8248755417862220E+00, 
    -0.8408410046283714E+00, -0.7367428078031174E+00, 
    -0.4590219429006754E+00, 0.9831224519425289E+00, 
    0.1099197404645035E+00, 0.2714915831470573E+00, 
    -0.4771292875453570E+00, 0.6410524347781664E+00, 
    -0.9808727711140862E+00, -0.1799508573122650E+00, 
    -0.8648319917262105E+00, -0.1821178559958520E+00, 
    -0.7478825522596463E+00, -0.4398682537286034E+00, 
    0.9929654026734756E+00, -0.9978596306712945E+00, 
    -0.9997568603576099E+00, 0.9997078506812050E+00, 
    -0.3559815310105911E+00, 0.9905772316230947E+00 };
  double ws[116] = {
    0.8171123714265456E-02, 0.2263229388891795E-01, 
    0.1818677478454230E-01, 0.2980129445309951E-02, 
    0.4365304132620811E-01, 0.3243988845839241E-01, 
    0.1486444991990349E-01, 0.1409533365332470E-01, 
    0.1188012333452950E-01, 0.6000058026819710E-02, 
    0.3512238833735714E-01, 0.3098299265015942E-01, 
    0.3654756930376477E-01, 0.1847118507157365E-01, 
    0.5391012235288447E-02, 0.4074719736210774E-01, 
    0.1096224300316088E-01, 0.4081648743391735E-01, 
    0.4154456995823457E-01, 0.1320278412821872E-01, 
    0.4290449708656973E-02, 0.3535406067309509E-01, 
    0.1105110007450795E-01, 0.3228649704663156E-02, 
    0.2983835839394490E-01, 0.2451444151701184E-01, 
    0.1248137360474183E-01, 0.2177139751257475E-01, 
    0.6405596175169058E-01, 0.5434653731649390E-01, 
    0.3076101092170598E-01, 0.4899566852520795E-01, 
    0.5225345103961254E-01, 0.7166342203347331E-01, 
    0.3920984475787756E-02, 0.4137916341593958E-01, 
    0.3600172354407817E-02, 0.4755213384458860E-01, 
    0.5073623632857935E-01, 0.3944403305910345E-01, 
    0.3387054595129749E-01, 0.6277110525242416E-01, 
    0.2166338571221270E-01, 0.3747883688963586E-01, 
    0.6163661481210551E-01, 0.1666526164418470E-01, 
    0.1395113427475628E-01, 0.2347682793313850E-01, 
    0.6306339243839931E-01, 0.1556571594055037E-01, 
    0.3387962780929576E-01, 0.3268836056103089E-01, 
    0.2278931777487743E-01, 0.1962500582231501E-01, 
    0.4456347447201976E-01, 0.5576706603430151E-02, 
    0.1594423386163151E-01, 0.2035708135880945E-01, 
    0.3103484011599505E-02, 0.9044075785104641E-02, 
    0.2689758906199793E-01, 0.1586801972478047E-01, 
    0.5176595113431329E-02, 0.3858960164224838E-01, 
    0.5739931404252700E-02, 0.7158286709420723E-01, 
    0.1871157285391988E-01, 0.8194850948531980E-02, 
    0.3678122219124356E-01, 0.4250046085039932E-01, 
    0.8513849699040082E-01, 0.4993413782858739E-01, 
    0.1233822876422000E-01, 0.3941817987846824E-01, 
    0.7113418117354728E-02, 0.1194652514212529E-01, 
    0.2508023865178747E-01, 0.1725533956116014E-01, 
    0.3737010567296990E-01, 0.3668781342194600E-01, 
    0.3343163204554410E-01, 0.4496577144136969E-01, 
    0.7409539266110132E-01, 0.1793429116081015E-01, 
    0.1653960634734396E-01, 0.9190332852452304E-02, 
    0.2809796363000793E-01, 0.1777992760135489E-01, 
    0.2391632963314971E-01, 0.5749099213189723E-02, 
    0.3433774681484043E-01, 0.2274673898154641E-01, 
    0.1536942439600953E-01, 0.1430939257823392E-01, 
    0.1282919496053713E-01, 0.1121380440794668E-01, 
    0.1184676390611188E-01, 0.7982226989051386E-02, 
    0.3477385886958077E-01, 0.1190123605011125E-01, 
    0.1168229078012223E-01, 0.1343804099688911E-01, 
    0.1228891741931907E-01, 0.1036505178017157E-01, 
    0.5937174344287885E-02, 0.4133423282546477E-02, 
    0.6371010607497968E-02, 0.1313727761639429E-01, 
    0.3669992824757512E-02, 0.1078959526158836E-01, 
    0.7230671036845260E-02, 0.6406235552064114E-02, 
    0.1096476832530329E-01, 0.7917718162474640E-02, 
    0.3572933439799740E-02, 0.1541652661200926E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[130] = {
    0.4870864101401451E+00, -0.2374760548181477E-01, 
    -0.7693492847298501E+00, 0.9412905403090202E+00, 
    0.6911155930532483E+00, 0.7898594135857243E+00, 
    -0.7034504617756053E+00, 0.1116726714264779E+00, 
    -0.6392758447394441E+00, -0.8512618100723829E+00, 
    -0.6550935879174717E+00, -0.2012982100456206E-01, 
    0.2379647994307194E+00, -0.8601912819787483E+00, 
    -0.9232873932846022E+00, -0.1900861086818608E+00, 
    0.1929412791278925E-01, 0.7258444409271682E+00, 
    -0.8143255678954792E+00, 0.2946629110201063E+00, 
    -0.4322866997565211E+00, 0.9734181389553355E+00, 
    -0.8655112784848059E+00, -0.9486962755224994E+00, 
    -0.8639618309045872E+00, 0.4192783763456998E+00, 
    -0.6493595219283761E-01, -0.3132633581182750E+00, 
    -0.7107686198927078E+00, 0.1584847796276057E+00, 
    0.7962238411361664E+00, 0.7050836730161649E+00, 
    0.8581636882302615E+00, 0.9437903485136390E+00, 
    -0.9469415138106567E+00, 0.3693463627427371E+00, 
    -0.4815016166149078E+00, 0.6225323271333136E+00, 
    -0.1216882752396847E+00, 0.4861399922133740E+00, 
    -0.6588662398008378E+00, 0.8733371851679250E+00, 
    -0.5553986064438454E+00, 0.7785003762786959E+00, 
    -0.9809869293544266E+00, -0.9200430029249042E+00, 
    0.3304920778295448E+00, 0.3661210439903886E+00, 
    0.6228621851378309E+00, -0.6948006452003259E+00, 
    -0.9310681159823312E+00, 0.1927229482124844E+00, 
    -0.9395909910668155E+00, -0.8918846544837020E+00, 
    -0.5046559139572119E+00, -0.5501818316586664E+00, 
    -0.9816538375537931E+00, 0.5810740673825343E-01, 
    -0.9917476914762476E+00, 0.1619246931269223E+00, 
    0.4044274445105526E+00, 0.9788832392715772E+00, 
    0.2896262066021671E+00, -0.9979092625656895E+00, 
    -0.9953891308905386E+00, -0.4870864101401453E+00, 
    0.2374760548181472E-01, 0.7693492847298502E+00, 
    -0.9412905403090200E+00, -0.6911155930532483E+00, 
    -0.7898594135857244E+00, 0.7034504617756053E+00, 
    -0.1116726714264779E+00, 0.6392758447394442E+00, 
    0.8512618100723829E+00, 0.6550935879174717E+00, 
    0.2012982100456194E-01, -0.2379647994307196E+00, 
    0.8601912819787482E+00, 0.9232873932846023E+00, 
    0.1900861086818609E+00, -0.1929412791278927E-01, 
    -0.7258444409271683E+00, 0.8143255678954792E+00, 
    -0.2946629110201063E+00, 0.4322866997565212E+00, 
    -0.9734181389553355E+00, 0.8655112784848059E+00, 
    0.9486962755224994E+00, 0.8639618309045871E+00, 
    -0.4192783763456998E+00, 0.6493595219283750E-01, 
    0.3132633581182749E+00, 0.7107686198927078E+00, 
    -0.1584847796276057E+00, -0.7962238411361664E+00, 
    -0.7050836730161649E+00, -0.8581636882302615E+00, 
    -0.9437903485136390E+00, 0.9469415138106567E+00, 
    -0.3693463627427371E+00, 0.4815016166149078E+00, 
    -0.6225323271333137E+00, 0.1216882752396847E+00, 
    -0.4861399922133739E+00, 0.6588662398008377E+00, 
    -0.8733371851679250E+00, 0.5553986064438454E+00, 
    -0.7785003762786958E+00, 0.9809869293544267E+00, 
    0.9200430029249042E+00, -0.3304920778295449E+00, 
    -0.3661210439903887E+00, -0.6228621851378309E+00, 
    0.6948006452003260E+00, 0.9310681159823312E+00, 
    -0.1927229482124843E+00, 0.9395909910668155E+00, 
    0.8918846544837020E+00, 0.5046559139572118E+00, 
    0.5501818316586664E+00, 0.9816538375537931E+00, 
    -0.5810740673825336E-01, 0.9917476914762476E+00, 
    -0.1619246931269224E+00, -0.4044274445105527E+00, 
    -0.9788832392715772E+00, -0.2896262066021671E+00, 
    0.9979092625656895E+00, 0.9953891308905385E+00 };
  double ys[130] = {
    -0.4125999616297640E+00, -0.3577512709695355E+00, 
    -0.4586571513966052E+00, -0.8809700870302296E+00, 
    0.6141953888878622E+00, -0.1098686336047354E-01, 
    -0.2992063293298189E+00, -0.5784543682045333E+00, 
    0.1546499308353232E+00, 0.4744749166521733E+00, 
    0.6850785787653841E+00, 0.3498227222335277E+00, 
    -0.3706543545438595E-01, 0.9489055210853593E+00, 
    -0.8258138282634399E+00, 0.8310925689243290E+00, 
    -0.2955546221774640E+00, -0.2274560569718798E+00, 
    -0.7326258705118668E+00, 0.6810827567696704E+00, 
    -0.1839671034672870E+00, 0.6471305653882479E+00, 
    0.4545843865988449E+00, 0.2460664244382031E-01, 
    0.1874945363707773E+00, -0.5593766578223933E+00, 
    -0.9791589336665545E+00, -0.6283651952819381E+00, 
    0.5755680670959202E+00, 0.2389500706336686E-01, 
    -0.8335047675899898E+00, 0.3839430357533925E+00, 
    -0.9053291152689182E+00, 0.3178523468094013E+00, 
    0.1336071580769665E-01, -0.2394172411300733E+00, 
    -0.9590729155294019E-01, 0.9385425560623889E+00, 
    0.7891616010496042E+00, -0.8504780607583017E+00, 
    -0.9475181655259435E+00, 0.3020909284966213E+00, 
    -0.9580268216765446E+00, -0.9705093056943978E+00, 
    0.7530029765856929E+00, -0.9658747597521009E+00, 
    0.8589417682173183E+00, 0.5576404948411705E+00, 
    -0.9193264723330097E+00, -0.7573749854479314E+00, 
    0.7034268760540773E+00, 0.9510178717882097E+00, 
    -0.9731097679191421E+00, -0.8322328408621713E+00, 
    0.6019176765989178E+00, -0.4092215952731691E+00, 
    -0.5199973715576629E+00, 0.8581532655611123E+00, 
    0.2298310606072314E+00, -0.8504134103774321E+00, 
    -0.9915729124045047E+00, 0.8271450907199669E+00, 
    -0.9986629495001677E+00, 0.5886511512141286E+00, 
    0.9903359117471948E+00, 0.4125999616297641E+00, 
    0.3577512709695354E+00, 0.4586571513966052E+00, 
    0.8809700870302297E+00, -0.6141953888878622E+00, 
    0.1098686336047354E-01, 0.2992063293298190E+00, 
    0.5784543682045332E+00, -0.1546499308353234E+00, 
    -0.4744749166521733E+00, -0.6850785787653841E+00, 
    -0.3498227222335276E+00, 0.3706543545438599E-01, 
    -0.9489055210853593E+00, 0.8258138282634397E+00, 
    -0.8310925689243290E+00, 0.2955546221774641E+00, 
    0.2274560569718798E+00, 0.7326258705118668E+00, 
    -0.6810827567696704E+00, 0.1839671034672870E+00, 
    -0.6471305653882478E+00, -0.4545843865988448E+00, 
    -0.2460664244382029E-01, -0.1874945363707773E+00, 
    0.5593766578223933E+00, 0.9791589336665545E+00, 
    0.6283651952819381E+00, -0.5755680670959203E+00, 
    -0.2389500706336671E-01, 0.8335047675899898E+00, 
    -0.3839430357533924E+00, 0.9053291152689182E+00, 
    -0.3178523468094013E+00, -0.1336071580769660E-01, 
    0.2394172411300734E+00, 0.9590729155294021E-01, 
    -0.9385425560623889E+00, -0.7891616010496042E+00, 
    0.8504780607583017E+00, 0.9475181655259436E+00, 
    -0.3020909284966213E+00, 0.9580268216765446E+00, 
    0.9705093056943976E+00, -0.7530029765856929E+00, 
    0.9658747597521009E+00, -0.8589417682173183E+00, 
    -0.5576404948411705E+00, 0.9193264723330097E+00, 
    0.7573749854479312E+00, -0.7034268760540774E+00, 
    -0.9510178717882097E+00, 0.9731097679191421E+00, 
    0.8322328408621713E+00, -0.6019176765989177E+00, 
    0.4092215952731690E+00, 0.5199973715576629E+00, 
    -0.8581532655611123E+00, -0.2298310606072313E+00, 
    0.8504134103774322E+00, 0.9915729124045047E+00, 
    -0.8271450907199669E+00, 0.9986629495001677E+00, 
    -0.5886511512141284E+00, -0.9903359117471947E+00 };
  double zs[130] = {
    0.2510431589001390E+00, -0.4730112580495535E+00, 
    -0.8296554957987250E+00, 0.6042839304817023E+00, 
    -0.4354086139851164E+00, 0.7512711569117528E+00, 
    0.2580885268523619E+00, 0.7181855545442590E+00, 
    0.7258539919433724E+00, 0.5394422830180516E+00, 
    -0.3665937773187872E+00, -0.4287198875778774E+00, 
    -0.1391272946798063E+00, -0.2939636504110460E+00, 
    -0.8894944904099870E-01, 0.6509390582294805E+00, 
    -0.8295699704108563E+00, 0.9037134226828011E-03, 
    0.7543687304209290E+00, 0.7789156932196367E+00, 
    0.6834706616686410E+00, -0.2649602815663506E+00, 
    -0.6640360044523207E+00, 0.3841143167045120E+00, 
    0.9462427733191926E+00, -0.3779968053092916E+00, 
    0.7871565548827959E+00, 0.8953414085684869E+00, 
    -0.9510332692943319E+00, -0.9747087096347502E+00, 
    -0.6936472418684265E+00, -0.9725884697577191E+00, 
    -0.9605395058198848E+00, -0.7896659742434409E+00, 
    -0.9585765042296354E+00, 0.8915927965605104E+00, 
    -0.5323083877061187E+00, -0.9089350785066914E+00, 
    -0.8104031924041570E-01, 0.7186845896034352E+00, 
    0.7635825638284459E-01, 0.2475045724812137E+00, 
    -0.8458514443440915E+00, 0.9116524510899568E+00, 
    -0.9030128536354431E+00, 0.5786900078194604E+00, 
    -0.4876904077471467E+00, 0.9522425964634561E-01, 
    -0.2288092843749514E+00, -0.4933688094107099E+00, 
    0.6791308826513548E-01, 0.3643415241007317E+00, 
    -0.5743612468117230E+00, -0.9410337888158190E+00, 
    0.9256166802221193E+00, -0.9794608522159433E+00, 
    -0.7050722652352325E+00, 0.9838945816912786E+00, 
    -0.3456879759150185E+00, 0.9876394980104540E+00, 
    -0.8354581552777850E+00, -0.9782936349079702E+00, 
    0.3136198827810360E+00, 0.8542937538183871E+00, 
    0.4938372507517853E+00, -0.2510431589001390E+00, 
    0.4730112580495534E+00, 0.8296554957987250E+00, 
    -0.6042839304817023E+00, 0.4354086139851165E+00, 
    -0.7512711569117528E+00, -0.2580885268523619E+00, 
    -0.7181855545442589E+00, -0.7258539919433724E+00, 
    -0.5394422830180515E+00, 0.3665937773187872E+00, 
    0.4287198875778774E+00, 0.1391272946798063E+00, 
    0.2939636504110459E+00, 0.8894944904099869E-01, 
    -0.6509390582294805E+00, 0.8295699704108564E+00, 
    -0.9037134226827184E-03, -0.7543687304209290E+00, 
    -0.7789156932196368E+00, -0.6834706616686410E+00, 
    0.2649602815663507E+00, 0.6640360044523208E+00, 
    -0.3841143167045120E+00, -0.9462427733191926E+00, 
    0.3779968053092917E+00, -0.7871565548827959E+00, 
    -0.8953414085684869E+00, 0.9510332692943317E+00, 
    0.9747087096347502E+00, 0.6936472418684265E+00, 
    0.9725884697577191E+00, 0.9605395058198848E+00, 
    0.7896659742434410E+00, 0.9585765042296354E+00, 
    -0.8915927965605104E+00, 0.5323083877061187E+00, 
    0.9089350785066915E+00, 0.8104031924041578E-01, 
    -0.7186845896034351E+00, -0.7635825638284455E-01, 
    -0.2475045724812137E+00, 0.8458514443440915E+00, 
    -0.9116524510899567E+00, 0.9030128536354431E+00, 
    -0.5786900078194604E+00, 0.4876904077471467E+00, 
    -0.9522425964634557E-01, 0.2288092843749515E+00, 
    0.4933688094107098E+00, -0.6791308826513545E-01, 
    -0.3643415241007316E+00, 0.5743612468117228E+00, 
    0.9410337888158191E+00, -0.9256166802221192E+00, 
    0.9794608522159431E+00, 0.7050722652352326E+00, 
    -0.9838945816912786E+00, 0.3456879759150185E+00, 
    -0.9876394980104540E+00, 0.8354581552777850E+00, 
    0.9782936349079702E+00, -0.3136198827810360E+00, 
    -0.8542937538183870E+00, -0.4938372507517853E+00 };
  double ws[130] = {
    0.3140524074284941E-01, 0.3981029681429966E-01, 
    0.1485446788498034E-01, 0.6244986380811999E-02, 
    0.2738837680317996E-01, 0.2157943912866181E-01, 
    0.3558093745023824E-01, 0.3212728314491558E-01, 
    0.3092422246858361E-01, 0.2343742079178917E-01, 
    0.3088702409027052E-01, 0.5359064880099586E-01, 
    0.6482379336742157E-01, 0.9118048913013883E-02, 
    0.1301789370977779E-01, 0.2655541823663328E-01, 
    0.3453297251160139E-01, 0.4411844784877157E-01, 
    0.1649537886932118E-01, 0.2909276637092253E-01, 
    0.4360788133113655E-01, 0.1053270755633742E-01, 
    0.2277866730141684E-01, 0.1989571475016660E-01, 
    0.1091026629930172E-01, 0.4918510213252013E-01, 
    0.8074020956961596E-02, 0.2335968213222367E-01, 
    0.1235756243300519E-01, 0.1499095494625247E-01, 
    0.1701688982544987E-01, 0.1042658914550453E-01, 
    0.4196082934549206E-02, 0.1415599321876098E-01, 
    0.6647555404704302E-02, 0.3129197633426045E-01, 
    0.5718706313414747E-01, 0.8214012330106898E-02, 
    0.4800378437465311E-01, 0.2511168810726739E-01, 
    0.1848081803186647E-01, 0.3560225217379276E-01, 
    0.9719778448629727E-02, 0.4568238333579292E-02, 
    0.3878678624072500E-02, 0.6252926394948910E-02, 
    0.3458846268652847E-01, 0.6444154518560745E-01, 
    0.2421665637314795E-01, 0.3403670963655146E-01, 
    0.2101561337278881E-01, 0.2321281908800168E-01, 
    0.4914553718764686E-02, 0.6665400363164401E-02, 
    0.2231229502988629E-01, 0.1251715201913546E-01, 
    0.9481650136297978E-02, 0.8059627544061664E-02, 
    0.1180468380164612E-01, 0.7820767064272208E-02, 
    0.6958745625819240E-02, 0.2234878561766285E-02, 
    0.1024235452666042E-01, 0.5000286844891411E-02, 
    0.2657409809447695E-02, 0.3140524074284939E-01, 
    0.3981029681429966E-01, 0.1485446788498035E-01, 
    0.6244986380811991E-02, 0.2738837680317996E-01, 
    0.2157943912866178E-01, 0.3558093745023824E-01, 
    0.3212728314491559E-01, 0.3092422246858360E-01, 
    0.2343742079178915E-01, 0.3088702409027051E-01, 
    0.5359064880099586E-01, 0.6482379336742154E-01, 
    0.9118048913013897E-02, 0.1301789370977778E-01, 
    0.2655541823663329E-01, 0.3453297251160137E-01, 
    0.4411844784877157E-01, 0.1649537886932117E-01, 
    0.2909276637092255E-01, 0.4360788133113656E-01, 
    0.1053270755633743E-01, 0.2277866730141684E-01, 
    0.1989571475016661E-01, 0.1091026629930171E-01, 
    0.4918510213252012E-01, 0.8074020956961608E-02, 
    0.2335968213222367E-01, 0.1235756243300520E-01, 
    0.1499095494625249E-01, 0.1701688982544985E-01, 
    0.1042658914550453E-01, 0.4196082934549209E-02, 
    0.1415599321876098E-01, 0.6647555404704300E-02, 
    0.3129197633426043E-01, 0.5718706313414745E-01, 
    0.8214012330106886E-02, 0.4800378437465309E-01, 
    0.2511168810726739E-01, 0.1848081803186646E-01, 
    0.3560225217379277E-01, 0.9719778448629732E-02, 
    0.4568238333579302E-02, 0.3878678624072491E-02, 
    0.6252926394948914E-02, 0.3458846268652847E-01, 
    0.6444154518560744E-01, 0.2421665637314794E-01, 
    0.3403670963655148E-01, 0.2101561337278882E-01, 
    0.2321281908800168E-01, 0.4914553718764693E-02, 
    0.6665400363164389E-02, 0.2231229502988629E-01, 
    0.1251715201913550E-01, 0.9481650136297967E-02, 
    0.8059627544061670E-02, 0.1180468380164611E-01, 
    0.7820767064272205E-02, 0.6958745625819234E-02, 
    0.2234878561766282E-02, 0.1024235452666043E-01, 
    0.5000286844891406E-02, 0.2657409809447699E-02 };


  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[172] = {
    0.3724711910241946E-01, 0.4142659262186393E+00, 
    0.5177982339713144E+00, 0.3604423376125923E+00, 
    0.8256117236986897E+00, 0.4265362215919388E+00, 
    -0.5205630040142383E+00, 0.5476096231885307E+00, 
    -0.3908413316164978E-01, -0.6869205136187482E+00, 
    -0.8269790632507413E+00, -0.3077403136300738E+00, 
    -0.4871197763974943E+00, -0.6602942921439789E+00, 
    -0.1996644422759564E+00, -0.2612979535504019E+00, 
    0.5138865516799961E+00, -0.8650532980097219E+00, 
    0.5098273324253541E+00, 0.7674795450461877E+00, 
    0.1907820276114263E+00, -0.7413413727558544E+00, 
    0.8578428692155201E+00, -0.9483245969337234E+00, 
    0.6622911479023382E+00, 0.1082533914200937E+00, 
    0.9727081876507619E+00, 0.5564069612990988E+00, 
    0.6174233759718846E+00, -0.7840704914833193E+00, 
    -0.3863554802732841E+00, -0.7642023424993235E+00, 
    -0.7532542829872437E+00, -0.8018433723045076E+00, 
    0.3770659555409062E+00, -0.9491540989399579E+00, 
    -0.9382558523915473E+00, -0.9808828605653619E+00, 
    -0.4844689691834917E+00, 0.9256039531559204E+00, 
    -0.8709952382834715E-01, -0.5396392548668583E+00, 
    0.1368025517956275E+00, 0.3374041492353974E+00, 
    0.7068715546921336E+00, -0.4085110348118505E+00, 
    -0.9681798506556836E+00, -0.1549674800732925E+00, 
    0.9150816905013831E+00, 0.5375252305267589E-01, 
    -0.5579823240421485E+00, 0.8969270740962718E-02, 
    0.7014120009018806E+00, 0.1912491108449427E+00, 
    0.5828147015903538E+00, 0.9521221772289220E+00, 
    -0.1994976199453465E+00, -0.8137190118036577E+00, 
    -0.5666026693931812E+00, -0.7926656084386846E+00, 
    -0.6502780637374249E+00, 0.1603598833824042E+00, 
    -0.9646753205843002E+00, 0.1874185821317708E+00, 
    0.8078893967819246E+00, 0.9821198889247245E+00, 
    -0.9134820612323216E+00, 0.9315563412192052E+00, 
    0.9673991240002898E-02, -0.5917684536649436E+00, 
    -0.5910952945368380E+00, -0.9004251112937950E+00, 
    0.4666785637176430E+00, -0.9279892538634408E+00, 
    0.5165215046136072E+00, -0.6522976391673276E+00, 
    -0.4635151732777533E+00, -0.2216103796843620E+00, 
    0.9720841352166502E+00, 0.8697669064712771E+00, 
    -0.7844145244283243E+00, 0.5653272793694807E+00, 
    0.5670262199937089E+00, -0.8760558996804383E+00, 
    -0.6568661047943578E+00, -0.6431016384228920E-01, 
    0.2699910056033774E+00, -0.3193920506608590E+00, 
    0.5594292856531002E+00, -0.8625783835676980E+00, 
    0.9642093306328916E+00, 0.9564210050320944E+00, 
    -0.5256008382082430E+00, 0.6831869769415302E+00, 
    -0.5564497515417263E+00, -0.8665833611119687E+00, 
    0.2119314545737717E+00, 0.6499134147160973E+00, 
    0.8416070492408827E+00, -0.8860664019820089E+00, 
    0.8630542601243912E+00, 0.4006206217350135E+00, 
    -0.3255422740709495E+00, 0.3889737632509294E+00, 
    -0.4031111086863209E+00, -0.1026555049543291E+00, 
    -0.2154587129488632E+00, -0.9587635011613380E+00, 
    -0.6781849005122116E+00, 0.6431694453965539E-02, 
    0.5990879315696672E+00, 0.2604246381903720E+00, 
    -0.2205299895937927E+00, 0.7360523872792837E+00, 
    0.9141664583516047E+00, 0.6460945155487459E+00, 
    0.8261067491745480E+00, 0.2582806673903185E+00, 
    0.7992011614665788E+00, 0.1952759431439480E+00, 
    -0.1404821582745676E+00, -0.9719800694682158E+00, 
    -0.9082932030606901E+00, -0.8069621339496692E+00, 
    0.5138739584252263E+00, -0.4508670463586423E+00, 
    0.7460163089423278E+00, 0.8559732256227004E+00, 
    0.2619517671437559E+00, -0.7178069404428022E+00, 
    0.7768409833572711E+00, 0.9566352889085462E+00, 
    -0.9129931899686473E+00, -0.2027074937002416E+00, 
    -0.9767056206307869E+00, -0.1953671488795780E+00, 
    0.9356675487713187E+00, 0.2733263440546579E-01, 
    0.8541785574421273E+00, -0.8509341178434021E+00, 
    -0.8389589448490465E+00, -0.9114666802875271E+00, 
    -0.2369922070544546E+00, 0.8352761372704256E+00, 
    0.6990375376349731E+00, -0.9825936001862040E+00, 
    0.9890254163849905E+00, 0.4299873858805139E+00, 
    0.9084398364888665E+00, -0.4942757972549205E+00, 
    0.1276694562023527E+00, 0.8433252017493162E+00, 
    0.9812363871492615E+00, -0.5112868612432091E+00, 
    -0.9919630673934499E+00, 0.9888971621532768E+00, 
    0.9906786057375682E+00, 0.9291016994895615E+00, 
    -0.9930351279012178E+00, -0.4368492693215694E+00, 
    -0.9926580930992508E+00, 0.9810341338091584E+00, 
    0.2829487311079641E+00, 0.9951522144360470E+00, 
    -0.9467788426833306E+00, 0.6231428267482866E+00, 
    -0.9966694864427889E+00, -0.5545516866973424E+00, 
    -0.2200546924714422E+00, 0.9855361391127346E+00, 
    0.1242997310669405E+00, -0.2411466505640982E-01 };
  double ys[172] = {
    0.4603796517412539E+00, 0.7850928815883216E+00, 
    0.5139455507370757E+00, 0.9223890225451771E+00, 
    0.8434522866134452E+00, 0.6558862622340167E-01, 
    -0.8559434795169205E+00, -0.8987656715485196E+00, 
    0.4503144430486827E+00, -0.2618729087942862E+00, 
    0.9452894304667214E+00, -0.7280261316218042E+00, 
    -0.2338055905041149E+00, -0.4867015686038832E+00, 
    -0.8999180244059493E+00, 0.8333786735111544E-01, 
    -0.8337854025340903E+00, -0.3449380786557670E+00, 
    -0.3067078127421212E+00, -0.3509360716139491E+00, 
    -0.8655219480345354E+00, 0.1796817314118926E-01, 
    -0.6893079251914401E+00, -0.1756439435663739E+00, 
    -0.1729573648651443E+00, -0.8312838577052795E+00, 
    0.8465219729055495E+00, -0.7218336696759129E+00, 
    -0.8968669124158696E-01, 0.1072951525625464E+00, 
    0.2674684823317927E+00, 0.4494767718386262E+00, 
    -0.8467786081052485E+00, 0.4339653736633297E+00, 
    0.3042380260051257E+00, 0.5590699107221699E+00, 
    -0.5068967389653932E+00, -0.2405947127531271E+00, 
    0.8437874725090614E+00, -0.9017396521006380E+00, 
    0.8482922518824456E+00, -0.6591372127582090E+00, 
    0.1703921286880628E+00, 0.9563525832577533E+00, 
    0.9739918780282115E+00, -0.7115110177234589E+00, 
    -0.9751728996744083E+00, -0.9631468160421225E+00, 
    0.6009999468106125E+00, -0.5288962426219667E+00, 
    0.4595578349581757E+00, 0.9432670812899058E+00, 
    0.7399828097883627E+00, -0.4141604710166898E+00, 
    -0.4084225192010409E+00, -0.3386216837958897E+00, 
    -0.6001665103475902E+00, -0.8022895195480112E+00, 
    -0.3398023662469815E+00, 0.8186709874482936E+00, 
    -0.9024700874428341E+00, -0.8294025892178960E+00, 
    0.3513024770242288E+00, 0.8650597368535108E+00, 
    -0.8463680274683520E+00, -0.7300261016180669E+00, 
    -0.9758837733429508E+00, 0.9751575986713164E+00, 
    0.7126310344823681E+00, 0.9451735803803867E+00, 
    0.1370096023427370E+00, -0.9445191743918876E+00, 
    -0.5896755668296407E+00, -0.1197749258942398E+00, 
    0.5242071257866383E+00, 0.9051990022441897E+00, 
    0.1957941852390964E+00, 0.5375124860468565E+00, 
    0.8696281579006412E+00, 0.9780419188336278E+00, 
    -0.4600808747557152E+00, -0.6775300402478708E+00, 
    0.2563064109512353E+00, -0.6053293404401969E+00, 
    0.7731104690643289E+00, -0.9394448657791105E+00, 
    -0.1604465004198996E+00, 0.9502436474270409E+00, 
    0.8261506546206309E+00, -0.8150926456217005E+00, 
    0.9776649598956791E-01, 0.8412957281351358E+00, 
    0.7214572271897897E+00, -0.8985831295061654E+00, 
    -0.9335211793245841E+00, 0.7067653633526974E-02, 
    -0.3615585073957835E+00, 0.2972904877318600E+00, 
    -0.4095443751975350E+00, 0.9558903952389940E+00, 
    0.8951478996227191E-01, 0.6611506169154207E+00, 
    0.6670036283630431E+00, -0.4785541651013403E-01, 
    0.5349559346436485E+00, 0.9899002101148441E+00, 
    -0.3900488015542639E+00, -0.7458885499515072E+00, 
    0.5665446610649311E+00, 0.1548891409840480E-01, 
    -0.9844705836707888E+00, -0.6746117169036595E+00, 
    -0.1470918231599819E+00, -0.5584752636915873E+00, 
    -0.6804678115064715E+00, 0.8909927345534838E+00, 
    0.6374291217859047E+00, 0.9675265148508291E+00, 
    0.4420117190898502E+00, 0.7027112199253245E+00, 
    -0.3799739422082143E+00, 0.1856369603681191E+00, 
    0.7746845502717523E+00, 0.1761706438449975E+00, 
    0.9850881697540705E+00, -0.1129086405431629E+00, 
    0.2271646891547787E+00, -0.1127987048407091E+00, 
    -0.9200033438731134E+00, -0.6626998668546163E+00, 
    0.9270956556038703E+00, -0.2470253489012092E+00, 
    0.7729357655654772E+00, -0.7356540617546378E+00, 
    -0.3212372190294076E+00, 0.8111308146975724E+00, 
    0.6299555026169626E+00, 0.3951231967843170E-03, 
    -0.8706689354769571E+00, -0.9692856262585002E+00, 
    0.9756175749168980E+00, 0.6261197595970289E+00, 
    0.9652057112084984E+00, -0.9811255607149989E+00, 
    0.8079753140237501E+00, 0.4087152028895200E+00, 
    0.4607837944330232E+00, -0.9858814414729271E+00, 
    0.1870441924567301E-01, 0.9880110117128622E+00, 
    0.4089802625067766E+00, -0.9609137755818135E+00, 
    0.2225216531350550E+00, -0.9861887397735654E+00, 
    0.8950245843584548E+00, -0.6572193683282491E+00, 
    -0.9304667586038193E+00, 0.4475608353766626E+00, 
    -0.7058486792297883E+00, -0.9958976269497638E+00, 
    -0.8836382486319199E+00, -0.9857261377173898E+00, 
    -0.9951450397546436E+00, -0.5152481827199886E+00, 
    0.9915918444310511E+00, -0.4698538990980113E+00, 
    0.8727201920285046E+00, -0.5170681183537231E+00, 
    -0.9631512995815722E+00, 0.9953119825440425E+00, 
    0.4988843272874602E+00, 0.3927236168362344E+00 };
  double zs[172] = {
    0.6938946186791187E+00, -0.5765923304339572E+00, 
    0.8051907870892863E+00, -0.8472224677887984E-01, 
    0.7318018811656403E-01, 0.5983817879906067E-01, 
    0.3384018846541420E+00, -0.4920664066906827E+00, 
    -0.7373417229245569E-01, 0.3303528888179985E+00, 
    0.1893355905396206E+00, -0.8380648944600999E+00, 
    0.9822767264128165E-01, -0.2594364391801900E+00, 
    0.1930055941915506E+00, -0.8498731486396089E+00, 
    -0.9480961804554179E-01, -0.8114075510008492E+00, 
    0.3445956986162864E+00, 0.6500251287969905E+00, 
    0.6911242761503570E+00, -0.5570201903906392E+00, 
    -0.4405333944883792E+00, -0.1725718027999678E+00, 
    -0.3345083448640582E+00, -0.2509356761788316E+00, 
    0.9317616130493932E+00, -0.9454774309944995E+00, 
    -0.7098026696648545E+00, -0.1134846343849767E+00, 
    -0.1315593065233067E+00, 0.2820827580654706E+00, 
    -0.8747191955455609E-01, -0.8330074038703938E+00, 
    -0.3660613547587568E+00, -0.9580088134576359E+00, 
    -0.4411464922737413E+00, 0.3949263385441928E+00, 
    -0.5640899102030010E+00, 0.2226335297251030E+00, 
    -0.9257134547555910E+00, 0.5715726078629155E+00, 
    -0.6345978485849848E+00, -0.9304344738485810E+00, 
    0.3590760844142656E+00, -0.4438907426169790E+00, 
    0.8671008381591335E+00, 0.6403756616117305E+00, 
    0.2778456794573282E+00, -0.7154487635298511E+00, 
    -0.5078493941700393E+00, -0.6363728062752509E+00, 
    -0.2513684642710016E+00, 0.8974010278661878E+00, 
    -0.9918409439634794E+00, -0.4534667592885467E+00, 
    0.4398586891128420E-01, -0.6130676703797119E+00, 
    -0.7477961790309134E+00, -0.7858333946639103E+00, 
    0.7974411188351147E+00, -0.9156598617795586E+00, 
    0.1715184891688095E+00, 0.1907051832001387E+00, 
    -0.7772414059750807E+00, -0.9276175191759890E-01, 
    -0.5405362134262204E+00, -0.1103237214049949E+00, 
    -0.3367426994768587E+00, 0.5093651265106862E+00, 
    -0.9492669502709199E+00, -0.9401434781115950E+00, 
    -0.5980941978161105E+00, -0.9625175885989743E+00, 
    0.1167123901439443E+00, 0.8885399133952907E+00, 
    0.5227273173456383E+00, -0.7706144348078099E+00, 
    0.5114844827530561E+00, 0.8060019555582266E+00, 
    0.8011735992610899E+00, 0.8195422038055056E+00, 
    0.9395527898900553E+00, 0.2015513588196933E+00, 
    -0.4186783562866535E-01, -0.5712296785602576E+00, 
    -0.9030794957455209E+00, -0.1456519314475667E+00, 
    0.5686925252996401E+00, 0.9600262609268753E+00, 
    -0.7217322859826423E+00, -0.5209570186170079E+00, 
    -0.9717413123988232E+00, 0.4303819820507280E+00, 
    -0.8147382303071704E+00, 0.6244560595709048E+00, 
    -0.1937189914196895E+00, -0.9083093391078847E+00, 
    -0.8428140765049459E+00, -0.9764924706512431E+00, 
    -0.9460095279803191E-01, -0.7952403780336077E+00, 
    0.2978906730661770E+00, 0.7125070355101365E+00, 
    0.9665991806647677E+00, 0.3839275268438931E+00, 
    -0.9734009359633906E+00, 0.6433849545611667E+00, 
    0.7568446904202943E+00, 0.2498470880629417E+00, 
    0.7852684316315199E+00, 0.3907525133961428E+00, 
    -0.4393237281615760E+00, 0.5260748276088058E-01, 
    0.6241597079858533E+00, 0.9654604035241454E+00, 
    0.7855716043610358E+00, 0.7924746367580968E+00, 
    -0.5575727107236519E+00, 0.9341485578568330E+00, 
    0.6036661107749994E+00, -0.6578753396306457E+00, 
    0.9656622221606493E+00, 0.9535888144290141E+00, 
    -0.3940334702415435E+00, 0.8677859635129408E+00, 
    0.4525286376195224E+00, 0.8943730814713799E+00, 
    0.9733506381066803E+00, -0.9666093420715280E+00, 
    -0.7599491907772208E+00, 0.2870616863165352E+00, 
    0.5428594142592057E+00, 0.9000765734384685E+00, 
    0.9346435403112985E+00, 0.7253765854245153E+00, 
    -0.8907140353727465E+00, 0.9838621660462418E+00, 
    0.9382671102570804E+00, 0.4142977100280188E+00, 
    -0.4245033033835321E+00, -0.3429310836114446E+00, 
    0.9810366085861507E+00, -0.2580286066432419E+00, 
    -0.9859800907904950E+00, 0.8105778659854276E+00, 
    -0.1748292017546381E+00, -0.8047213401011649E+00, 
    -0.9868201233515799E+00, -0.8647620143857092E+00, 
    -0.9830541038366302E+00, -0.9825555113262869E+00, 
    0.6976904274215842E+00, -0.2151083099716469E+00, 
    0.4203550437156745E-01, -0.9338232453847534E+00, 
    -0.6702199613543625E+00, 0.9924356236360254E+00, 
    -0.8478229607550496E+00, 0.9583412329573594E+00, 
    -0.7904934590950125E-01, 0.7258320418569250E+00, 
    0.1347146763176373E+00, 0.9283115095600205E+00, 
    0.7983610958403323E+00, 0.9982250382393918E+00, 
    -0.7346265359815700E+00, 0.9993883513767151E+00, 
    -0.9964956614827338E+00, -0.9525278817342810E+00, 
    0.6030693170605685E+00, 0.7499893148930252E+00 };
  double ws[172] = {
    -0.1426541755470291E+00, 0.1385765810780117E-01, 
    0.1307994851680021E-01, 0.1073487590359840E-01, 
    0.1051260418442039E-01, 0.3275658093319255E-01, 
    0.1563154221413661E-01, 0.1207051149698361E-01, 
    0.3472763431163371E-01, 0.2607207312794260E-01, 
    0.7063494865440152E-02, 0.1455303522761317E-01, 
    0.3637446239280831E-01, 0.2741442134745744E-01, 
    0.1814853970042907E-01, 0.2211942357145999E-01, 
    0.2063146835573402E-01, 0.1187881032517930E-01, 
    0.3411740640767904E-01, 0.2044586424163537E-01, 
    0.1584398903577958E-01, 0.2546900742398995E-01, 
    0.1514164315096591E-01, 0.1383700180253657E-01, 
    0.3244335479410856E-01, 0.2480839462235662E-01, 
    0.1901942550392541E-02, 0.8435911241596891E-02, 
    0.2592855621096610E-01, 0.2917638594610470E-01, 
    0.4254452270715581E-01, 0.2672040615070632E-01, 
    0.1679948502767675E-01, 0.1441579204032182E-01, 
    0.4055353005194201E-01, 0.3562579203350458E-02, 
    0.1304869501588879E-01, 0.7974389607375001E-02, 
    0.1942169211659009E-01, 0.7696214255731893E-02, 
    0.9786582979889373E-02, 0.2656456655588142E-01, 
    0.3887064238866873E-01, 0.4865147372961578E-02, 
    0.7182928711371533E-02, 0.2992836230232769E-01, 
    0.1246769069332907E-02, 0.1018121368113786E-01, 
    0.1583310102000952E-01, 0.3101719294140454E-01, 
    0.3352533262979018E-01, 0.1304614239329544E-01, 
    0.2424339545109358E-01, 0.2062227856017715E-01, 
    0.4634332034618291E-02, 0.1316721066365703E-01, 
    0.4197768965655505E-01, 0.1445189210678322E-01, 
    0.2742773234129495E-01, 0.1147786426088701E-01, 
    0.1051388351726134E-01, 0.1177681391912994E-01, 
    0.1280321776997680E-01, 0.2638280742080748E-01, 
    0.1039146162515972E-01, 0.6350797568985297E-02, 
    0.3731595769989320E-02, 0.4022230391599942E-02, 
    0.3700138698657982E-01, 0.1223831728952482E-01, 
    0.1361018792408060E-01, 0.2432559937042062E-02, 
    0.3173506278445155E-01, 0.5299185434865185E-02, 
    0.4119660596888559E-01, 0.8140787806672345E-02, 
    0.4245764049999567E-01, 0.2966712362252400E-01, 
    0.5251158142618067E-02, 0.3117784323015366E-02, 
    0.1887240736415798E-01, 0.1982870468110188E-01, 
    0.1525468960328201E-01, 0.2145937561316863E-01, 
    0.2757411556050665E-01, 0.1580892051031026E-01, 
    0.2339199782020639E-01, 0.1652864001989000E-01, 
    0.2227125289030323E-01, 0.4542047728295454E-02, 
    0.1029753747302055E-01, 0.7588833799248887E-02, 
    0.7729962070932306E-02, 0.1686200103735939E-01, 
    0.9855191228454099E-02, 0.2298392421656214E-01, 
    0.5390281878584920E-01, 0.1777695356252309E-01, 
    0.1543307662800650E-01, 0.1565350686857245E-02, 
    0.2990218331164426E-01, 0.2492829949352890E-01, 
    0.4068595698602832E-01, 0.3904492542455568E-01, 
    0.1145081666124139E-01, 0.7199499019479311E-02, 
    0.1180596844175838E-01, 0.8490203829980697E-02, 
    0.2407429627173780E-01, 0.6013060579704863E-01, 
    0.4782368835976528E-02, 0.4070082427805708E-01, 
    0.5421074129521822E-01, 0.3495817844335809E-01, 
    0.1411530671886532E-01, 0.5360156659230995E-02, 
    0.1684523030402797E-01, 0.8957506702452017E-02, 
    0.2843883212961609E-01, 0.1548562016507650E-01, 
    0.4755338701620224E-01, 0.1061956205236735E-01, 
    0.4143328303679084E-02, 0.1091266243253667E-01, 
    0.7928114591801557E-02, 0.2853726563554603E-01, 
    0.3789221238377655E-01, 0.1479871521595558E-01, 
    0.5265072839453472E-02, 0.8371632675069558E-02, 
    0.9778033410801315E-02, 0.1741500947112345E-01, 
    0.1416545817648641E-01, 0.1911128041899456E-01, 
    0.4409869231823535E-02, 0.2671133974518937E-01, 
    0.8042675590771235E-02, 0.1102479574395968E-01, 
    0.5810239704806290E-02, 0.7613319884605023E-02, 
    0.6957655736787461E-02, 0.2091424497221789E-01, 
    0.2998965152610065E-02, 0.6589150225887293E-02, 
    0.4511407449002006E-02, 0.6488450322177047E-02, 
    0.8640589231658198E-02, 0.5930519400069490E-02, 
    0.4449535891375782E-02, 0.4479599530453893E-02, 
    0.1130073799626979E-01, 0.1782480610982645E-02, 
    0.9295402400959387E-02, 0.9494483239679656E-02, 
    0.3985977390728825E-02, 0.2723087629850607E-02, 
    0.2577391123867757E-02, 0.3091874244052285E-02, 
    0.3705458814905999E-02, 0.2140475474717141E-02, 
    0.4774963198255120E-02, 0.1630610059874353E-02, 
    0.9085452239351545E-02, 0.2955971697674354E-02, 
    0.2082908274049303E-02, 0.5720118558010348E-02, 
    0.2931930843304469E-02, 0.5718185266004892E-02, 
    0.2289842864717349E-02, 0.6408796456581044E-03, 
    0.1069581923634952E+00, 0.1069506278589204E+00 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
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
//    09 July 2014
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
//    Output, double X[3*N], the coordinates of the nodes.
//
//    Output, double W[N], the weights.
//
{
  double xs[190] = {
    0.4846214646897683E+00, 0.6027166661210369E-01, 
    0.5119120893947191E+00, -0.5797527815798613E+00, 
    0.8916317278515248E+00, 0.4427375637743738E+00, 
    -0.2409001206926237E+00, 0.4213881124359606E+00, 
    0.2454110415378226E+00, -0.1642469924955429E+00, 
    -0.5907939646559498E+00, -0.6671386201069288E+00, 
    0.7552290798850946E+00, 0.1001864121384495E+00, 
    -0.7219847932977880E+00, 0.2982096063117719E+00, 
    0.5073944886321697E+00, -0.4863692703613542E+00, 
    0.4639404445534597E+00, -0.8081387191848195E+00, 
    -0.5896154562963155E+00, 0.1349771790076746E+00, 
    -0.1459114981397958E+00, -0.7771044530879095E+00, 
    0.3049106464545794E+00, 0.7478398889947702E+00, 
    0.2593964431357709E+00, -0.8826393735985893E+00, 
    0.5355784681079047E+00, -0.1464672879292506E+00, 
    0.6066659900076390E+00, -0.8588774763360313E+00, 
    -0.9697378209073355E+00, -0.3207448401664478E+00, 
    0.3238825883494653E+00, -0.8911881463557618E+00, 
    0.9588421432293531E+00, -0.9459659824876482E-01, 
    0.4326888655727574E+00, -0.4031429838762605E+00, 
    -0.4784772616280669E+00, 0.8992518592510214E+00, 
    0.6586422425595307E+00, 0.8204427703547663E+00, 
    -0.5614110159343418E+00, 0.9520903109905504E-01, 
    0.5640256443183821E+00, -0.8042276926515493E+00, 
    0.6821201856074234E+00, -0.7682858721251195E+00, 
    0.7082541697630407E+00, -0.4618970778119202E+00, 
    -0.2168388877255139E+00, 0.5371020397039006E+00, 
    0.8885379409916848E+00, 0.6251053221154703E+00, 
    0.8755744976159847E+00, -0.1289396220206777E+00, 
    0.9646166796039118E+00, 0.2607898402928713E+00, 
    0.9703940236997658E+00, 0.7108751014443095E-01, 
    0.7712967479205184E+00, 0.8842270989720199E+00, 
    0.7991197473160693E+00, -0.3410291594283488E+00, 
    -0.4192756150578085E+00, 0.7230818664415763E+00, 
    0.1436938638381775E+00, -0.9757899278543728E+00, 
    0.8447822406339708E+00, 0.9613637788842541E+00, 
    0.8507427476805113E+00, -0.8848833562540662E+00, 
    0.9720070913947075E+00, -0.6259228134698737E+00, 
    0.9441521606893316E+00, 0.2886043547146998E+00, 
    0.9428658265051933E+00, 0.9499838153205579E+00, 
    -0.1175912491528422E+00, -0.9810321115887223E+00, 
    0.9305760666668829E+00, -0.9783274226607617E+00, 
    -0.9834890551889414E+00, 0.4453699774671322E-01, 
    -0.5306871050560035E+00, 0.9868369576014261E+00, 
    -0.5232830891127319E+00, -0.8918143354641236E+00, 
    0.8676179401449021E+00, 0.9964872079697851E+00, 
    0.9918319243897599E+00, 0.7865264660879382E+00, 
    -0.9241660476457195E-01, -0.4846214646897686E+00, 
    -0.6027166661210399E-01, -0.5119120893947190E+00, 
    0.5797527815798613E+00, -0.8916317278515249E+00, 
    -0.4427375637743741E+00, 0.2409001206926235E+00, 
    -0.4213881124359606E+00, -0.2454110415378225E+00, 
    0.1642469924955428E+00, 0.5907939646559499E+00, 
    0.6671386201069287E+00, -0.7552290798850947E+00, 
    -0.1001864121384497E+00, 0.7219847932977881E+00, 
    -0.2982096063117720E+00, -0.5073944886321696E+00, 
    0.4863692703613542E+00, -0.4639404445534597E+00, 
    0.8081387191848196E+00, 0.5896154562963155E+00, 
    -0.1349771790076746E+00, 0.1459114981397958E+00, 
    0.7771044530879095E+00, -0.3049106464545795E+00, 
    -0.7478398889947702E+00, -0.2593964431357710E+00, 
    0.8826393735985893E+00, -0.5355784681079047E+00, 
    0.1464672879292507E+00, -0.6066659900076391E+00, 
    0.8588774763360314E+00, 0.9697378209073355E+00, 
    0.3207448401664478E+00, -0.3238825883494653E+00, 
    0.8911881463557617E+00, -0.9588421432293531E+00, 
    0.9459659824876483E-01, -0.4326888655727574E+00, 
    0.4031429838762604E+00, 0.4784772616280668E+00, 
    -0.8992518592510215E+00, -0.6586422425595307E+00, 
    -0.8204427703547664E+00, 0.5614110159343418E+00, 
    -0.9520903109905507E-01, -0.5640256443183821E+00, 
    0.8042276926515494E+00, -0.6821201856074234E+00, 
    0.7682858721251195E+00, -0.7082541697630407E+00, 
    0.4618970778119202E+00, 0.2168388877255139E+00, 
    -0.5371020397039005E+00, -0.8885379409916848E+00, 
    -0.6251053221154704E+00, -0.8755744976159846E+00, 
    0.1289396220206776E+00, -0.9646166796039118E+00, 
    -0.2607898402928713E+00, -0.9703940236997658E+00, 
    -0.7108751014443089E-01, -0.7712967479205184E+00, 
    -0.8842270989720198E+00, -0.7991197473160693E+00, 
    0.3410291594283487E+00, 0.4192756150578084E+00, 
    -0.7230818664415765E+00, -0.1436938638381774E+00, 
    0.9757899278543728E+00, -0.8447822406339707E+00, 
    -0.9613637788842542E+00, -0.8507427476805113E+00, 
    0.8848833562540663E+00, -0.9720070913947076E+00, 
    0.6259228134698737E+00, -0.9441521606893316E+00, 
    -0.2886043547146998E+00, -0.9428658265051933E+00, 
    -0.9499838153205579E+00, 0.1175912491528422E+00, 
    0.9810321115887223E+00, -0.9305760666668829E+00, 
    0.9783274226607617E+00, 0.9834890551889414E+00, 
    -0.4453699774671330E-01, 0.5306871050560035E+00, 
    -0.9868369576014261E+00, 0.5232830891127319E+00, 
    0.8918143354641237E+00, -0.8676179401449020E+00, 
    -0.9964872079697851E+00, -0.9918319243897600E+00, 
    -0.7865264660879382E+00, 0.9241660476457186E-01 };
  double ys[190] = {
    -0.2262155685815833E+00, -0.5802224500744104E+00, 
    -0.8812418910175070E+00, 0.4577182827511190E-01, 
    0.7310521233845833E+00, 0.4830906487655276E+00, 
    -0.8061293440666314E+00, -0.1078387681721650E+00, 
    0.8648936466868999E+00, 0.1116488353049886E+00, 
    -0.9117169599889896E+00, 0.7671850694134570E+00, 
    0.8098232031170816E+00, 0.9416109016304105E+00, 
    -0.8571607992505791E+00, -0.8294951885681802E+00, 
    -0.3287584215285467E-01, -0.9716511178926364E+00, 
    0.4508573465614439E+00, -0.9829067455694369E+00, 
    -0.3886494734052406E+00, 0.6210687867431628E+00, 
    0.7058206696430964E+00, 0.8874742900432545E+00, 
    0.2936260801274535E-01, 0.4728249466868350E+00, 
    -0.7999771197371979E+00, 0.9406572379268412E-01, 
    -0.4746762454784017E+00, 0.5158509951163104E+00, 
    -0.8256099261094548E+00, -0.5791703340444312E+00, 
    0.4466546017661202E+00, -0.1927320854548640E+00, 
    -0.4345261425945900E+00, 0.5425166397866776E-02, 
    0.9375941025615940E+00, -0.3420778054235573E+00, 
    0.7075968971274725E+00, 0.6496964078624228E+00, 
    0.4004462354823037E+00, -0.6968996840182875E+00, 
    0.3055994655170601E+00, -0.5169213379049986E+00, 
    0.9549786512838463E+00, -0.8791296893067777E+00, 
    0.8303000567130658E+00, -0.6377483111957468E+00, 
    0.8739025562234195E+00, -0.1693477428011365E+00, 
    -0.2351350601654628E+00, -0.6299185210533546E+00, 
    0.1604152098962463E+00, -0.5905825901125751E+00, 
    0.1734591472408230E+00, 0.2417993452840939E+00, 
    -0.7167382075250317E+00, -0.5500881365309197E+00, 
    -0.2400335193040073E+00, 0.7704776594613698E+00, 
    -0.5549501648642681E+00, 0.2871796808122397E+00, 
    -0.2817721114790634E+00, 0.9549722793117502E+00, 
    -0.9499810867427769E+00, 0.9910163892213160E+00, 
    -0.9801994414073382E+00, -0.4857741672686406E+00, 
    -0.9028324174014650E+00, -0.7689895270308303E+00, 
    -0.8885921033548654E+00, 0.2255470425271027E+00, 
    -0.7650697965238176E-01, 0.9735656776897391E+00, 
    -0.2305869558790821E+00, 0.7981924334040106E+00, 
    -0.6544414972640588E+00, -0.9790786374374271E+00, 
    0.9251057242786117E+00, 0.6489915714062549E+00, 
    -0.9627047612899647E+00, 0.9597706404726861E+00, 
    0.5373877582566533E+00, 0.8266106540930876E+00, 
    -0.8547964831867243E+00, 0.9558026916663981E+00, 
    0.9856879028440860E+00, -0.9065724188604554E+00, 
    -0.9845747454374670E+00, 0.9891294476252082E+00, 
    0.9881191368012432E+00, 0.3020236952873330E+00, 
    0.3856116113348522E+00, 0.7292572061308893E+00, 
    -0.1882646447022738E+00, 0.2262155685815840E+00, 
    0.5802224500744105E+00, 0.8812418910175072E+00, 
    -0.4577182827511186E-01, -0.7310521233845831E+00, 
    -0.4830906487655275E+00, 0.8061293440666314E+00, 
    0.1078387681721650E+00, -0.8648936466868998E+00, 
    -0.1116488353049886E+00, 0.9117169599889896E+00, 
    -0.7671850694134570E+00, -0.8098232031170816E+00, 
    -0.9416109016304105E+00, 0.8571607992505792E+00, 
    0.8294951885681802E+00, 0.3287584215285466E-01, 
    0.9716511178926363E+00, -0.4508573465614438E+00, 
    0.9829067455694369E+00, 0.3886494734052407E+00, 
    -0.6210687867431628E+00, -0.7058206696430964E+00, 
    -0.8874742900432544E+00, -0.2936260801274522E-01, 
    -0.4728249466868351E+00, 0.7999771197371979E+00, 
    -0.9406572379268432E-01, 0.4746762454784017E+00, 
    -0.5158509951163104E+00, 0.8256099261094548E+00, 
    0.5791703340444312E+00, -0.4466546017661200E+00, 
    0.1927320854548639E+00, 0.4345261425945901E+00, 
    -0.5425166397866769E-02, -0.9375941025615940E+00, 
    0.3420778054235572E+00, -0.7075968971274725E+00, 
    -0.6496964078624228E+00, -0.4004462354823038E+00, 
    0.6968996840182876E+00, -0.3055994655170602E+00, 
    0.5169213379049985E+00, -0.9549786512838463E+00, 
    0.8791296893067777E+00, -0.8303000567130657E+00, 
    0.6377483111957469E+00, -0.8739025562234195E+00, 
    0.1693477428011367E+00, 0.2351350601654628E+00, 
    0.6299185210533547E+00, -0.1604152098962464E+00, 
    0.5905825901125752E+00, -0.1734591472408231E+00, 
    -0.2417993452840938E+00, 0.7167382075250318E+00, 
    0.5500881365309197E+00, 0.2400335193040072E+00, 
    -0.7704776594613698E+00, 0.5549501648642682E+00, 
    -0.2871796808122396E+00, 0.2817721114790634E+00, 
    -0.9549722793117502E+00, 0.9499810867427769E+00, 
    -0.9910163892213160E+00, 0.9801994414073382E+00, 
    0.4857741672686406E+00, 0.9028324174014649E+00, 
    0.7689895270308303E+00, 0.8885921033548654E+00, 
    -0.2255470425271028E+00, 0.7650697965238172E-01, 
    -0.9735656776897391E+00, 0.2305869558790820E+00, 
    -0.7981924334040106E+00, 0.6544414972640588E+00, 
    0.9790786374374271E+00, -0.9251057242786117E+00, 
    -0.6489915714062550E+00, 0.9627047612899647E+00, 
    -0.9597706404726861E+00, -0.5373877582566532E+00, 
    -0.8266106540930876E+00, 0.8547964831867243E+00, 
    -0.9558026916663980E+00, -0.9856879028440859E+00, 
    0.9065724188604554E+00, 0.9845747454374671E+00, 
    -0.9891294476252082E+00, -0.9881191368012432E+00, 
    -0.3020236952873331E+00, -0.3856116113348522E+00, 
    -0.7292572061308894E+00, 0.1882646447022739E+00 };
  double zs[190] = {
    0.9718923364663833E+00, -0.6995279408119103E+00, 
    0.4077821712986755E+00, -0.1984488505968536E+00, 
    -0.1635978697790265E+00, -0.9519012083224356E+00, 
    -0.7709973566567310E+00, 0.6575891838559990E+00, 
    -0.2431469413756804E+00, -0.2154126099125328E+00, 
    0.2917708722416409E+00, -0.7837093593401955E+00, 
    0.5503667466636896E-01, 0.6034320247288605E+00, 
    -0.5348830267053428E+00, 0.6967890074634113E+00, 
    -0.6004372124809728E+00, -0.7978753021283712E+00, 
    0.7111633203919522E+00, -0.4796592851863407E+00, 
    -0.1066936069842310E+00, -0.4937724895804262E-02, 
    -0.3617523885223834E+00, 0.7864082151686236E+00, 
    0.8896593937798930E+00, -0.3003206340764916E+00, 
    -0.8604298301675264E+00, -0.9764052940319004E+00, 
    0.4627660680406495E+00, -0.8466728944778892E+00, 
    -0.9693428983943551E+00, 0.7807520536361380E+00, 
    0.2410080019231437E+00, 0.2275963934164933E+00, 
    -0.1058345123272911E+00, -0.4676766976465098E-01, 
    0.1025372637039018E+00, 0.6347152609726029E+00, 
    0.3752067881462010E+00, 0.5361771132313945E+00, 
    0.9695076036670979E+00, 0.6398892617052058E+00, 
    0.9539776797271410E+00, 0.1946562772645561E+00, 
    0.5161190440947852E+00, -0.2830148609612029E+00, 
    0.9517662443841286E+00, -0.7967109780188084E+00, 
    -0.7304370458992153E+00, -0.5338593513474772E+00, 
    -0.3179475898860780E+00, 0.5547973598962450E+00, 
    0.8313518081338683E+00, 0.9654245596905445E+00, 
    -0.5628516000227574E+00, -0.8613302773849665E+00, 
    -0.4483881659919932E+00, -0.9517697705546500E+00, 
    -0.7234814778587770E+00, -0.8875078228929659E+00, 
    0.9243612035652934E+00, 0.4601935017159294E+00, 
    0.8090143745218911E+00, 0.9163662035802967E+00, 
    0.4382382798746217E+00, 0.8459991721751295E+00, 
    -0.1366146439527482E+00, -0.7720184760904137E+00, 
    0.9707987948940664E+00, 0.5206831670219514E+00, 
    0.9588210246069581E+00, 0.8185288301633218E+00, 
    -0.9577089524402139E+00, 0.1925548617062123E+00, 
    0.4997201390064559E+00, 0.5406623725365253E-01, 
    -0.9375810397798050E+00, 0.1306044624286102E+00, 
    -0.8929122835431057E+00, 0.9831724010468225E+00, 
    0.6385058885805324E+00, -0.7720170702032852E+00, 
    0.3549666743399563E+00, -0.1714024163539602E+00, 
    -0.6784170912845590E+00, 0.9880880116453774E+00, 
    -0.8468765296700308E+00, -0.6791166853273773E+00, 
    0.9378710595411145E+00, 0.9586435791263295E+00, 
    -0.4583535837605982E+00, -0.1539477252636465E+00, 
    -0.9251749829870214E+00, -0.9949152686131382E+00, 
    0.9985713487179375E+00, -0.9718923364663833E+00, 
    0.6995279408119103E+00, -0.4077821712986756E+00, 
    0.1984488505968536E+00, 0.1635978697790266E+00, 
    0.9519012083224355E+00, 0.7709973566567310E+00, 
    -0.6575891838559990E+00, 0.2431469413756804E+00, 
    0.2154126099125329E+00, -0.2917708722416409E+00, 
    0.7837093593401955E+00, -0.5503667466636887E-01, 
    -0.6034320247288605E+00, 0.5348830267053428E+00, 
    -0.6967890074634113E+00, 0.6004372124809728E+00, 
    0.7978753021283711E+00, -0.7111633203919521E+00, 
    0.4796592851863406E+00, 0.1066936069842310E+00, 
    0.4937724895804275E-02, 0.3617523885223834E+00, 
    -0.7864082151686236E+00, -0.8896593937798931E+00, 
    0.3003206340764915E+00, 0.8604298301675265E+00, 
    0.9764052940319005E+00, -0.4627660680406495E+00, 
    0.8466728944778892E+00, 0.9693428983943552E+00, 
    -0.7807520536361380E+00, -0.2410080019231437E+00, 
    -0.2275963934164933E+00, 0.1058345123272911E+00, 
    0.4676766976465104E-01, -0.1025372637039019E+00, 
    -0.6347152609726029E+00, -0.3752067881462010E+00, 
    -0.5361771132313945E+00, -0.9695076036670979E+00, 
    -0.6398892617052058E+00, -0.9539776797271410E+00, 
    -0.1946562772645561E+00, -0.5161190440947853E+00, 
    0.2830148609612028E+00, -0.9517662443841286E+00, 
    0.7967109780188085E+00, 0.7304370458992153E+00, 
    0.5338593513474772E+00, 0.3179475898860781E+00, 
    -0.5547973598962450E+00, -0.8313518081338683E+00, 
    -0.9654245596905445E+00, 0.5628516000227572E+00, 
    0.8613302773849664E+00, 0.4483881659919933E+00, 
    0.9517697705546500E+00, 0.7234814778587769E+00, 
    0.8875078228929659E+00, -0.9243612035652934E+00, 
    -0.4601935017159294E+00, -0.8090143745218911E+00, 
    -0.9163662035802967E+00, -0.4382382798746217E+00, 
    -0.8459991721751295E+00, 0.1366146439527482E+00, 
    0.7720184760904137E+00, -0.9707987948940665E+00, 
    -0.5206831670219514E+00, -0.9588210246069581E+00, 
    -0.8185288301633218E+00, 0.9577089524402139E+00, 
    -0.1925548617062122E+00, -0.4997201390064560E+00, 
    -0.5406623725365260E-01, 0.9375810397798049E+00, 
    -0.1306044624286103E+00, 0.8929122835431057E+00, 
    -0.9831724010468225E+00, -0.6385058885805324E+00, 
    0.7720170702032852E+00, -0.3549666743399564E+00, 
    0.1714024163539601E+00, 0.6784170912845590E+00, 
    -0.9880880116453774E+00, 0.8468765296700308E+00, 
    0.6791166853273775E+00, -0.9378710595411145E+00, 
    -0.9586435791263295E+00, 0.4583535837605982E+00, 
    0.1539477252636465E+00, 0.9251749829870214E+00, 
    0.9949152686131382E+00, -0.9985713487179375E+00 };
  double ws[190] = {
    0.2931082855526895E-02, 0.1466168295291819E-01, 
    0.1008190603381851E-01, 0.2521840249289902E-01, 
    0.9545986541148931E-02, 0.7815725861454997E-02, 
    0.1225157612225792E-01, 0.2516512639883486E-01, 
    0.1749166437590727E-01, 0.4056629885298555E-01, 
    0.1199350194114567E-01, 0.1133390863336471E-01, 
    0.1485090749543295E-01, 0.1026773216326104E-01, 
    0.1201022409690237E-01, 0.1547758419389476E-01, 
    0.2898163426866567E-01, 0.5077671719424524E-02, 
    0.2462022242995652E-01, 0.3796424546988349E-02, 
    0.3310665553884919E-01, 0.3508140615002558E-01, 
    0.2950704004490954E-01, 0.8017748658106175E-02, 
    0.1975278275621659E-01, 0.2559496685146834E-01, 
    0.1349736154064143E-01, 0.4377789154192862E-02, 
    0.3252004086526810E-01, 0.2138516935436983E-01, 
    0.5015100019990201E-02, 0.1268265369185493E-01, 
    0.9980860281711945E-02, 0.4271671806230670E-01, 
    0.4274198226264674E-01, 0.2231057111860769E-01, 
    0.4634811584165995E-02, 0.3628925840326913E-01, 
    0.2963883283828190E-01, 0.3004053280506377E-01, 
    0.9656867842652010E-02, 0.1208027996271507E-01, 
    0.1067097545113014E-01, 0.2483163417257641E-01, 
    0.1052923525628832E-01, 0.2351898637367898E-01, 
    0.7064042544802274E-02, 0.1458074471394978E-01, 
    0.1255548203713305E-01, 0.2826856991152390E-01, 
    0.3484550434776072E-01, 0.3029815093584674E-01, 
    0.2852474545169975E-01, 0.9146009437075924E-02, 
    0.1990115803594484E-01, 0.2069177392023117E-01, 
    0.1613360621469677E-01, 0.1339344756044318E-01, 
    0.9274781087469087E-02, 0.1530433781694090E-01, 
    0.3917613624484169E-02, 0.4743890758461208E-01, 
    0.2001420743439168E-01, 0.2832855798485957E-02, 
    0.9111027190743101E-02, 0.3380728915063295E-02, 
    0.9357526881990973E-02, 0.2181586678118927E-01, 
    0.5497095242202157E-02, 0.6401916148758926E-02, 
    0.3861833303734157E-02, 0.8708380742935694E-02, 
    0.8616911070460580E-02, 0.5860715399693166E-02, 
    0.1134352478246400E-01, 0.2833705081651916E-01, 
    0.4870087478700985E-02, 0.1085506021481164E-01, 
    0.3085543406158564E-02, 0.2320131607643815E-02, 
    0.1217435492763373E-01, 0.1876868717373962E-02, 
    0.1764889584855365E-01, 0.6738482260511419E-02, 
    0.4000917701843583E-02, 0.2530717864603938E-02, 
    0.4454840968844286E-02, 0.2917862788830831E-02, 
    0.3031055423950877E-02, 0.1091305108499325E-02, 
    0.4268181631337476E-02, 0.6826137445231240E-02, 
    0.3247017230887978E-02, 0.3899933571265500E-02, 
    0.7578127425388513E-02, 0.2931082855526887E-02, 
    0.1466168295291819E-01, 0.1008190603381851E-01, 
    0.2521840249289902E-01, 0.9545986541148931E-02, 
    0.7815725861455011E-02, 0.1225157612225792E-01, 
    0.2516512639883484E-01, 0.1749166437590725E-01, 
    0.4056629885298554E-01, 0.1199350194114567E-01, 
    0.1133390863336472E-01, 0.1485090749543296E-01, 
    0.1026773216326105E-01, 0.1201022409690238E-01, 
    0.1547758419389476E-01, 0.2898163426866567E-01, 
    0.5077671719424526E-02, 0.2462022242995652E-01, 
    0.3796424546988351E-02, 0.3310665553884920E-01, 
    0.3508140615002557E-01, 0.2950704004490954E-01, 
    0.8017748658106186E-02, 0.1975278275621657E-01, 
    0.2559496685146833E-01, 0.1349736154064142E-01, 
    0.4377789154192863E-02, 0.3252004086526811E-01, 
    0.2138516935436982E-01, 0.5015100019990189E-02, 
    0.1268265369185493E-01, 0.9980860281711927E-02, 
    0.4271671806230667E-01, 0.4274198226264675E-01, 
    0.2231057111860768E-01, 0.4634811584165991E-02, 
    0.3628925840326913E-01, 0.2963883283828189E-01, 
    0.3004053280506377E-01, 0.9656867842652003E-02, 
    0.1208027996271506E-01, 0.1067097545113014E-01, 
    0.2483163417257640E-01, 0.1052923525628832E-01, 
    0.2351898637367896E-01, 0.7064042544802277E-02, 
    0.1458074471394977E-01, 0.1255548203713306E-01, 
    0.2826856991152389E-01, 0.3484550434776074E-01, 
    0.3029815093584674E-01, 0.2852474545169975E-01, 
    0.9146009437075919E-02, 0.1990115803594483E-01, 
    0.2069177392023119E-01, 0.1613360621469677E-01, 
    0.1339344756044317E-01, 0.9274781087469082E-02, 
    0.1530433781694089E-01, 0.3917613624484167E-02, 
    0.4743890758461208E-01, 0.2001420743439166E-01, 
    0.2832855798485960E-02, 0.9111027190743106E-02, 
    0.3380728915063295E-02, 0.9357526881990973E-02, 
    0.2181586678118927E-01, 0.5497095242202148E-02, 
    0.6401916148758928E-02, 0.3861833303734160E-02, 
    0.8708380742935682E-02, 0.8616911070460571E-02, 
    0.5860715399693163E-02, 0.1134352478246400E-01, 
    0.2833705081651915E-01, 0.4870087478700985E-02, 
    0.1085506021481165E-01, 0.3085543406158568E-02, 
    0.2320131607643819E-02, 0.1217435492763373E-01, 
    0.1876868717373960E-02, 0.1764889584855365E-01, 
    0.6738482260511409E-02, 0.4000917701843592E-02, 
    0.2530717864603941E-02, 0.4454840968844280E-02, 
    0.2917862788830830E-02, 0.3031055423950871E-02, 
    0.1091305108499324E-02, 0.4268181631337478E-02, 
    0.6826137445231243E-02, 0.3247017230887967E-02, 
    0.3899933571265500E-02, 0.7578127425388515E-02 };

  r8mat_row_copy ( 3, n, 0, xs, x );
  r8mat_row_copy ( 3, n, 1, ys, x );
  r8mat_row_copy ( 3, n, 2, zs, x );
  r8vec_copy ( n, ws, w );

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

  fprintf ( stdout, "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}


