// File recommented by recomment.C
// on Jul  7 2014 at 20:37:06.
//
# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "cube_exactness.hpp"

int main ( );
void test01 ( );
void test02 ( );
void legendre_3d_set ( double a[], double b[], int nx, int ny, int nz, 
  double x[], double y[], double z[], double w[] );
void legendre_set ( int n, double x[], double w[] );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CUBE_EXACTNESS_PRB.
//
//  Discussion:
//
//    CUBE_EXACTNESS_PRB tests the CUBE_EXACTNESS library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "CUBE_EXACTNESS_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the CUBE_EXACTNESS library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CUBE_EXACTNESS_PRB\n";
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
//    TEST01 tests product Gauss-Legendre rules for the Legendre 3D integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[3];
  double b[3];
  int l;
  int n;
  int nx;
  int ny;
  int nz;
  int p_max;
  int t;
  double *w;
  double *x;
  double *y;
  double *z;

  a[0] = -1.0;
  a[1] = -1.0;
  a[2] = -1.0;
  b[0] = +1.0;
  b[1] = +1.0;
  b[2] = +1.0;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Product Gauss-Legendre rules for the 3D Legendre integral.\n";
  cout << "  Density function rho(x) = 1.\n";
  cout << "  Region: -1 <= x <= +1.\n";
  cout << "          -1 <= y <= +1.\n";
  cout << "          -1 <= z <= +1.\n";
  cout << "  Level: L\n";
  cout << "  Exactness: 2*L+1\n";
  cout << "  Order: N = (L+1)*(L+1)*(L+1)\n";

  for ( l = 0; l <= 5; l++ )
  {
    nx = l + 1;
    ny = l + 1;
    nz = l + 1;
    n = nx * ny * nz;
    t = 2 * l + 1;

    x = new double[n];
    y = new double[n];
    z = new double[n];
    w = new double[n];

    legendre_3d_set ( a, b, nx, ny, nz, x, y, z, w );

    p_max = t + 1;
    legendre_3d_exactness ( a, b, n, x, y, z, w, p_max );

    delete [] x;
    delete [] y;
    delete [] z;
    delete [] w;
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests product Gauss-Legendre rules for the Legendre 3D integral.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2014
//
//  Author:
//
//    John Burkardt
//
{
  double a[3];
  double b[3];
  int i;
  int l;
  int n;
  int nx;
  int ny;
  int nz;
  int p_max;
  int t;
  double *w;
  double *x;
  double *y;
  double *z;

  a[0] = -1.0;
  a[1] = -1.0;
  a[2] = -1.0;
  b[0] = +1.0;
  b[1] = +1.0;
  b[2] = +1.0;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Product Gauss-Legendre rules for the 3D Legendre integral.\n";
  cout << "  Density function rho(x) = 1.\n";
  cout << "  Region: -1 <= x <= +1.\n";
  cout << "          -1 <= y <= +1.\n";
  cout << "          -1 <= z <= +1.\n";
  cout << "  Exactness: 3 = 2 * min ( 2, 3, 4 ) - 1\n";
  cout << "  Order: N = 2 * 3 * 4\n";

  nx = 2;
  ny = 3;
  nz = 4;
  n = nx * ny * nz;

  x = new double[n];
  y = new double[n];
  z = new double[n];
  w = new double[n];

  legendre_3d_set ( a, b, nx, ny, nz, x, y, z, w );

  p_max = 4;

  legendre_3d_exactness ( a, b, n, x, y, z, w, p_max );

  delete [] x;
  delete [] y;
  delete [] z;
  delete [] w;

  return;
}
//****************************************************************************80

void legendre_3d_set ( double a[], double b[], int nx, int ny, int nz, 
  double x[], double y[], double z[], double w[] )

//****************************************************************************80
//
//  Purpose:
//
//    LEGENDRE_3D_SET: set a 3D Gauss-Legendre quadrature rule.
//
//  Discussion:
//
//    The integral:
//
//      integral ( z1 <= z <= z2 )
//               ( y1 <= y <= y2 ) 
//               ( x1 <= x <= x2 ) f(x,y,z) dx dy dz
//
//    The quadrature rule:
//
//      sum ( 1 <= i <= n ) w(i) * f ( x(i),y(i),z(i) )
//
//    where n = nx * ny * nz.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 August 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double A[3], B[3], the lower and upper integration
//    limits.
//
//    Input, int NX, NY, NZ, the orders in the X and Y directions.
//    These orders must be between 1 and 10.
//
//    Output, double X[N], Y[N], Z[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  int i;
  int j;
  int k;
  int l;
  double *wx;
  double *wy;
  double *wz;
  double *xx;
  double *yy;
  double *zz;
//
//  Get the rules for [-1,+1].
//
  xx = new double[nx];
  wx = new double[nx];
  legendre_set ( nx, xx, wx );

  yy = new double[ny];
  wy = new double[ny];
  legendre_set ( ny, yy, wy );

  zz = new double[nz];
  wz = new double[nz];
  legendre_set ( nz, zz, wz );
//
//  Adjust from [-1,+1] to [A,B].
//
  for ( i = 0; i < nx; i++ )
  {
    xx[i] = ( ( 1.0 - xx[i] ) * a[0]   
            + ( 1.0 + xx[i] ) * b[0] ) 
            /   2.0;
    wx[i] = wx[i] * ( b[0] - a[0] ) / 2.0;
  }

  for ( j = 0; j < ny; j++ )
  {
    yy[j] = ( ( 1.0 - yy[j] ) * a[1]   
            + ( 1.0 + yy[j] ) * b[1] ) 
            /   2.0;
    wy[j] = wy[j] * ( b[1] - a[1] ) / 2.0;
  }

  for ( k = 0; k < nz; k++ )
  {
    zz[k] = ( ( 1.0 - zz[k] ) * a[2]   
            + ( 1.0 + zz[k] ) * b[2] ) 
            /   2.0;
    wz[k] = wz[k] * ( b[2] - a[2] ) / 2.0;
  }
//
//  Compute the product rule.
//
  l = 0;
  for ( k = 0; k < nz; k++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      for ( i = 0; i < nx; i++ )
      {
        x[l] = xx[i];
        y[l] = yy[j];
        z[l] = zz[k];
        w[l] = wx[i] * wy[j] * wz[k];
        l = l + 1;
      }
    }
  }

  delete [] wx;
  delete [] wy;
  delete [] wz;
  delete [] xx;
  delete [] yy;
  delete [] zz;

  return;
}
//****************************************************************************80

void legendre_set ( int n, double x[], double w[] )

//****************************************************************************80
//
//  Purpose:
//  
//    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
//  
//  Discussion:
//    
//    The integral:
//  
//      Integral ( -1 <= X <= 1 ) F(X) dX
//  
//    Quadrature rule:
//  
//      Sum ( 1 <= I <= N ) W(I) * F ( X(I) )
//  
//    The quadrature rule is exact for polynomials through degree 2*N-1.
//  
//    The abscissas are the zeroes of the Legendre polynomial P(ORDER)(X).
//  
//    Mathematica can compute the abscissas and weights of a Gauss-Legendre
//    rule of order N for the interval [A,B] with P digits of precision
//    by the commands:
//
//    Needs["NumericalDifferentialEquationAnalysis`"]
//    GaussianQuadratureWeights [n, a, b, p ]
//  
//  Licensing:
//  
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//  
//    20 April 2010
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
//    Vladimir Krylov,
//    Approximate Calculation of Integrals,
//    Dover, 2006,
//    ISBN: 0486445798,
//    LC: QA311.K713.
//  
//    Arthur Stroud, Don Secrest,
//    Gaussian Quadrature Formulas,
//    Prentice Hall, 1966,
//    LC: QA299.4G3S7.
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
//    CRC Press, 1996,
//    ISBN: 0-8493-2479-3,
//    LC: QA47.M315.
//  
//  Parameters:
//  
//    Input, int N, the order.
//    N must be between 1 and 33 or 63/64/65, 127/128/129, 
//    255/256/257.
//
//    Output, double X[N], the abscissas.
//
//    Output, double W[N], the weights.
//
{
  if ( n == 1 )
  {
    x[0] = 0.000000000000000000000000000000;

    w[0] = 2.000000000000000000000000000000;
  }
  else if ( n == 2 )
  {
    x[0] = -0.577350269189625764509148780502;
    x[1] = 0.577350269189625764509148780502;

    w[0] = 1.000000000000000000000000000000;
    w[1] = 1.000000000000000000000000000000;
  }
  else if ( n == 3 )
  {
    x[0] = -0.774596669241483377035853079956;
    x[1] = 0.000000000000000000000000000000;
    x[2] = 0.774596669241483377035853079956;

    w[0] = 0.555555555555555555555555555556;
    w[1] = 0.888888888888888888888888888889;
    w[2] = 0.555555555555555555555555555556;
  }
  else if ( n == 4 )
  {
    x[0] = -0.861136311594052575223946488893;
    x[1] = -0.339981043584856264802665759103;
    x[2] = 0.339981043584856264802665759103;
    x[3] = 0.861136311594052575223946488893;

    w[0] = 0.347854845137453857373063949222;
    w[1] = 0.652145154862546142626936050778;
    w[2] = 0.652145154862546142626936050778;
    w[3] = 0.347854845137453857373063949222;
  }
  else if ( n == 5 )
  {
    x[0] = -0.906179845938663992797626878299;
    x[1] = -0.538469310105683091036314420700;
    x[2] = 0.000000000000000000000000000000;
    x[3] = 0.538469310105683091036314420700;
    x[4] = 0.906179845938663992797626878299;

    w[0] = 0.236926885056189087514264040720;
    w[1] = 0.478628670499366468041291514836;
    w[2] = 0.568888888888888888888888888889;
    w[3] = 0.478628670499366468041291514836;
    w[4] = 0.236926885056189087514264040720;
  }
  else if ( n == 6 )
  {
    x[0] = -0.932469514203152027812301554494;
    x[1] = -0.661209386466264513661399595020;
    x[2] = -0.238619186083196908630501721681;
    x[3] = 0.238619186083196908630501721681;
    x[4] = 0.661209386466264513661399595020;
    x[5] = 0.932469514203152027812301554494;

    w[0] = 0.171324492379170345040296142173;
    w[1] = 0.360761573048138607569833513838;
    w[2] = 0.467913934572691047389870343990;
    w[3] = 0.467913934572691047389870343990;
    w[4] = 0.360761573048138607569833513838;
    w[5] = 0.171324492379170345040296142173;
  }
  else if ( n == 7 )
  {
    x[0] = -0.949107912342758524526189684048;
    x[1] = -0.741531185599394439863864773281;
    x[2] = -0.405845151377397166906606412077;
    x[3] = 0.000000000000000000000000000000;
    x[4] = 0.405845151377397166906606412077;
    x[5] = 0.741531185599394439863864773281;
    x[6] = 0.949107912342758524526189684048;

    w[0] = 0.129484966168869693270611432679;
    w[1] = 0.279705391489276667901467771424;
    w[2] = 0.381830050505118944950369775489;
    w[3] = 0.417959183673469387755102040816;
    w[4] = 0.381830050505118944950369775489;
    w[5] = 0.279705391489276667901467771424;
    w[6] = 0.129484966168869693270611432679;
  }
  else if ( n == 8 )
  {
    x[0] = -0.960289856497536231683560868569;
    x[1] = -0.796666477413626739591553936476;
    x[2] = -0.525532409916328985817739049189;
    x[3] = -0.183434642495649804939476142360;
    x[4] = 0.183434642495649804939476142360;
    x[5] = 0.525532409916328985817739049189;
    x[6] = 0.796666477413626739591553936476;
    x[7] = 0.960289856497536231683560868569;

    w[0] = 0.101228536290376259152531354310;
    w[1] = 0.222381034453374470544355994426;
    w[2] = 0.313706645877887287337962201987;
    w[3] = 0.362683783378361982965150449277;
    w[4] = 0.362683783378361982965150449277;
    w[5] = 0.313706645877887287337962201987;
    w[6] = 0.222381034453374470544355994426;
    w[7] = 0.101228536290376259152531354310;
  }
  else if ( n == 9 )
  {
    x[0] = -0.968160239507626089835576203;
    x[1] = -0.836031107326635794299429788;
    x[2] = -0.613371432700590397308702039;
    x[3] = -0.324253423403808929038538015;
    x[4] = 0.000000000000000000000000000;
    x[5] = 0.324253423403808929038538015;
    x[6] = 0.613371432700590397308702039;
    x[7] = 0.836031107326635794299429788;
    x[8] = 0.968160239507626089835576203;

    w[0] = 0.081274388361574411971892158111;
    w[1] = 0.18064816069485740405847203124;
    w[2] = 0.26061069640293546231874286942;
    w[3] = 0.31234707704000284006863040658;
    w[4] = 0.33023935500125976316452506929;
    w[5] = 0.31234707704000284006863040658;
    w[6] = 0.26061069640293546231874286942;
    w[7] = 0.18064816069485740405847203124;
    w[8] = 0.081274388361574411971892158111;
  }
  else if ( n == 10 )
  {
    x[0] = -0.973906528517171720077964012;
    x[1] = -0.865063366688984510732096688;
    x[2] = -0.679409568299024406234327365;
    x[3] = -0.433395394129247190799265943;
    x[4] = -0.148874338981631210884826001;
    x[5] = 0.148874338981631210884826001;
    x[6] = 0.433395394129247190799265943;
    x[7] = 0.679409568299024406234327365;
    x[8] = 0.865063366688984510732096688;
    x[9] = 0.973906528517171720077964012;

    w[0] = 0.066671344308688137593568809893;
    w[1] = 0.14945134915058059314577633966;
    w[2] = 0.21908636251598204399553493423;
    w[3] = 0.26926671930999635509122692157;
    w[4] = 0.29552422471475287017389299465;
    w[5] = 0.29552422471475287017389299465;
    w[6] = 0.26926671930999635509122692157;
    w[7] = 0.21908636251598204399553493423;
    w[8] = 0.14945134915058059314577633966;
    w[9] = 0.066671344308688137593568809893;
  }
  else
  {
    cerr << "\n";
    cerr << "LEGENDRE_SET - Fatal error\n";
    cerr << "  Illegal value of N = " << n << "\n";
    cerr << "  Legal values are 1:10.\n";
    exit ( 1 );
  }
  return;
}