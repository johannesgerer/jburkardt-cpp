# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>

using namespace std;

# include "toms886.hpp"

int main ( );
double sigma1 ( double t1, double t2, double u1, double u2, double v1, 
  double v2, double w1, double w2 );
double isigm1 ( double sigma1, double sigma2, double u1, double u2, double v1, 
  double v2, double w1, double w2 );
double sigma2 ( double t1, double t2, double u1, double u2, double v1, 
  double v2, double w1, double w2 );
double isigm2 ( double sigma1, double sigma2, double u1, double u2, double v1, 
  double v2, double w1, double w2 );
void target ( double u1, double u2, double v1, double v2, double w1, double w2, 
  int ntg1, int ntgmax, double tg1[], double tg2[], int &ntg );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TRIANGLE.
//
//  Discussion:
//
//    This driver computes the interpolation of the Franke function
//    on the triangle T(U,V,W) with vertices U=(U1,U2)=(0,0), 
//    V=(V1,V2)=(1,0) and W=(W1,W2)=(0,1) (unit triangle) 
//    at the first family of Padua points. 
//
//    The degree of interpolation is DEG = 60 and the number of target 
//    points is NTG = NTG1 ** 2 - NTG1 + 1, NTG1 = 100.
//
//    The maps from the reference square [-1,1]^2 to the triangle
//    are SIGMA1 and SIGMA2 with inverses ISIGM1 and ISIGM2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Local, int DEGMAX, the maximum degree of interpolation.
//
//    Local, int NPDMAX, the maximum number of Padua points
//    = (DEGMAX + 1) * (DEGMAX + 2) / 2.
//
//    Local, int NTG1MX, the maximum value of the parameter determining 
//    the number of target points.
//
//    Local, int NTGMAX, the maximum number of target points,
//    dependent on NTG1MX.
//
//    Local, int DEG, the degree of interpolation.
//
//    Local, int NTG1, the parameter determining the number 
//   of target points.
//
//    Local, int NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
//
//    Local, int NTG, the number of target points, dependent on NTG1.
//
//    Local, double PD1(NPDMAX), the first coordinates of 
//    the Padua points.
//
//    Local, double PD2(NPDMAX), the second coordinates of the 
//    Padua points.
//
//    Local, double WPD(NPDMAX), the weights.
//
//    Local, double FPD(NPDMAX), the function at the Padua points.
//
//    Workspace, double RAUX1(DEGMAX+1)*(DEGMAX+2)).
//
//    Workspace, double RAUX2(DEGMAX+1)*(DEGMAX+2)).
//
//    Local, double C0(0:DEGMAX+1,0:DEGMAX+1), the coefficient matrix.
//
//    Local, double TG1(NTGMAX), the first coordinates of the 
//    target points.
//
//    Local, double TG2(NTGMAX), the second coordinates of the 
//    target points.
//
//    Local, double INTFTG(NTGMAX), the values of the 
//    interpolated function.
//
//    Local, double MAXERR, the maximum norm of the error at target 
//    points.
//
//    Local, double ESTERR, the estimated error.
//
{
# define DEGMAX 60
# define NTG1MX 100
# define NPDMAX ( ( DEGMAX + 1 ) * ( DEGMAX + 2 ) / 2 )
# define NTGMAX ( NTG1MX * NTG1MX - NTG1MX + 1 )

  double c0[(DEGMAX+2)*(DEGMAX+2)];
  int deg;
  int degmax = DEGMAX;
  double esterr;
  int family;
  string filename;
  double fmax;
  double fmin;
  double fpd[NPDMAX];
  double fxy;
  int i;
  double intftg[NTGMAX];
  double ixy;
  double maxdev;
  double maxerr;
  double mean;
  int npd;
  int npdmax = NPDMAX;
  int ntg;
  int ntg1;
  int ntg1mx = NTG1MX;
  int ntgmax = NTGMAX;
  ofstream output;
  double pd1[NPDMAX];
  double pd2[NPDMAX];
  double raux1[(DEGMAX+1)*(DEGMAX+2)];
  double raux2[(DEGMAX+1)*(DEGMAX+2)];
  double tg1[NTGMAX];
  double tg2[NTGMAX];
  double u1;
  double u2;
  double v1;
  double v2;
  double w1;
  double w2;
  double wpd[NPDMAX];
  double x;
  double y;

  u1 = 0.0;
  u2 = 0.0;
  v1 = 1.0;
  v2 = 0.0;
  w1 = 0.0;
  w2 = 1.0;
  family = 1;
  deg = 60;
  ntg1 = 100;

  timestamp ( );
  cout << "\n";
  cout << "TRIANGLE:\n";
  cout << "  C version\n";
  cout << "  Interpolation of the Franke function\n";
  cout << "  on the unit triangle T((0,0),(1,0),(0,1))\n";
  cout << "  at degree = " << deg << "\n";

  if ( degmax < deg )
  {
    cerr << "\n";
    cerr << "TRIANGLE - Fatal error!\n";
    cerr << "  DEGMAX < DEG.\n";
    cerr << "  DEG =    " << deg << "\n";
    cerr << "  DEGMAX = " << degmax << "\n";
    exit ( 1 );
  }
//
//  Build the first family of Padua points in the square [-1,1]^2
// 
  pdpts ( deg, pd1, pd2, wpd, npd );
// 
//  Compute the Franke function at Padua points mapped to T(U,V,W).
//
  for ( i = 0; i < npd; i++ )
  {
    x = sigma1 ( pd1[i], pd2[i], u1, u2, v1, v2, w1, w2 );
    y = sigma2 ( pd1[i], pd2[i], u1, u2, v1, v2, w1, w2 );
    fpd[i] = franke ( x, y );
  }
//
//  Write X, Y, F(X,Y) to a file.
//
  filename = "triangle_fpd.txt";
  output.open ( filename.c_str ( ) );
  for ( i = 0; i < npd; i++ )
  {
    x = sigma1 ( pd1[i], pd2[i], u1, u2, v1, v2, w1, w2 );
    y = sigma2 ( pd1[i], pd2[i], u1, u2, v1, v2, w1, w2 );
    output << x
           << "  " << y
           << "  " << fpd[i] << "\n";
  }
  output.close ( );
  cout << "\n";
  cout << "  Wrote F(x,y) at Padua points in '" << filename << "'\n";
//    
//  Compute the matrix C0 of the coefficients in the bivariate
//  orthonormal Chebyshev basis
//
  padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr );
//    
//  Build the set of target points on T(U,V,W)
//
  target ( u1, u2, v1, v2, w1, w2, ntg1, ntgmax, tg1, tg2, ntg );
//    
//  Evaluate the interpolant at the target points.
//
  for ( i = 0; i < ntg; i++ )
  {
    x = isigm1 ( tg1[i], tg2[i], u1, u2, v1, v2, w1, w2 );
    y = isigm2 ( tg1[i], tg2[i], u1, u2, v1, v2, w1, w2 );
    intftg[i] = pd2val ( deg, degmax, c0, x, y );
  }
//
//  Write the function value at target points to a file.
//
  filename = "triangle_ftg.txt";
  output.open ( filename.c_str ( ) );
  for ( i = 0; i < ntg; i++ )
  {
    output << tg1[i]
           << "  " << tg2[i]
           << "  " << franke ( tg1[i], tg2[i] ) << "\n";
  }
  output.close ( );
  cout << "  Wrote F(x,y) at target points in '" << filename << "'\n";
//
//  Write the interpolated function value at target points to a file.
//
  filename = "triangle_itg.txt";
  output.open ( filename.c_str ( ) );
  for ( i = 0; i < ntg; i++ )
  {
    output << tg1[i]
           << "  " << tg2[i]
           << "  " << intftg[i] << "\n";
  }
  output.close ( );
  cout << "  Wrote I(F)(x,y) at target points in '" << filename << "'\n";
//
//  Compute the error relative to the max deviation from the mean.
//
  maxerr = 0.0;
  mean = 0.0;
  fmax = - r8_huge ( );
  fmin = + r8_huge ( );

  for ( i = 0; i < ntg; i++ )
  {
    fxy = franke ( tg1[i], tg2[i] );
    ixy = intftg[i];
    maxerr = r8_max ( maxerr, fabs ( fxy - ixy ) );
    mean = mean + fxy;
    fmax = r8_max ( fmax, fxy );
    fmin = r8_min ( fmin, fxy );
  }
 
  if ( fmax == fmin )
  {
    maxdev = 1.0;
  }
  else
  {
    mean = mean / ( double ) ( ntg );
    maxdev = r8_max ( fmax - mean, mean - fmin );
  }
//
//  Print error ratios.
//
  cout << "\n";
  cout << "  Estimated error:  " << esterr / maxdev << "\n";
  cout << "  Actual error:     " << maxerr / maxdev << "\n";
  cout << "  Expected error:   " << 0.1226E-09 << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "TRIANGLE:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;

# undef DEGMAX
# undef NTG1MX
# undef NPDMAX
# undef NTGMAX
}
//****************************************************************************80

double sigma1 ( double t1, double t2, double u1, double u2, double v1, 
  double v2, double w1, double w2 )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA1 maps first coordinate from square to triangle.
//
//  Discussion:
//
//    This function returns the first component of the map
//    from the square [-1,1]^2 to the triangle T(U,V,W). 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, double T1, T2, the coordinates of a point in the square.
//
//    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
//    vertices of the triangle.
//
//    Output, double SIGMA1, the X coordinate of the corresponding
//    point in the triangle.
//
{
  double value;

  value = ( v1 - u1 ) * ( 1.0 + t1 ) 
    * ( 1.0 - t2 ) / 4.0 
    + ( w1 - u1 ) * ( 1.0 + t2 ) / 2.0 + u1;

  return value;
}
//****************************************************************************80

double isigm1 ( double sigma1, double sigma2, double u1, double u2, double v1, 
  double v2, double w1, double w2 )

//****************************************************************************80
//
//  Purpose:
//
//    ISIGM1 maps first coordinate from triangle to the square.
//
//  Discussion:
//
//    This functions returns the first component of the map
//    from the the triangle T(U,V,W) to the square [-1,1]^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, double SIGMA1, SIGMA2, the coordinates of a point 
//    in the triangle.
//
//    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
//    vertices of the triangle.
//
//    Output, double ISIGM1, the X coordinate of the corresponding
//    point in the square.
//
{
  double rho1;
  double rho2;
  double value;

  rho1 = ( sigma1 * ( w2 - u2 ) - sigma2 * ( w1 - u1 ) 
    + ( w1 - u1 ) * u2 - ( w2 - u2 ) * u1 ) / 
    ( ( v1 - u1 ) * ( w2 - u2 ) - ( v2 - u2 ) * ( w1 - u1 ) );

  rho2 = ( sigma1 * ( v2 - u2 ) - sigma2 * ( v1 - u1 ) 
    + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / 
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) );

  if ( rho2 == 1.0 )
  {
    value = 0.0;
  }
  else
  {
    value = 2.0 * rho1 / ( 1.0 - rho2 ) - 1.0;
  }  
 
  return value;
}
//****************************************************************************80

double sigma2 ( double t1, double t2, double u1, double u2, double v1, 
  double v2, double w1, double w2 )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA2 maps the second coordinate from square to triangle.
//
//  Discussion:
//
//    This functions returns the second component of the map
//    from the square [-1,1]^2 to the triangle T(U,V,W).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, double T1, T2, the coordinates of a point in the square.
//
//    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
//    vertices of the triangle.
//
//    Output, double SIGMA2, the Y coordinate of the corresponding
//    point in the triangle.
//
{
  double value;

  value = ( v2 - u2 ) * ( 1.0 + t1 ) 
    * ( 1.0 - t2 ) / 4.0 + ( w2 - u2 ) 
    * ( 1.0 + t2 ) / 2.0 + u2;

  return value;
}
//****************************************************************************80

double isigm2 ( double sigma1, double sigma2, double u1, double u2, double v1, 
  double v2, double w1, double w2 )

//****************************************************************************80
//
//  Purpose:
//
//    ISIGM2 maps second coordinate from triangle to the square.
//
//  Discussion:
//
//    This functions returns the second component of the map
//    from the the triangle T(U,V,W) to the square [-1,1]^2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, double SIGMA1, SIGMA2, the coordinates of a point 
//    in the ellipse.
//
//    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
//    vertices of the triangle.
//
//    Output, double ISIGM2, the Y coordinate of the corresponding
//    point in the triangle.
//
{
  double rho2;
  double value;

  rho2 = ( sigma1 * ( v2 - u2 ) - 
    sigma2 * ( v1 - u1) + ( v1 - u1 ) * u2 - ( v2 - u2 ) * u1 ) / 
    ( ( w1 - u1 ) * ( v2 - u2 ) - ( w2 - u2 ) * ( v1 - u1 ) );

  if ( rho2 == 1.0 )
  {
    value = 1.0;
  }
  else
  {
    value = 2.0 * rho2 - 1.0;
  } 
  
  return value;
}
//****************************************************************************80

void target ( double u1, double u2, double v1, double v2, double w1, double w2, 
  int ntg1, int ntgmax, double tg1[], double tg2[], int &ntg )

//****************************************************************************80
//
//  Purpose:
//
//    TARGET returns the target points on the triangle.
//
//  Discussion:
//
//    Target points on the triangle T(U,V,W).
//    The number of target points is NTG = NTG1^2 - NGT1 + 1.
//
//  Licensing:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C version by John Burkardt.
//  
//  Modified:
//
//    16 February 2014
//
//  Author:
//
//    Marco Caliari, Stefano De Marchi, Marco Vianello
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Parameters:
//
//    Input, double U1, U2, V1, V2, W1, W2, the coordinates of the 
//    vertices of the triangle.
//
//    Input, int NTG1, a parameter determining the number 
//    of target points
//
//    Input, int NTGMAX, the maximum number of target points.
//
//    Output, double TG1[NTG], TG2[NTG], the X and Y coordinates
//    of the target points.
//
//    Output, int &NTG, the number of target points computed.
//
{
  int i;
  int j;

  if ( ntg1 < 2 )
  {
    cout << "\n";
    cout << "TARGET - Fatal error!\n";
    cout << "  NTG1 < 2\n";
    cout << "  NTG1 = " << ntg1 << "\n";
    exit ( 1 );
  }  
   
  if ( ntgmax < ntg1 * ntg1- ntg1 + 1 )
  {
    cout << "\n";
    cout << "TARGET - Fatal error!\n";
    cout << "  NTGMAX < NTG1 * NTG1 - NTG1 + 1.\n";
    cout << "  NTG1 = " << ntg1 << "\n";
    cout << "  NTGMAX = " << ntgmax << "\n";
    exit ( 1 );
  }    

  ntg = 0;

  for ( i = 1; i <= ntg1 - 1; i++ )
  {
    for ( j = 1; j <= ntg1; j++ )
    {
      tg1[ntg] = ( v1 - u1 ) * ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) 
        + ( w1 - u1 ) * ( ( double ) ( j - 1 ) / ( double ) ( ntg1 - 1 ) ) 
        * ( 1.0 - ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) ) + u1;

      tg2[ntg] = ( v2 - u2 ) * ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) 
        + ( w2 - u2 ) * ( ( double ) ( j - 1 ) / ( double ) ( ntg1 - 1 ) ) 
        * ( 1.0 - ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) ) + u2;

      ntg = ntg + 1;
    }
  }

  i = ntg1;
  j = 1;

  tg1[ntg] = ( v1 - u1 ) * ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) 
    + ( w1 - u1 ) * ( ( double ) ( j - 1 ) / ( double ) ( ntg1 - 1 ) ) 
    * ( 1.0 - ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) ) + u1;

  tg2[ntg] = ( v2 - u2 ) * ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) 
    + ( w2 - u2 ) * ( ( double ) ( j - 1 ) / ( double ) ( ntg1 - 1 ) ) 
    * ( 1.0 - ( double ) ( i - 1 ) / ( double ) ( ntg1 - 1 ) ) + u2;

  ntg = ntg + 1;

  return;
}
