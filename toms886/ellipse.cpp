# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cmath>
# include <cstring>

using namespace std;

# include "toms886.hpp"

int main ( );
double sigma1 ( double t1, double t2, double c1, double c2, double alpha, 
  double beta );
double isigm1 ( double sigma1, double sigma2, double c1, double c2, 
  double alpha, double beta );
double sigma2 ( double t1, double t2, double c1, double c2, double alpha, 
  double beta );
double isigm2 ( double sigma1, double sigma2, double c1, double c2, 
  double alpha, double beta );
double phi ( double x );
double iphi ( double x );
void target ( double c1, double c2, double alpha, double beta, int ntg1, 
  int ntgmax, double tg1[], double tg2[], int &ntg );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for ELLIPSE.
//
//  Discussion:
//
//    This driver computes the interpolation of the Franke function
//    on the ellipse E((C1,C2),ALPHA,BETA) = E((0.5,0.5),0.5,0.5)  
//    at the first family of Padua points. 
//
//    The ellipse has the equation:
//
//      ( ( X - C1 ) / ALPHA )^2 + ( ( Y - C2 ) / BETA )^2 = 1
//
//    The degree of interpolation DEG = 60 and the number of target 
//    points is NTG = NTG1 ^ 2 - 2 * NTG1 + 2, NTG1 = 100.  
//
//    The maps from the reference square [-1,1]^2 to the current domain 
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
//    This C++ version by John Burkardt.
//
//  Reference:
//
//    Marco Caliari, Stefano de Marchi, Marco Vianello,
//    Algorithm 886:
//    Padua2D: Lagrange Interpolation at Padua Points on Bivariate Domains,
//    ACM Transactions on Mathematical Software,
//    Volume 35, Number 3, October 2008, Article 21, 11 pages.
//
//  Local Parameters:
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
//    Local, int NTG1, the parameter determining the number of target points.
//
//    Local, int NPD, the number of Padua points = (DEG + 1) * (DEG + 2) / 2.
//
//    Local, int NTG, the number of target points, dependent on NTG1.
//
//    Local, double PD1[NPDMAX], the first coordinates of 
//    the Padua points.
//
//    Local, double PD2[NPDMAX], the second coordinates of the 
//    Padua points.
//
//    Local, double WPD[NPDMAX], the weights.
//
//    Local, double FPD[NPDMAX], the function at the Padua points.
//
//    Workspace, double RAUX1[(DEGMAX+1)*(DEGMAX+2)].
//
//    Workspace, double RAUX2[(DEGMAX+1)*(DEGMAX+2)].
//
//    Local, double C0[(0:DEGMAX+1)*(0:DEGMAX+1)], the coefficient matrix.
//
//    Local, double TG1[NTGMAX], the first coordinates of the 
//    target points.
//
//    Local, double TG2[NTGMAX], the second coordinates of the 
//    target points.
//
//    Local, double INTFTG[NTGMAX], the values of the 
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
# define NTGMAX ( NTG1MX * NTG1MX - 2 * NTG1MX + 2 )

  double alpha;
  double beta;
  double c0[(DEGMAX+2)*(DEGMAX+2)];
  double c1;
  double c2;
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
  double wpd[NPDMAX];
  double x;
  double y;

  alpha = 0.5;
  beta = 0.5;
  c1 = 0.5;
  c2 = 0.5;
  family = 1;
  deg = 60;
  ntg1 = 100;

  timestamp ( );
  cout << "\n";
  cout << "ELLIPSE:\n";
  cout << "  C++ version\n";
  cout << "  Interpolation of the Franke function\n";
  cout << "  on the disk with center = (0.5,0.5) and radius = 0.5\n";
  cout << "  of degree = " << deg << "\n";

  if ( degmax < deg )
  {
    cerr << "\n";
    cerr << "ELLIPSE - Fatal error!\n";
    cerr << "  DEGMAX < DEG.\n";
    cerr << "  DEG =    " << deg << "\n";
    cerr << "  DEGMAX = " << degmax << "\n";
    exit ( 1 );
  }
//   
//  Build the first family of Padua points in the square [-1,1]^2.
//
  pdpts ( deg, pd1, pd2, wpd, npd );
//    
//  Compute the Franke function at Padua points mapped to the region.
//
  for ( i = 0; i < npd; i++ )
  {
    x = sigma1 ( pd1[i], pd2[i], c1, c2, alpha, beta );
    y = sigma2 ( pd1[i], pd2[i], c1, c2, alpha, beta );
    fpd[i] = franke ( x, y );
  }
//
//  Write X, Y, F(X,Y) to a file.
//
  filename = "ellipse_fpd.txt";
  output.open ( filename.c_str ( ) );
  for ( i = 0; i < npd; i++ )
  {
    x = sigma1 ( pd1[i], pd2[i], c1, c2, alpha, beta );
    y = sigma2 ( pd1[i], pd2[i], c1, c2, alpha, beta );
    output << x
           << "  " << y
           << "  " << fpd[i] << "\n";
  }
  output.close ( );
  cout << "\n";
  cout << "  Wrote F(x,y) at Padua points in '" << filename << "'\n";
//
//  Compute the matrix C0 of the coefficients in the bivariate
//  orthonormal Chebyshev basis.
//
  padua2 ( deg, degmax, npd, wpd, fpd, raux1, raux2, c0, esterr );
//    
//  Evaluate the target points in the region.
//
  target ( c1, c2, alpha, beta, ntg1, ntgmax, tg1, tg2, ntg );
//
//  Evaluate the interpolant at the target points.
//
  for ( i = 0; i < ntg; i++ ) 
  {
    x = isigm1 ( tg1[i], tg2[i], c1, c2, alpha, beta );
    y = isigm2 ( tg1[i], tg2[i], c1, c2, alpha, beta );
    intftg[i] = pd2val ( deg, degmax, c0, x, y );
  }
//
//  Write the function value at target points to a file.
//
  filename = "ellipse_ftg.txt";
  output.open ( filename.c_str ( ) );
  for ( i = 0; i < ntg; i++ )
  {
    output <<  tg1[i]
           << "  " << tg2[i]
           << "  " << franke ( tg1[i], tg2[i] ) << "\n";
  }
  output.close ( );
  cout << "  Wrote F(x,y) at target points in '" << filename << "'\n";
//
//  Write the interpolated function value at target points to a file.
//
  filename = "ellipse_itg.txt";
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
  cout << "  Expected error:   " << 0.1769E-09 << "\n";
//
//  Terminate.
//
  cout << "\n";
  cout << "ELLIPSE:\n";
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

double sigma1 ( double t1, double t2, double c1, double c2, double alpha, 
  double beta )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA1 maps first coordinate from square to ellipse.
//
//  Discussion:
//
//    This function returns the first component of the map 
//    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
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
//    This C++ version by John Burkardt.
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
//    Input, double C1, C2, ALPHA, BETA, the center and scale
//    parameters of the ellipse.
//
//    Output, double SIGMA1, the X coordinate of the corresponding
//    point in the ellipse.
//
{
  double value;

  value = c1 - alpha * t2 * sin ( phi ( t1 ) );

  return value;
}
//****************************************************************************80

double isigm1 ( double sigma1, double sigma2, double c1, double c2, 
  double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    ISIGM1 maps the first coordinate from the ellipse to the square.
//
//  Discussion:
//
//    This function returns the first component of the map 
//    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//  
//  Modified:
//
//    09 February 2014
//
//  Author:
//
//    Original FORTRAN77 version by Marco Caliari, Stefano De Marchi, 
//    Marco Vianello.
//    This C++ version by John Burkardt.
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
//    Input, double C1, C2, ALPHA, BETA, the center and scale
//    parameters of the ellipse.
//
//    Output, double ISIGM1, the X coordinate of the corresponding
//    point in the square.
//
{
  double value;

  if ( sigma2 == c2 )
  {
    value = 1.0;
  }
  else
  {
    value = iphi ( atan ( beta * ( c1 - sigma1 ) / 
      ( alpha * ( sigma2 - c2 ) ) ) );
  }   

  return value;
}
//****************************************************************************80

double sigma2 ( double t1, double t2, double c1, double c2, double alpha, 
  double beta )

//****************************************************************************80
//
//  Purpose:
//
//    SIGMA2 maps the second coordinate from square to ellipse.
//
//  Discussion:
//
//    This function returns the second component of the map 
//    from the square [-1,1]^2 to the ellipse E((C1,C2),ALPHA,BETA).
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
//    This C++ version by John Burkardt.
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
//    Input, double C1, C2, ALPHA, BETA, the center and scale
//    parameters of the ellipse.
//
//    Output, double SIGMA2, the Y coordinate of the corresponding
//    point in the ellipse.
//
{
  double value;

  value = c2 + beta * t2 * cos ( phi ( t1 ) );

  return value;
}
//****************************************************************************80

double isigm2 ( double sigma1, double sigma2, double c1, double c2, 
  double alpha, double beta )

//****************************************************************************80
//
//  Purpose:
//
//    ISIGM2 maps second coordinate from ellipse to the square.
//
//  Discussion:
//
//    This function returns the second component of the map 
//    from the ellipse E((C1,C2),ALPHA,BETA) to the square [-1,1]^2. 
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
//    This C++ version by John Burkardt.
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
//    Input, double C1, C2, ALPHA, BETA, the center and scale
//    parameters of the ellipse.
//
//    Output, double ISIGM2, the Y coordinate of the corresponding
//    point in the square.
//
{
  double value;

  if ( sigma2 == c2 )
  {
    value = ( c1 - sigma1 ) / alpha;
  }
  else
  {
    value = sqrt ( beta * beta * pow ( c1 - sigma1, 2 ) + 
                   alpha * alpha * pow ( c2 - sigma2, 2 ) ) 
      / ( alpha * beta ) * r8_sign ( sigma2 - c2 );
  }   

  return value;
}
//****************************************************************************80

double phi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    PHI maps from [-1,+1] to [-pi/2,+pi/2].
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
//    This C++ version by John Burkardt.
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
//    Input, double X, a point in [-1,+1];
//
//    Output, double PHI, a corresponding point in [-pi/2,+pi/2].
//
{
  const double pi = 3.1415926535897931;
  double value;

  value = pi * x / 2.0;

  return value;
}
//****************************************************************************80

double iphi ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    IPHI maps from [-pi/2,+pi/2] to [-1,+1].
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
//    This C++ version by John Burkardt.
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
//    Input, double X, a point in [-pi/2,+pi/2].
//
//    Output, double IPHI, a corresponding point in [-1,+1].
//
{
  const double pi = 3.1415926535897931;
  double value;

  value = 2.0 * x / pi;

  return value;
}
//****************************************************************************80

void target ( double c1, double c2, double alpha, double beta, int ntg1, 
  int ntgmax, double tg1[], double tg2[], int &ntg )

//****************************************************************************80
//
//  Purpose:
//
//    TARGET returns the target points on the ellipse.
//
//  Discussion:
//
//    Target points on the ellipse E((C1,C2),ALPHA,BETA).
//    The number of target points is NTG = NTG1^2 - 2 * NTG1 + 2.
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
//    This C++ version by John Burkardt.
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
//    Input, double C1, C2, ALPHA, BETA, the center and scale
//    parameters of the ellipse.
//
//    Input, int NTG1, a parameter determining the number 
//    of target points.  2 <= NTG1.
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
  double t;

  if ( ntg1 < 2 )
  {
    cout << "\n";
    cout << "TARGET - Fatal error!\n";
    cout << "  NTG1 < 2\n";
    cout << "  NTG1 = " << ntg1 << "\n";
    exit ( 1 );
  }

  if ( ntgmax < ntg1 * ntg1 - 2 * ntg1 + 2 )
  {
    cout << "\n";
    cout << "TARGET - Fatal error!\n";
    cout << "  NTGMAX < NTG1 * NTG1 - 2 * NTG1 + 2.\n";
    cout << "  NTG1 = " << ntg1 << "\n";
    cout << "  NTGMAX = " << ntgmax << "\n";
    exit ( 1 );
  }      

  i = 1;
  j = 1;
  ntg = 0;

  tg1[ntg] = alpha * ( - 1.0 + ( double ) ( i - 1 ) * 2.0 
    / ( double ) ( ntg1 - 1 ) ) + c1;

  t = - 1.0 + ( double ) ( i - 1 ) * 2.0 / ( double ) ( ntg1 - 1 );

  tg2[ntg] =  beta * ( - 1.0 + ( double ) ( j - 1 ) * 2.0 
    / ( double ) ( ntg1 - 1 ) ) * sqrt ( 1.0 - t * t ) + c2;

  ntg = ntg + 1;

  for ( i = 2; i <= ntg1 - 1; i++ )
  {
    for ( j = 1; j <= ntg1; j++ )
    {
      tg1[ntg] = alpha * ( - 1.0 + ( double ) ( i - 1 ) * 2.0 
        / ( double ) ( ntg1 - 1 ) ) + c1;

      t = - 1.0 + ( double ) ( i - 1 ) * 2.0 / ( double ) ( ntg1 - 1 );

      tg2[ntg] =  beta * ( - 1.0 + ( double ) ( j - 1 ) * 2.0 
        / ( double ) ( ntg1 - 1 ) ) * sqrt ( 1.0 - t * t ) + c2;

      ntg = ntg + 1;
    }
  }

  i = ntg1;
  j = 1;

  tg1[ntg] = alpha * ( -1.0 + ( double ) ( i - 1 ) * 2.0 
    / ( double ) ( ntg1 - 1 ) ) + c1;

  t = - 1.0 + ( double ) ( i - 1 ) * 2.0 / ( double ) ( ntg1 - 1 );

  tg2[ntg] = beta * ( -1.0 + ( double ) ( j - 1 ) * 2.0 
    / ( double ) ( ntg1 - 1 ) ) * sqrt ( 1.0 - t * t ) + c2;

  ntg = ntg + 1;

  return;
}
