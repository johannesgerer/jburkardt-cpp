# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

# include "TasmanianSparseGrid.hpp"

using namespace std;

int main ( int argc, const char ** argv );
void tsg_test01 ( );
void tsg_test02 ( );
void tsg_test03 ( );
void tsg_test04 ( );
void tsg_test05 ( );
void tsg_test06 ( );
void timestamp ( );

//****************************************************************************80

int main ( int argc, const char ** argv )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for TSG_PRB.
//
//  Discussion:
//
//    TSG_PRB tests the TSG library.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{
  timestamp ( );
  cout << endl;
  cout << "TSG_PRB:" << endl;
  cout << "  C++ version" << endl;
  cout << "  Test the TSG library." << endl;
//
//  Call the individual test functions.
//
  tsg_test01 ( );
  tsg_test02 ( );
  tsg_test03 ( );
  tsg_test04 ( );
  tsg_test05 ( );
  tsg_test06 ( );
//
//  Terminate.
//
  cout << endl;
  cout << "TSG_PRB:" << endl;
  cout << "  Normal end of execution." << endl;
  cout << endl;
  timestamp ( );

  return 0;
}
//****************************************************************************80

void tsg_test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST01 estimates an integral using a CC Smolyak grid.
//
//  Discussion:
//
//    Make a classical Smolyak grid using Clenshaw-Curtis quadrature.
//    Integrate the function f(x,y) = exp ( -x^2 ) * cos ( y ) over [-1,+1]^2.
//
//    The exact answer is 2.513723354063905e+00.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{    
  cout << endl;
  cout << "TSG_TEST01" << endl;
  cout << "  Make a classical 2D Smolyak grid using the Clenshaw-Curtis rule." 
       << endl;
  cout << "  Integrate f(x,y) = exp ( -x^2 ) * cos ( y ) over [-1,+1]^2"
       << endl;
//
//  Create the grid.
//
  const int dimension = 2;
  const int outputs = 0;
  const int level = 7;
  TasGrid::TasmanianSparseGrid grid;

  cout << "  The grid is of level " << level << endl;

  grid.makeGlobalGrid ( dimension, outputs, level, TasGrid::type_level,
    TasGrid::rule_clenshawcurtis );

  cout << "  The grid uses " << grid.getNumPoints() << " points." << endl;
//
//  Get the points and weights of the grid.
//              
  double *points = 0;
  double *weights = 0;

  grid.getPoints ( points );
  grid.getWeights ( weights );
//
//  Estimate the integral.
//
  double sum = 0.0;
  double x;
  double y;

  for ( int i = 0; i < grid.getNumPoints(); i++ )
  {
    x = points[2*i];
    y = points[2*i+1];
    sum += weights[i] * exp ( - x * x ) * cos ( y );
  }
//
//  Report.
//
  double exact = 2.513723354063905e+00;

  cout.precision(17);
  cout << "  The estimated integral is: " << std::scientific << sum << endl;
  cout << "  The exact integral is:     " << std::scientific << exact << endl;
  cout.precision(4);
  cout << "  The error is:              " << std::scientific 
       << fabs ( sum - exact ) << endl;
//
//  Free memory.
//
  delete [] points;
  delete [] weights;

  return;
}
//****************************************************************************80

void tsg_test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST02 uses a CC Smolyak grid for interpolation and integration.
//
//  Discussion:
//
//    Make a classical Smolyak grid using the Clenshaw-Curtis rule.
//    The grid should exactly interpolate polynomials of total degree 10
//    or less.
//    Compare f(x,y) to its interpolant.
//    Compare the integral of f(x,y) to the estimated integral.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{
  cout << endl;
  cout << "TSG_TEST02" << endl;
  cout << "  Make a 2D Smolyak grid using the Clenshaw-Curtis rule" << endl;
  cout << "  over the region [-1,+1]^2." << endl;
  cout << "  This grid interpolates polynomials of total degree 10 or less." 
       << endl;
  cout << "  Consider f(x,y) = exp ( -x^2 ) * cos ( y )" << endl;
  cout << "  Compare f(x,y) to the interpolant I(f)(x,y)." << endl;
  cout << "  Compare the integral of f(x,y) to the estimated integral." << endl;
//
//  Create the grid.
//
  const int dimension = 2;
  const int outputs = 1;
  const int precision = 10;
  TasGrid::TasmanianSparseGrid grid;

  grid.makeGlobalGrid ( dimension, outputs, precision, TasGrid::type_basis,
    TasGrid::rule_clenshawcurtis );
//
//  Report the number of points in the grid.
//  Also report the number of points at which function values are needed.
//  For this example, these values should be the same.
//
  int N = grid.getNumNeededPoints();

  cout << endl;
  cout << "  The grid has a total of " << grid.getNumPoints() 
       << " points." << endl;
  cout << "  For interpolation, function values are still needed for " 
       << N << " grid points." << endl;
//
//  Get the points at which the function must be evaluated.
//              
  double *points = 0;

  grid.getNeededPoints( points );
//
//  Evaluate the function to supply the needed information.
//
  double *values = new double[N];
  double x;
  double y;

  for ( int i = 0; i < N; i++ )
  {
    x = points[2*i];
    y = points[2*i+1];
    values[i] = exp ( - x * x ) * cos ( y );
  }

  grid.loadNeededPoints ( values );
//
//  Compare f(x,y) to its interpolant at some random points in [-1,+1]^2.
//
  double estimate;
  double exact;
  double xy[2] = { 0.3, 0.7 };

  cout << endl;
  cout << "  Compare F(X,Y) and its interpolant." << endl;

  for ( int test = 0; test < 5; test++ )
  {
    x = 2.0 * drand48 ( ) - 1.0;
    y = 2.0 * drand48 ( ) - 1.0;
    xy[0] = x;
    xy[1] = y;
    grid.evaluate ( xy, &estimate );
    exact = exp ( - x * x ) * cos ( y );
    cout.precision(17);
    cout << endl;
    cout << "  (X,Y) = ( " << x << ", " << y << " )." << endl;
    cout << "  F(X,Y) =              " << exact << endl;
    cout << "  Interpolant(F)(X,Y) = " << estimate << endl;
    cout.precision(4);
    cout << "  Error =               " << std::scientific 
         << fabs ( estimate - exact ) << endl;
  }
//
//  Compare the integral of f(x,y) to the estimated integral.
//
  exact = 2.513723354063905e+00;

  grid.integrate( &estimate );
  cout.precision(17);
  cout << endl;
  cout << "  The estimated integral is: " 
       << std::scientific << estimate << endl;
  cout << "  The exact integral is:     " << std::scientific << exact << endl;
  cout.precision(4);
  cout << "  The error is: " << std::scientific 
       << fabs ( estimate - exact ) << endl;
//
//  Free memory.
//
  delete [] points;
  delete [] values;

  return;
}
//****************************************************************************80

void tsg_test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST03 creates a Gauss-Legendre sparse grid in 2D.
//
//  Discussion:
//
//    Create a sparse grid quadrature rule from the Gauss-Legendre rule,
//    with the property that it is exact for polynomials up to degree 10.
//
//    Estimate the integral of f(x,y) = exp ( -x^2 ) * cos ( y ).
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{
  cout << endl;
  cout << "TSG_TEST03" << endl;
  cout << "  Make a 2D Smolyak grid using the Gauss-Legendre rule" << endl;
  cout << "  over the region [-1,+1]^2." << endl;
  cout << "  This grid interpolates polynomials of total degree 10 or less." 
       << endl;
  cout << "  Consider f(x,y) = exp ( -x^2 ) * cos ( y )" << endl;
  cout << "  Compare the integral of f(x,y) to the estimated integral." << endl;
//
//  Create the grid.
//
  const int dimension = 2;
  const int outputs = 0;
  const int precision = 10;
  TasGrid::TasmanianSparseGrid grid;

  grid.makeGlobalGrid ( dimension, outputs, precision, TasGrid::type_basis,
    TasGrid::rule_gausslegendre );

  cout << "  " << endl;
  cout << "  The grid uses " << grid.getNumPoints() << " points." << endl;
//
//  Get the points and weights of the grid.
//
  double *points = 0;
  double *weights = 0;

  grid.getPoints ( points );
  grid.getWeights ( weights );
//
//  Estimate the integral.
//
  double sum = 0.0;
  double x;
  double y;

  for ( int i = 0; i < grid.getNumPoints(); i++ )
  {
    x = points[2*i];
    y = points[2*i+1];
    sum += weights[i] * exp ( - x * x ) * cos ( y );
  }
//
//  Report.
//
  double exact = 2.513723354063905e+00;

  cout.precision(17);
  cout << "  The estimated integral is: " << std::scientific << sum << endl;
  cout << "  The exact integral is:     " << std::scientific << exact << endl;
  cout.precision(4);
  cout << "  The error is:              " << std::scientific 
       << fabs ( sum - exact ) << endl;
//
//  Free memory.
//
  delete [] points;
  delete [] weights;

  return;
}
//****************************************************************************80

void tsg_test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST04 makes an anisotropic Gauss-Gegenbauer sparse grid.
//
//  Discussion:
//
//    Make a sparse grid quadrature rule that uses the Gauss-Gegenbauer points.
//
//    The grid should be anisotropic, using 8 times as many points in the
//    y direction as in the x direction.
//
//    Use the rule to estimate the integral of the function 
//      f(x,y) = ( x - 2 )^3 * exp ( -y^2 )
//    with the Gegenbauer-type product weight
//      rho(x,y) = ( 1 - x^2 )^0.4 * ( 1 - y^2 )^0.4.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    18 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{
  cout << endl;
  cout << "TSG_TEST04" << endl;
  cout << "  Make a 2D Smolyak grid using the Gauss-Gegenbauer rule" << endl;
  cout << "  over the region [-1,+1]^2." << endl;
  cout << "  This grid should be 8 times as dense in the Y direction as in X." 
       << endl;
  cout << "  Consider f(x,y) = exp ( x - 2 )^3 * exp ( -y^2 )" << endl;
  cout << "  Use the Gegenbauer-type product weight" << endl;
  cout << "    rho(x,y) = ( 1 - x^2 )^0.4 * ( 1 - y^2 )^0.4" << endl;
  cout << "  Compare the integral of f(x,y) * rho(x,y) to the estimate." 
       << endl;
//
//  Create the grid.
//
  const int dimension = 2;
  const int outputs = 0;
  const int depth = 16;
  const int anisotropic_weights[2] = { 8, 1 };
  const double alpha = 0.4;
  TasGrid::TasmanianSparseGrid grid;

  grid.makeGlobalGrid ( dimension, outputs, depth, TasGrid::type_level,
    TasGrid::rule_gaussgegenbauer, anisotropic_weights, &alpha );

  cout << "  " << endl;
  cout << "  The grid uses " << grid.getNumPoints() << " points." << endl;
//
//  Get the points and weights of the grid.
//
  double *points = 0;
  double *weights = 0;

  grid.getPoints ( points );
  grid.getWeights ( weights );
//
//  Estimate the integral.
//
  double sum = 0.0;
  double x;
  double y;

  for ( int i = 0; i < grid.getNumPoints(); i++ )
  {
    x = points[2*i];
    y = points[2*i+1];
    sum += weights[i] * ( x - 2.0 ) * ( x - 2.0 ) * ( x - 2.0 ) 
      * exp ( - y * y );
  }
//
//  Report.
//
  double exact = -2.029979511486524e+01;

  cout.precision(17);
  cout << "  The estimated integral is: " << std::scientific << sum << endl;
  cout << "  The exact integral is:     " << std::scientific << exact << endl;
  cout.precision(4);
  cout << "  The error is:              " << std::scientific 
       << fabs ( sum - exact ) << endl;
//
//  Free memory.
//
  delete [] points;
  delete [] weights;

  return;
}
//****************************************************************************80

void tsg_test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST05 interpolates a function using an adaptive piecewise polynomial.
//
//  Discussion:
//
//    Interpolate the function f(x,y) = exp ( -x^2 ) * cos ( y ) 
//    using adaptive piecewise local quadratic polynomials over [0,1]^2.
//
//    This requires transforming the grid from [-1,+1]^2 to [0,1]^2,
//
//    Evaluate the interpolant at (0.3,07).
//
//    The exact value of f(x,y) there is 6.990131267703512e-01.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    Miro Stoyanov
//
{
  cout << endl;
  cout << "TSG_TEST05" << endl;
  cout << "  Interpolate f(x,y) = exp ( -x^2 ) * cos ( y ) " << endl;
  cout << "  using adaptive piecewise local quadratic polynomials over [0,1]^2."
       << endl;
  cout << "  (This requires transform the grid from [-1,+1]^2 to [0,1]^2.)" 
       << endl;
  cout << "  Evaluate the interpolant at (0.3,07)." << endl;
  cout << "  The exact value of f(x,y) there is 6.990131267703512e-01." << endl;
//
//  Create the grid.
//
  const int dimension = 2;
  const int outputs = 1;
  const int initial_level = 4;
  const int order = 2;
  TasGrid::TasmanianSparseGrid grid;

  grid.makeLocalPolynomialGrid ( dimension, outputs, initial_level, order,
    TasGrid::rule_pwpolynomial );
//
//  The grid is initially defined on [-1,+1]^2.
//  Transform it to [0,1]^2.
//
  const double a[2] = { 0.0, 0.0 };
  const double b[2] = { 1.0, 1.0 };

  grid.setTransformAB( a, b );
//
//  Carry out the adaptive interpolation repeatedly, until the tolerance
//  is satisfied.
//
  double exact = 6.990131267703512e-01;
  const double tolerance = 1.0E-6;              
  int iteration = 0;
  double x;
  double xy[2] = { 0.3, 0.7 };
  double y;

  cout << endl;
  cout << "  Iteration  Samples  f(x,y)                   I(f)(x,y)";
  cout << "                error" << endl;
  cout << endl;

  while ( grid.getNumNeededPoints() > 0 )
  {
    int N = grid.getNumNeededPoints();                    
    double *points = 0;
    double *values = new double[N];
//
//  Request the locations of points where function values are needed.
//  Evaluate the function at those points.
//  Return the values.
//
    grid.getNeededPoints ( points );
                        
    for ( int i = 0; i < N; i++ )
    {
      x = points[2*i];
      y = points[2*i+1];
      values[i] = exp ( - x * x ) * cos ( y );
    }

    grid.loadNeededPoints ( values );
//
//  Evaluate the interpolant at an arbitrary point XY.
//
    double I;

    grid.evaluate ( xy, &I );

    cout.precision(17);
    cout << "  " << setw(9) << iteration 
         << "  " << setw(7) << grid.getNumPoints()
         << "  " << setw(10) << exact 
         << "  " << setw(9) << I;
    cout.precision(4);
    cout << "  " << setw(10) << fabs ( exact - I ) << endl;
//
//  If the tolerance was not reached, allow the grid to take
//  another refinement step.
//
    grid.setRefinement ( tolerance, TasGrid::refine_classic );

    iteration++;
//
//  Free memory.
//
    delete [] points;
    delete [] values;
  }

  return;
}
//****************************************************************************80

void tsg_test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TSG_TEST06 compares 2D CC sparse grids using different growth rules.
//
//  Licensing:
//
//    This code is distributed under the GNU GPL license.
//
//  Modified:
//
//    16 November 2013
//
//  Author:
//
//    John Burkardt
//
{
  int dim = 2;
  TasGrid::TasmanianSparseGrid grid;
  int level;
  int n1;
  int n2;
  int n3;
  int outputs = 0;

  cout << "\n";
  cout << "TSG_TEST06\n";
  cout << "  Generate Clenshaw-Curtis sparse grids in 2D.\n";
  cout << "  For the growth argument, compare all three options.\n";
  cout << "\n";
  cout << "  HYPER: Hyperbolic cross\n";
  cout << "  LEVEL: Smolyak level\n";
  cout << "  BASIS: Polynomial exactness criterion.\n";
  cout << "\n";
  cout << "  Level     HYPER     LEVEL     BASIS\n";
  cout << "\n";

  for ( level = 0; level <= 11; level++ )
  {
    if ( level == 0 )
    {
      grid.makeGlobalGrid ( dim, outputs, level, TasGrid::type_basis,
        TasGrid::rule_clenshawcurtis );
      n3 = grid.getNumPoints ( );

      cout << "  " << setw(5) << level
           << "  " << "       *"
           << "  " << "       *"
           << "  " << setw(8) << n3 << "\n";
    }
    else
    {
      grid.makeGlobalGrid ( dim, outputs, level, TasGrid::type_hyperbolic,
        TasGrid::rule_clenshawcurtis );
      n1 = grid.getNumPoints ( );

      grid.makeGlobalGrid ( dim, outputs, level, TasGrid::type_level,
        TasGrid::rule_clenshawcurtis );
      n2 = grid.getNumPoints ( );

      grid.makeGlobalGrid ( dim, outputs, level, TasGrid::type_basis,
        TasGrid::rule_clenshawcurtis );
      n3 = grid.getNumPoints ( );

      cout << "  " << setw(5) << level
           << "  " << setw(8) << n1
           << "  " << setw(8) << n2
           << "  " << setw(8) << n3 << "\n";
    }
  }

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

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", 
    tm_ptr );

  std::cout << time_buffer << endl;

  return;
# undef TIME_SIZE
}
