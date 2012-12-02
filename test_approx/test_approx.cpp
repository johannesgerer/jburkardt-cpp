# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cstring>
# include <ctime>
# include <fstream>

using namespace std;

# include "test_approx.hpp"

//****************************************************************************80

void p00_dat ( int prob, int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P00_DAT returns the data vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Input, int DATA_NUM, the number of data points,
//    as specified by P00_DATA_NUM.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  if ( prob == 1 )
  {
    p01_dat ( data_num, x, y );
  }
  else if ( prob == 2 )
  {
    p02_dat ( data_num, x, y );
  }
  else if ( prob == 3 )
  {
    p03_dat ( data_num, x, y );
  }
  else if ( prob == 4 )
  {
    p04_dat ( data_num, x, y );
  }
  else if ( prob == 5 )
  {
    p05_dat ( data_num, x, y );
  }
  else if ( prob == 6 )
  {
    p06_dat ( data_num, x, y );
  }
  else if ( prob == 7 )
  {
    p07_dat ( data_num, x, y );
  }
  else if ( prob == 8 )
  {
    p08_dat ( data_num, x, y );
  }
  else if ( prob == 9 )
  {
    p09_dat ( data_num, x, y );
  }
  else if ( prob == 10 )
  {
    p10_dat ( data_num, x, y );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_DAT - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }
  return;
}
//****************************************************************************80

int p00_data_num ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_DATA_NUM returns the dimension of the data vector for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, int P00_DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  if ( prob == 1 )
  {
    data_num = p01_data_num ( );
  }
  else if ( prob == 2 )
  {
    data_num = p02_data_num ( );
  }
  else if ( prob == 3 )
  {
    data_num = p03_data_num ( );
  }
  else if ( prob == 4 )
  {
    data_num = p04_data_num ( );
  }
  else if ( prob == 5 )
  {
    data_num = p05_data_num ( );
  }
  else if ( prob == 6 )
  {
    data_num = p06_data_num ( );
  }
  else if ( prob == 7 )
  {
    data_num = p07_data_num ( );
  }
  else if ( prob == 8 )
  {
    data_num = p08_data_num ( );
  }
  else if ( prob == 9 )
  {
    data_num = p09_data_num ( );
  }
  else if ( prob == 10 )
  {
    data_num = p10_data_num ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_DATA_NUM - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }

  return data_num;
}
//****************************************************************************80

int p00_prob_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P00_PROB_NUM returns the number of test problems.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int PROB_NUM, the number of test problems.
//
{
  int prob_num;

  prob_num = 10;

  return prob_num;
}
//****************************************************************************80

void p00_story ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_STORY prints the "story" for any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
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
  if ( prob == 1 )
  {
    p01_story ( );
  }
  else if ( prob == 2 )
  {
    p02_story ( );
  }
  else if ( prob == 3 )
  {
    p03_story ( );
  }
  else if ( prob == 4 )
  {
    p04_story ( );
  }
  else if ( prob == 5 )
  {
    p05_story ( );
  }
  else if ( prob == 6 )
  {
    p06_story ( );
  }
  else if ( prob == 7 )
  {
    p07_story ( );
  }
  else if ( prob == 8 )
  {
    p08_story ( );
  }
  else if ( prob == 9 )
  {
    p09_story ( );
  }
  else if ( prob == 10 )
  {
    p10_story ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_STORY - Fatal error!\n";
    cerr << "  Unexpected input value of PROB.\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

string p00_title ( int prob )

//****************************************************************************80
//
//  Purpose:
//
//    P00_TITLE returns the title of any problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int PROB, the number of the desired test problem.
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  if ( prob == 1 )
  {
    title = p01_title ( );
  }
  else if ( prob == 2 )
  {
    title = p02_title ( );
  }
  else if ( prob == 3 )
  {
    title = p03_title ( );
  }
  else if ( prob == 4 )
  {
    title = p04_title ( );
  }
  else if ( prob == 5 )
  {
    title = p05_title ( );
  }
  else if ( prob == 6 )
  {
    title = p06_title ( );
  }
  else if ( prob == 7 )
  {
    title = p07_title ( );
  }
  else if ( prob == 8 )
  {
    title = p08_title ( );
  }
  else if ( prob == 9 )
  {
    title = p09_title ( );
  }
  else if ( prob == 10 )
  {
    title = p10_title ( );
  }
  else
  {
    cerr << "\n";
    cerr << "P00_TITLE - Fatal error!\n";
    cerr << "  Illegal problem number = " << prob << "\n";
    exit ( 1 );
  }

  return title;
}
//****************************************************************************80

void p01_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P01_DAT returns the data vector for problem 1.
//
//  Discussion:
//
//    The X data is measured in days, and the Y data represents the
//    observed position of Mars in a heliocentric coordinate system.
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Samuel Conte, Carl deBoor,
//    Elementary Numerical Analysis,
//    McGraw Hill, 1972, page 217.
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[10] = {
    1250.5, 1260.5, 1270.5, 1280.5, 1290.5, 
    1300.5, 1310.5, 1320.5, 1330.5, 1340.5 };
  double y_save[10] = {
    1.39140, 1.37696, 1.34783, 1.30456, 1.24787, 
    1.17862, 1.09776, 1.00636, 0.90553, 0.79642 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p01_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_DATA_NUM returns the dimension of the data vector for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 10;

  return data_num;
}
//****************************************************************************80

void p01_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_STORY prints the "story" for problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This example is due to deBoor.\n";
  cout << "  For this example, X is measured in days,\n";
  cout << "  and Y records the observed position of Mars in a heliocentric\n";
  cout << "  coordinate system.\n";

  return;
}
//****************************************************************************80

string p01_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P01_TITLE returns the title of problem 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "DeBoor example, Mars position";

  return title;
}
//****************************************************************************80

void p02_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P02_DAT returns the data vector for problem 2.
//
//  Discussion:
//
//    The data lies roughly along a straight line.  Polynomial
//    interpolation is inappropriate.  Instead, a least squares
//    approximation should be sought, of the form:
//
//      F(X) = A + B * X
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Samuel Conte, Carl deBoor,
//    Elementary Numerical Analysis,
//    McGraw Hill, 1972, page 217.
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[11] = {
    1.0, 2.0, 3.0, 4.0,  5.0, 
    6.0, 7.0, 8.0, 9.0, 10.0, 
   11.0 };
  double y_save[11] = {
    0.00, 0.60, 1.77, 1.92, 3.31, 
    3.52, 4.59, 5.31, 5.79, 7.06, 
    7.17 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p02_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_DATA_NUM returns the dimension of the data vector for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 11;

  return data_num;
}
//****************************************************************************80

void p02_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_STORY prints the "story" for problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This example is due to deBoor.\n";
  cout << "  The data lies roughly along a straight line.  Polynomial\n";
  cout << "  interpolation is inappropriate.  Instead, a least squares\n";
  cout << "  approximation should be sought, of the form:\n";
  cout << "    F(X) = A + B * X\n";

  return;
}
//****************************************************************************80

string p02_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P02_TITLE returns the title of problem 2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "DeBoor example, roughly linear data";

  return title;
}
//****************************************************************************80

void p03_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P03_DAT returns the data vector for problem 3.
//
//  Discussion:
//
//    The data is all zero except for a single value of 1 in the center.
//    This data set is interesting because an interpolation method that
//    is "local" will produce an interpolating curve that is exactly
//    zero over most of the outlying intervals, whereas a nonlocal
//    interpolation method may produce a curve that "wiggles" over the
//    entire interpolation interval.
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[11] = {
     0.0, 1.0, 2.0, 3.0, 4.0, 
     5.0, 6.0, 7.0, 8.0, 9.0, 
    10.0 };
  double y_save[11] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 
    1.0, 0.0, 0.0, 0.0, 0.0, 
    0.0 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p03_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_DATA_NUM returns the dimension of the data vector for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 11;

  return data_num;
}
//****************************************************************************80

void p03_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_STORY prints the "story" for problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  The data is all zero except for a single value of 1 in the center.\n";
  cout << "  This data set is interesting because an interpolation method that\n";
  cout << "  is \"local\" will produce an interpolating curve that is exactly\n";
  cout << "  zero over most of the outlying intervals, whereas a nonlocal\n";
  cout << "  interpolation method may produce a curve that \"wiggles\" over the\n";
  cout << "  entire interpolation interval.\n";

  return;
}
//****************************************************************************80

string p03_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P03_TITLE returns the title of problem 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "The pulse data, 0 0 0 0 0 1 0 0 0 0 0";

  return title;
}
//****************************************************************************80

void p04_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P04_DAT returns the data vector for problem 4.
//
//  Discussion:
//
//    Theoretically, the data is a step, 0 to the left of 5, and 1
//    to the right.  To keep things simple, the data is defined
//    to be 0 up to 5 - RADIUS, 1/2 at 5, 1 at 5 + RADIUS and beyond,
//    with RADIUS set to a "small" value, currently 0.01.
//    Some interpolation methods will violently overreact to this
//    jump.
//
//    The X data is NOT equally spaced because of the handling of the pulse.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[13] = {
    0.0, 1.0, 2.0, 3.0, 4.0, 
    4.99, 5.00, 5.01, 6.0, 7.0, 
    8.0, 9.0, 10.0 };
  double y_save[13] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 
    0.0, 0.5, 1.0, 1.0, 1.0, 
    1.0, 1.0, 1.0 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p04_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_DATA_NUM returns the dimension of the data vector for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 13;

  return data_num;
}
//****************************************************************************80

void p04_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_STORY prints the "story" for problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  Theoretically, the data is a step, 0 to the left of 5, and 1\n";
  cout << "  to the right.  To keep things simple, the data is defined\n";
  cout << "  to be 0 up to 5 - RADIUS, 1/2 at 5, 1 at 5 + RADIUS and beyond,\n";
  cout << "  with RADIUS set to a \"small\" value, currently 0.01.\n";
  cout << "  Some interpolation methods will violently overreact to this jump.\n";

  return;
}
//****************************************************************************80

string p04_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P04_TITLE returns the title of problem 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "The jump data, 0 0 0 0 0 1/2 1 1 1 1 1";

  return title;
}
//****************************************************************************80

void p05_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P05_DAT returns the data vector for problem 5.
//
//  Discussion:
//
//    This example is due to deBoor.
//    This data represents a property of titanium as a function of temperature.
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Carl deBoor,
//    A Practical Guide to Splines,
//    Springer Verlag.
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[49] = {
    595.0, 
    605.0,  615.0,  625.0,  635.0,  645.0, 
    655.0,  665.0,  675.0,  685.0,  695.0, 
    705.0,  715.0,  725.0,  735.0,  745.0, 
    755.0,  765.0,  775.0,  785.0,  795.0, 
    805.0,  815.0,  825.0,  835.0,  845.0, 
    855.0,  865.0,  875.0,  885.0,  895.0, 
    905.0,  915.0,  925.0,  935.0,  945.0, 
    955.0,  965.0,  975.0,  985.0,  995.0, 
   1005.0, 1015.0, 1025.0, 1035.0, 1045.0, 
   1055.0, 1065.0, 1075.0 };
  double y_save[49] = {
    0.644, 
    0.622, 0.638, 0.649, 0.652, 0.639, 
    0.646, 0.657, 0.652, 0.655, 0.644, 
    0.663, 0.663, 0.668, 0.676, 0.676, 
    0.686, 0.679, 0.678, 0.683, 0.694, 
    0.699, 0.710, 0.730, 0.763, 0.812, 
    0.907, 1.044, 1.336, 1.881, 2.169, 
    2.075, 1.598, 1.211, 0.916, 0.746, 
    0.672, 0.627, 0.615, 0.607, 0.606, 
    0.609, 0.603, 0.601, 0.603, 0.601, 
    0.611, 0.601, 0.608 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p05_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_DATA_NUM returns the dimension of the data vector for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 49;

  return data_num;
}
//****************************************************************************80

void p05_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_STORY prints the "story" for problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This example is due to deBoor.\n";
  cout << "  This data represents a temperature dependent property of titanium.\n";

  return;
}
//****************************************************************************80

string p05_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P05_TITLE returns the title of problem 5.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "DeBoor's Titanium property data";

  return title;
}
//****************************************************************************80

void p06_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P06_DAT returns the data vector for problem 6.
//
//  Discussion:
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  int i;
  int j;
  int n;
  int num_int = 5;

  n = 0;
  x[n] = 0.0;
  y[n] = 0.0;

  for ( i = 1; i <= num_int; i++ )
  {
    for ( j = 1; j <= i; j++ )
    {
      n = n + 1;
      x[n] = ( double ) ( i - 1 ) + 0.5 * ( double ) ( j ) / ( double ) ( i );
      y[n] = ( double ) ( j ) / ( double ) ( i );
    }

    for ( j = 1; j <= i; j++ )
    {
      n = n + 1;
      x[n] = ( double ) ( i - 1 ) + 0.5
        + 0.5 * ( double ) ( j ) / ( double ) ( i );
      y[n] = 1.0 - ( double ) ( j ) / ( double ) ( i );
    }

  }

  return;
}
//****************************************************************************80

int p06_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_DATA_NUM returns the dimension of the data vector for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;
  int num_int = 5;

  data_num = 1 + num_int * ( num_int + 1 );

  return data_num;
}
//****************************************************************************80

void p06_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_STORY prints the "story" for problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This is a data vector.\n";

  return;
}
//****************************************************************************80

string p06_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P06_TITLE returns the title of problem 6.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "The Sawtooth data";

  return title;
}
//****************************************************************************80

void p07_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P07_DAT returns the data vector for problem 7.
//
//  Discussion:
//
//    The X data is NOT equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[9] = {
    0.0,  0.1,  0.2,  0.3, 0.4, 
    0.5,  0.6,  0.8,  1.0 };
  double y_save[9] = {
     0.0,  0.9,  0.95, 0.9, 0.1, 
     0.05, 0.05, 0.2,  1.0 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p07_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_DATA_NUM returns the dimension of the data vector for problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 9;

  return data_num;
}
//****************************************************************************80

void p07_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_STORY prints the "story" for problem 7
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This is a data vector.\n";

  return;
}
//****************************************************************************80

string p07_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P07_TITLE returns the title of problem 7.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "Concavity test data";

  return title;
}
//****************************************************************************80

void p08_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P08_DAT returns the data vector for problem 8.
//
//  Discussion:
//
//    This example is due to Pierre Blais.
//
//    Data is only available over the interval [0, 238], but extrapolation
//    is to be used to extend the approximate function to a maximum argument
//    of 1023.  The behavior of the extrapolated curve is of great interest.
//
//    The X data is NOT equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[12] = {
      0.0,  71.0,  104.0,  135.0, 145.0, 
    160.0, 181.0,  193.0,  205.0, 215.0, 
    225.0, 238.0 };
  double y_save[12] = {
      0.0000,   7.7554,  19.7062,  35.53786,  42.91537, 
     54.7752,  66.75865, 78.49286, 89.76833, 101.746, 
    113.4824, 135.4566 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p08_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_DATA_NUM returns the dimension of the data vector for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 12;

  return data_num;
}
//****************************************************************************80

void p08_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_STORY prints the "story" for problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This example is due to Pierre Blais.\n";
  cout << "  Data is only available over the interval [0, 238], but extrapolation\n";
  cout << "  is to be used to extend the approximate function to a maximum argument\n";
  cout << "  of 1023.  The behavior of the extrapolated curve is of great interest.\n";

  return;
}
//****************************************************************************80

string p08_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P08_TITLE returns the title of problem 8.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "Extrapolation test data";

  return title;
}
//****************************************************************************80

void p09_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P09_DAT returns the data vector for problem 9.
//
//  Discussion:
//
//    This example is due to Max Waldmeier.
//
//    This data represents a measure of sunspot activity over the
//    years 1700 to 1960.  The X value is the year, and the Y value
//    is a measure of the sunspot activity, which is usually, but
//    not always, an integer.
//
//    The X data is equally spaced.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Max Waldmeier,
//    The Sunspot-Activity in the Years 1610-1960,
//    Shulthess, Zurich, 1961.
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  int i;
  double y_save[261] = {
    5.0,  11.0,  16.0,  23.0,  36.0, 
   58.0,  29.0,  20.0,  10.0,   8.0, 
    3.0,   0.0,   0.0,   2.0,  11.0, 
   27.0,  47.0,  63.0,  60.0,  39.0, 
   28.0,  26.0,  22.0,  11.0,  21.0, 
   40.0,  78.0, 122.0, 103.0,  73.0, 
   47.0,  35.0,  11.0,   5.0,  16.0, 
   34.0,  70.0,  81.0, 111.0, 101.0, 
   73.0,  40.0,  20.0,  16.0,   5.0, 
   11.0,  22.0,  40.0,  60.0,  80.9, 
   83.4,  47.7,  47.8,  30.7,  12.2, 
    9.6,  10.2,  32.4,  47.6,  54.0, 
   62.9,  85.9,  61.2,  45.1,  36.4, 
   20.9,  11.4,  37.8,  69.8, 106.1, 
  100.8,  81.6,  66.5,  34.8,  30.6, 
    7.0,  19.8,  92.5, 154.4, 125.9, 
   84.8,  68.1,  38.5,  22.8,  10.2,
   24.1,  82.9, 132.0, 130.9, 118.1, 
   89.9,  66.6,  60.0,  46.9,  41.0, 
   21.3,  16.0,   6.4,   4.1,   6.8, 
   14.5,  34.0,  45.0,  43.1,  47.5, 
   42.2,  28.1,  10.1,   8.1,   2.5, 
    0.0,   1.4,   5.0,  12.2,  13.9, 
   35.4,  45.8,  41.1,  30.1,  23.9, 
   15.6,   6.6,   4.0,   1.8,   8.5, 
   16.6,  36.3,  49.6,  64.2,  67.0, 
   70.9,  47.8,  27.5,   8.5,  13.2, 
   56.9, 121.5, 138.3, 103.2,  85.7, 
   64.6,  36.7,  24.2,  10.7,  15.0, 
   40.1,  61.5,  98.5, 124.7,  96.3, 
   66.6,  64.5,  54.1,  39.0,  20.6, 
    6.7,   4.3,  22.7,  54.8,  93.8, 
   95.8,  77.2,  59.1,  44.0,  47.0, 
   30.5,  16.3,   7.3,  37.6,  74.0, 
  139.0, 111.2, 101.6,  66.2,  44.7, 
   17.0,  11.3,  12.4,   3.4,   6.0, 
   32.3,  54.3,  59.7,  63.7,  63.5, 
   52.2,  25.4,  13.1,   6.8,   6.3, 
    7.1,  35.6,  73.0,  85.1,  78.0, 
   64.0,  41.8,  26.2,  26.7,  12.1, 
    9.5,   2.7,   5.0,  24.4,  42.0, 
   63.5,  53.8,  62.0,  48.5,  43.9, 
   18.6,   5.7,   3.6,   1.4,   9.6, 
   47.4,  57.1, 103.9,  80.6,  63.6, 
   37.6,  26.1,  14.2,   5.8,  16.7, 
   44.3,  63.9,  69.0,  77.8,  64.9, 
   35.7,  21.2,  11.1,   5.7,   8.7, 
   36.1,  79.7, 114.4, 109.6,  88.8, 
   67.8,  47.5,  30.6,  16.3,   9.6, 
   33.2,  92.6, 151.6, 136.3, 134.7, 
   83.9,  69.4,  31.5,  13.9,   4.4, 
   38.0, 141.7, 190.2, 184.8, 159.0, 
  112.3 };

  for ( i = 0; i < 261; i++ )
  {
    x[i] = ( double ) ( 1700 + i );
  }

  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p09_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_DATA_NUM returns the dimension of the data vector for problem 9.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 261;

  return data_num;
}
//****************************************************************************80

void p09_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_STORY prints the "story" for problem 09
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
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
  cout << "\n";
  cout << "  This example is due to Max Waldmeier.\n";
  cout << "  This data represents a measure of sunspot activity over the\n";
  cout << "  years 1700 to 1960.  The X value is the year, and the Y value\n";
  cout << "  is a measure of the sunspot activity, which is usually, but\n";
  cout << "  not always, an integer.\n";

  return;
}
//****************************************************************************80

string p09_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P09_TITLE returns the title of problem 09.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 August 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "Sunspot data, 1700-1960.";

  return title;
}
//****************************************************************************80

void p10_dat ( int data_num, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    P10_DAT returns the data vector for problem 10.
//
//  Discussion:
//
//    100 uniformly random X values between -2 and 5 were selected,
//    and the formula Y = 2 + 5 * X + 10 * N(0,1) was evaluated, where
//    N(0,1) represents random normal values with 0 mean and unit variance.
//
//    The original data was unsorted, but this caused problems for various
//    approximation codes, so the data has now been sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DATA_NUM, the dimension of the data vector.
//
//    Output, double X[DATA_NUM], the abscissa data.
//
//    Output, double Y[DATA_NUM], the ordinate data.
//
{
  double x_save[100] = {
    -1.9935, -1.8884, -1.8813, -1.8786, -1.8714, 
    -1.8108, -1.6614, -1.5991, -1.4423, -1.3441, 
    -1.1730, -1.1575, -0.9886, -0.8196, -0.7771, 
    -0.7706, -0.7270, -0.7243, -0.7228, -0.6235, 
    -0.6023, -0.1650, -0.1117, -0.0570,  0.0923, 
     0.1132,  0.3157,  0.3343,  0.3356,  0.3827, 
     0.4019,  0.4570,  0.4866,  0.6041,  0.6803, 
     0.8395,  0.8734,  0.9115,  0.9190,  1.0716, 
     1.1649,  1.2261,  1.3118,  1.3484,  1.3698, 
     1.5721,  1.6258,  1.7749,  1.8234,  1.9732, 
     1.9784,  2.0022,  2.0231,  2.0398,  2.0401, 
     2.0757,  2.0824,  2.2051,  2.2634,  2.3076, 
     2.3660,  2.3931,  2.4175,  2.5928,  2.6235, 
     2.8109,  2.8898,  2.9561,  2.9669,  3.1161, 
     3.1430,  3.1749,  3.2373,  3.4166,  3.6131, 
     3.6635,  3.6669,  3.7486,  3.7530,  3.7862, 
     3.8021,  3.8298,  3.8696,  3.9051,  3.9509, 
     4.0678,  4.1038,  4.1557,  4.2456,  4.2564, 
     4.3228,  4.3403,  4.4157,  4.4998,  4.5767, 
     4.6312,  4.7029,  4.7717,  4.8424,  4.8915  };
  double y_save[261] = {
     -1.8918,  -8.6202,  -0.4152,  -4.6967,  -2.4144, 
    -21.8853, -16.5097, -10.4657,  -4.1150,   6.5665, 
     -6.7649,   8.8276,   1.8109,   9.6428,  -0.6165, 
     -8.4213, -16.4491,  -0.0667,   6.5713,  -4.0436, 
     -6.4194,  -1.9114,  -9.5246,  -3.2154,   0.6541, 
      3.0247,   2.9410,   9.7848,   4.7713,  22.0539, 
      7.1301,  22.3302,  -2.7979,  10.2864,   2.7993, 
     12.1989,  12.3065, -15.3026,  -6.6749,  -7.0521, 
     11.8429,  22.8325,   5.2912,  16.8654,  14.3045, 
     -0.6553,  14.1038,   3.3556,  26.2797,  11.5406, 
     28.2522,   7.7605,  18.0102,  11.5710,  -8.0188, 
      2.5576,  18.5371,  12.4771,   2.1298,   7.2745, 
     16.3256,   4.0354,  23.8374,   8.5572,  33.2063, 
      5.2561,  18.4409,   1.5702,   9.5982,  11.6483, 
     21.7285,  27.2958,  21.1916,  15.3527,  28.2207, 
     28.3065,  21.5368,  26.4556,  24.8933,  11.0617, 
     28.6064,  14.5770,  15.3088,  23.2952,  18.6796, 
     21.0208,  28.4730,  33.2469,  21.2484,  26.5592, 
     21.2314,  25.9979,  28.4785,  18.3307,  27.6318, 
     31.1673,  26.4379,  43.1573,  20.1264,  19.0873 };

  r8vec_copy ( data_num, x_save, x );
  r8vec_copy ( data_num, y_save, y );

  return;
}
//****************************************************************************80

int p10_data_num ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_DATA_NUM returns the dimension of the data vector for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, int DATA_NUM, the dimension of the data vector.
//
{
  int data_num;

  data_num = 100;

  return data_num;
}
//****************************************************************************80

void p10_story ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_STORY prints the "story" for problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
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
  cout << "\n";
  cout << "  100 uniformly random X values between -2 and 5 were selected, and\n";
  cout << "  the formula Y = 2 + 5 * X + 10 * N(0,1) was evaluated, where N(0,1)\n";
  cout << "  represents random normal values with 0 mean and unit variance.\n";

  return;
}
//****************************************************************************80

string p10_title ( )

//****************************************************************************80
//
//  Purpose:
//
//    P10_TITLE returns the title of problem 10.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 February 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, string TITLE, the title of the problem.
//
{
  string title;

  title = "Y = 2 + 5*X + 10*N(0,1).";

  return title;
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
//    Output, double A2[N], the copy of A1.
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

void r8vec2_print ( int n, double a1[], double a2[], string title )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_PRINT prints an R8VEC2.
//
//  Discussion:
//
//    An R8VEC2 is a dataset consisting of N pairs of real values, stored
//    as two separate vectors A1 and A2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 November 2002
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of components of the vector.
//
//    Input, double A1[N], double A2[N], the vectors to be printed.
//
//    Input, string TITLE, a title.
//
{
  int i;

  cout << "\n";
  cout << title << "\n";
  cout << "\n";
  for ( i = 0; i <= n - 1; i++ )
  {
    cout << setw(6)  << i
         << ": " << setw(14) << a1[i]
         << "  " << setw(14) << a2[i] << "\n";
  }

  return;
}
//****************************************************************************80

void r8vec2_write ( string output_filename, int n, double x[], double y[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC2_WRITE writes an R8VEC2 file.
//
//  Discussion:
//
//    An R8VEC2 is a pair of vectors of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string OUTPUT_FILENAME, the output filename.
//
//    Input, int N, the number of points.
//
//    Input, double X[N], Y[N], the data.
//
{
  int j;
  ofstream output;
//
//  Open the file.
//
  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "R8VEC2_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    exit ( 1 );
  }
//
//  Write the data.
//
  for ( j = 0; j < n; j++ )
  {
    output << "  " << setw(24) << setprecision(16) << x[j] 
           << "  " << setw(24) << setprecision(16) << y[j] << "\n";
  }
//
//  Close the file.
//
  output.close ( );

  return;
}
