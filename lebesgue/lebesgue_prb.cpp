# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "lebesgue.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );
void test05 ( );
void test06 ( );
void test07 ( );
void test08 ( );
void test09 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LEBESGUE_PRB.
//
//  Discussion:
//
//    LEBESGUE_PRB tests the LEBESGUE library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "LEBESGUE_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the LEBESGUE library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
  test05 ( );
  test06 ( );
  test07 ( );
  test08 ( );
  test09 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "LEBESGUE_PRB\n";
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
//    LEBESGUE_TEST01 looks at Chebyshev1 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "chebyshev1";
  double *l;
  string label = "Chebyshev1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST01:\n";
  cout << "  Analyze Chebyshev1 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l,
    "  Chebyshev1 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = chebyshev1 ( n );
  r8vec_print ( n, x, "  Chebyshev1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST02 looks at Chebyshev2 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "chebyshev2";
  double *l;
  string label = "Chebyshev2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST02:\n";
  cout << "  Analyze Chebyshev2 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l,
    "  Chebyshev2 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = chebyshev2 ( n );
  r8vec_print ( n, x, "  Chebyshev2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST03 looks at Chebyshev3 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2014
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "chebyshev3";
  double *l;
  string label = "Chebyshev3 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST03:\n";
  cout << "  Analyze Chebyshev3 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev3 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l, 
    "  Chebyshev3 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = chebyshev3 ( n );
  r8vec_print ( n, x, "  Chebyshev3 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST04 looks at Chebyshev4 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "chebyshev4";
  double *l;
  string label = "Chebyshev4 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST04:\n";
  cout << "  Analyze Chebyshev4 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = chebyshev4 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l, 
    "  Chebyshev4 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = chebyshev4 ( n );
  r8vec_print ( n, x, "  Chebyshev4 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test05 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST05 looks at Equidistant1 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "equidistant1";
  double *l;
  string label = "Equidistant1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST05:\n";
  cout << "  Analyze Equidistant1 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l, 
    "  Equidistant1 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = equidistant1 ( n );
  r8vec_print ( n, x, "  Equidistant1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test06 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST06 looks at Equidistant2 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "equidistant2";
  double *l;
  string label = "Equidistant2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST06:\n";
  cout << "  Analyze Equidistant2 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l, 
    "  Equidistant2 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = equidistant2 ( n );
  r8vec_print ( n, x, "  Equidistant2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test07 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST07 looks at Equidistant3 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "equidistant3";
  double *l;
  string label = "Equidistant3 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST07:\n";
  cout << "  Analyze Equidistant3 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = equidistant3 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l,
    "  Equidistant3 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = equidistant3 ( n );
  r8vec_print ( n, x, "  Equidistant3 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test08 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST08 looks at Fejer 1 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "fejer1";
  double *l;
  string label = "Fejer1 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST08:\n";
  cout << "  Analyze Fejer1 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = fejer1 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l,
    "  Fejer1 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = fejer1 ( n );
  r8vec_print ( n, x, "  Fejer1 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}
//****************************************************************************80

void test09 ( )

//****************************************************************************80
//
//  Purpose:
//
//    LEBESGUE_TEST09 looks at Fejer2 points.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    04 March 2015
//
//  Author:
//
//    John Burkardt
//
{
  string filename = "fejer2";
  double *l;
  string label = "Fejer2 points for N = 11";
  int n;
  int n_max = 11;
  int nfun = 501;
  double *x;
  double *xfun;

  cout << "\n";
  cout << "LEBESGUE_TEST09:\n";
  cout << "  Analyze Fejer2 points.\n";

  xfun = r8vec_linspace_new ( nfun, -1.0, +1.0 );

  l = new double[nfun];

  for ( n = 1; n <= n_max; n++ )
  {
    x = fejer2 ( n );
    l[n-1] = lebesgue_constant ( n, x, nfun, xfun );
    delete [] x;
  }

  r8vec_print ( n_max, l,
    "  Fejer2 Lebesgue constants for N = 1 to 11:" );
//
//  Examine one case more closely.
//
  n = 11;
  x = fejer2 ( n );
  r8vec_print ( n, x, "  Fejer2 points for N = 11" );

  lebesgue_plot ( n, x, nfun, xfun, label, filename );

  delete [] l;
  delete [] x;
  delete [] xfun;

  return;
}

