# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "vandermonde.hpp"

int main ( );
void bivand1_test ( );
void bivand2_test ( );
void dvand_test ( );
void dvandprg_test ( );
void pvand_test ( );
void pvandprg_test ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
/*
  Purpose:

    MAIN is the main program for VANDERMONDE_PRB.

  Discussion:

    VANDERMONDE_TEST tests the VANDERMONDE library.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    01 May 2014

  Author:

    John Burkardt
*/
{
  timestamp ( );
  cout << "\n";
  cout << "VANDERMONDE_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the VANDERMONDE library.\n";

  bivand1_test ( );
  bivand2_test ( );
  dvand_test ( );
  dvandprg_test ( );
  pvand_test ( );
  pvandprg_test ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "VANDERMONDE_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void bivand1_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BIVAND1_TEST tests BIVAND1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *a;
  double alpha[N] = { 1.0, 2.0, 3.0 };
  double beta[N] = { 10.0, 20.0, 30.0 };
  int n = N;
  int n2;

  cout << "\n";
  cout << "BIVAND1_TEST:\n";
  cout << "  Compute a bidimensional Vandermonde matrix\n";
  cout << "  associated with polynomials of\n";
  cout << "  total degree less than N.\n";

  r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );
  r8vec_print ( n, beta, "  Vandermonde vector BETA:" );

  a = bivand1 ( n, alpha, beta );

  n2 = ( n * ( n + 1 ) ) / 2;
  r8mat_print ( n2, n2, a, "  Bidimensional Vandermonde matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void bivand2_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    BIVAND2_TEST tests BIVAND2.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 May 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 3

  double *a;
  double alpha[N] = { 1.0, 2.0, 3.0 };
  double beta[N] = { 10.0, 20.0, 30.0 };
  int n = N;
  int n2;

  cout << "\n";
  cout << "BIVAND2_TEST:\n";
  cout << "  Compute a bidimensional Vandermonde matrix\n";
  cout << "  associated with polynomials of maximum degree less than N.\n";

  r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );
  r8vec_print ( n, beta, "  Vandermonde vector BETA:" );

  a = bivand2 ( n, alpha, beta );

  n2 = n * n;
  r8mat_print ( n2, n2, a, "  Bidimensional Vandermonde matrix:" );

  delete [] a;

  return;
# undef N
}
//****************************************************************************80

void dvand_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DVAND_TEST tests DVAND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  int n = N;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  cout << "\n";
  cout << "DVAND_TEST:\n";
  cout << "  Solve a Vandermonde linear system A'*x=b\n";

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 ) 
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      alpha = r8vec_uniform_01_new ( n, seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mtv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    delete [] x;

    x = dvand ( n, alpha, b );
    r8vec_print ( n, x, "  Solution X:" );

    delete [] a;
    delete [] alpha;
    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void dvandprg_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    DVANDPRG_TEST tests DVANDPRG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  double *c;
  double *m;
  int n = N;
  int nsub;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  cout << "\n";
  cout << "DVANDPRG_TEST:\n";
  cout << "  Solve a Vandermonde linear system A'*x=b\n";
  cout << "  progressively.\n";
  cout << "  First we use ALPHA = 0, 1, 2, 3, 4.\n";
  cout << "  Then we choose ALPHA at random.\n";

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 ) 
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mtv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    delete [] x;

    x = new double[n];
    c = new double[n];
    m = new double[n];

    for ( nsub = 1; nsub <= n; nsub++ )
    {
      dvandprg ( nsub, alpha, b, x, c, m );
      r8vec_print ( nsub, x, "  Solution X:" );
    }

    delete [] a;
    delete [] alpha;
    delete [] b;
    delete [] c;
    delete [] m;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void pvand_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PVAND_TEST tests PVAND.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 February 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  int n = N;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  cout << "\n";
  cout << "PVAND_TEST:\n";
  cout << "  Solve a Vandermonde linear system A*x=b\n";

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    delete [] x;

    x = pvand ( n, alpha, b );
    r8vec_print ( n, x, "  Solution X:" );

    delete [] a;
    delete [] alpha;
    delete [] b;
    delete [] x;
  }

  return;
# undef N
}
//****************************************************************************80

void pvandprg_test ( )

//****************************************************************************80
//
//  Purpose:
//
//    PVANDPRG_TEST tests PVANDPRG.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 April 2014
//
//  Author:
//
//    John Burkardt
//
{
# define N 5

  double *a;
  double *alpha;
  double alpha1[N] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double *b;
  double *c;
  double *m;
  int n = N;
  int nsub;
  int seed;
  int test;
  double *x;
  double x1[N] = { 5.0, 3.0, 4.0, 1.0, 2.0 };

  cout << "\n";
  cout << "PVANDPRG_TEST:\n";
  cout << "  Solve a Vandermonde linear system A*x=b\n";

  for ( test = 1; test <= 2; test++ )
  {
    if ( test == 1 )
    {
      alpha = r8vec_copy_new ( n, alpha1 );
    }
    else if ( test == 2 )
    {
      seed = 123456789;
      alpha = r8vec_uniform_01_new ( n, seed );
    }

    r8vec_print ( n, alpha, "  Vandermonde vector ALPHA:" );

    a = vand1 ( n, alpha );

    x = r8vec_copy_new ( n, x1 );
    b = r8mat_mv_new ( n, n, a, x );
    r8vec_print ( n, b, "  Right hand side B:" );
    delete [] x;

    x = new double[n];
    c = new double[n];
    m = new double[n];

    for ( nsub = 1; nsub <= n; nsub++ )
    {
      pvandprg ( nsub, alpha, b, x, c, m );
      r8vec_print ( nsub, x, "  Solution X:" );
    }

    delete [] a;
    delete [] alpha;
    delete [] b;
    delete [] c;
    delete [] m;
    delete [] x;
  }

  return;
# undef N
}

