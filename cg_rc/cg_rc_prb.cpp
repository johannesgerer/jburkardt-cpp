# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

# include "cg_rc.hpp"

int main ( );
void test01 ( );
void test02 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CG_RC_PRB.
//
//  Discussion:
//
//    CG_RC_PRB tests the CG_RC library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "CG_RC_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the CG_RC library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "CG_RC_PRB:\n";
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
//    TEST01 uses CG_RC for the simple 1, -2, 1 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double angle;
  double *b;
  double bnrm2;
  double err;
  int i;
  int it;
  int it_max;
  int job;
  int n = 21;
  double *p;
  double pi = 3.141592653589793;
  double *q;
  double *r;
  double rnrm2;
  double t;
  double tol;
  double *x;
  double *x_exact;
  double *z;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Use CG_RC on the 1, -2, 1 matrix.\n";
//
//  In order to specify the right hand side, pick an exact solution,
//  and multiply by the matrix.
//
  x_exact = new double[n];
  for ( i = 0; i < n; i++ )
  {
    angle = 2.0 * pi * ( double ) ( i ) / ( double ) ( n - 1 );
    x_exact[i] = sin ( angle );
  }

  b = new double[n];

  for ( i = 0; i < n; i++ )
  {
    b[i] = - 2.0 * x_exact[i];
  }
  for ( i = 0; i < n - 1; i++ )
  {
    b[i] = b[i] + x_exact[i+1];
  }
  for ( i = 1; i < n; i++ )
  {
    b[i] = b[i] + x_exact[i-1];
  }
//
//  Here is the initial guess for the solution.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }
//
//  Parameters for the stopping test.
//
  it = 0;
  it_max = 30;
  tol = 1.0E-05;
  bnrm2 = 0.0;
  for ( i = 0; i < n; i++ )
  {
    bnrm2 = bnrm2 + b[i] * b[i];
  }
  bnrm2 = sqrt ( bnrm2 );
//
//  Set parameters for the CG_RC code.
//
  r = new double[n];
  z = new double[n];
  p = new double[n];
  q = new double[n];

  job = 1;
//
//  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
//
  for ( ; ; )
  {
    job = cg_rc ( n, b, x, r, z, p, q, job );
//
//  Compute q = A * p.
//
    if ( job == 1 )
    {
      for ( i = 0; i < n; i++ )
      {
        q[i] = - 2.0 * p[i];
      }
      for ( i = 0; i < n - 1; i++ )
      {
        q[i] = q[i] + p[i+1];
      }
      for ( i = 1; i < n; i++ )
      {
        q[i] = q[i] + p[i-1];
      }
    }
//
//  Solve M * z = r.
//
    else if ( job == 2 )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = r[i] / ( - 2.0 );
      }
    }
//
//  Compute r = r - A * x.
//
    else if ( job == 3 )
    {
      for ( i = 0; i < n; i++ )
      {
        r[i] = r[i] + 2.0 * x[i];
      }
      for ( i = 0; i < n - 1; i++ )
      {
        r[i] = r[i] - x[i+1];
      }
      for ( i = 1; i < n; i++ )
      {
        r[i] = r[i] - x[i-1];
      }
    }
//
//  Stopping test on R.
//
    else if ( job == 4 )
    {
      rnrm2 = 0.0;
      for ( i = 0; i < n; i++ )
      {
        rnrm2 = rnrm2 + r[i] * r[i];
      }
      rnrm2 = sqrt ( rnrm2 );

      if ( bnrm2 == 0.0 )
      {
        if ( rnrm2 <= tol )
        {
          break;
        }
      }
      else
      {
        if ( rnrm2 <= tol * bnrm2 )
        {
          break;
        }
      }

      it = it + 1;

      if ( it_max <= it )
      {
        cout << "\n";
        cout << "  Iteration limit exceeded.\n";
        cout << "  Terminating early.\n";
        break;
      }
    }
    job = 2;
  }
  
  cout << "\n";
  cout << "  Number of iterations was " << it << "\n";
  cout << "  Estimated error is " << rnrm2 << "\n";
  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = fabs ( x_exact[i] - x[i] );
    if ( err < t )
    {
      err = t;
    }
  }
  cout << "  Loo error is " << err << "\n";

  cout << "\n";
  cout << "     I      X(I)         X_EXACT(I)        B(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << x_exact[i]
         << "  " << setw(14) << b[i] << "\n";
  }

  delete [] b;
  delete [] p;
  delete [] q;
  delete [] r;
  delete [] x;
  delete [] x_exact;
  delete [] z;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests CG_RC with the Wathen matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 January 2013
//
//  Author:
//
//    John Burkardt
//
{
  double *a;
  double *ax;
  double *b;
  double bnrm2;
  double err;
  int i;
  int ii;
  int it;
  int it_max;
  int job;
  int n;
  int nx;
  int ny;
  double *p;
  double *q;
  double *r;
  double rnrm2;
  int seed;
  double t;
  double tol;
  double *x;
  double *x_exact;
  double *z;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Use CG_RC to solve a linear system\n";
  cout << "  involving the Wathen matrix.\n";

  nx = 5;
  ny = 4;
 
  cout << "\n";
  cout << "  NX = " << nx << "\n";
  cout << "  NY = " << ny << "\n";
  cout << "  N  = " << n << "\n";

  n = wathen_order ( nx, ny );

  a = wathen ( nx, ny, n );

  seed = 123456789;
  x_exact = r8vec_uniform_01_new ( n, seed );

  b = new double[n];
  r8mat_mv ( n, n, a, x_exact, b );
//
//  Here is the initial guess for the solution.
//
  x = new double[n];
  for ( i = 0; i < n; i++ )
  {
    x[i] = 0.0;
  }

  ax = new double[n];
//
//  Parameters for the stopping test.
//
  it = 0;
  it_max = 30;
  tol = 1.0E-05;
  bnrm2 = 0.0;
  for ( i = 0; i < n; i++ )
  {
    bnrm2 = bnrm2 + b[i] * b[i];
  }
  bnrm2 = sqrt ( bnrm2 );
//
//  Set parameters for the CG_RC code.
//
  r = new double[n];
  z = new double[n];
  p = new double[n];
  q = new double[n];
  job = 1;
//
//  Repeatedly call the CG_RC code, and on return, do what JOB tells you.
//
  for ( ; ; )
  {
    job = cg_rc ( n, b, x, r, z, p, q, job );
//
//  Compute q = A * p.
//
    if ( job == 1 )
    {
      r8mat_mv ( n, n, a, p, q );
    }
//
//  Solve M * z = r.
//
    else if ( job == 2 )
    {
      for ( i = 0; i < n; i++ )
      {
        z[i] = r[i] / a[i+i*n];
      }
    }
//
//  Compute r = r - A * x.
//
    else if ( job == 3 )
    {
      r8mat_mv ( n, n, a, x, ax );
      for ( i = 0; i < n; i++ )
      {
        r[i] = r[i] - ax[i];
      }
    }
//
//  Stopping test.
//
    else if ( job == 4 )
    {
      rnrm2 = 0.0;
      for ( i = 0; i < n; i++ )
      {
        rnrm2 = rnrm2 + r[i] * r[i];
      }
      rnrm2 = sqrt ( rnrm2 );

      if ( bnrm2 == 0.0 )
      {
        if ( rnrm2 <= tol )
        {
          break;
        }
      }
      else
      {
        if ( rnrm2 <= tol * bnrm2 )
        {
          break;
        }
      }

      it = it + 1;

      if ( it_max <= it )
      {
        cout << "\n";
        cout << "  Iteration limit exceeded.\n";
        cout << "  Terminating early.\n";
        break;
      }
    }
    job = 2;
  }
  
  cout << "\n";
  cout << "  Number of iterations was " << it << "\n";
  cout << "  Estimated error is " << rnrm2 << "\n";
  err = 0.0;
  for ( i = 0; i < n; i++ )
  {
    t = fabs ( x_exact[i] - x[i] );
    if ( err < t )
    {
      err = t;
    }
  }
  cout << "  Loo error is " << err << "\n";

  cout << "\n";
  cout << "     I      X(I)         X_EXACT(I)        B(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i]
         << "  " << setw(14) << x_exact[i]
         << "  " << setw(14) << b[i] << "\n";
  }

  delete [] a;
  delete [] ax;
  delete [] b;
  delete [] p;
  delete [] q;
  delete [] r;
  delete [] x;
  delete [] x_exact;
  delete [] z;

  return;
}

