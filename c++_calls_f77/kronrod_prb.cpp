# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

using namespace std;

int main ( );
void test01 ( );
void test02 ( );
//
//  Interfaces to FORTRAN77 routines.
//  The names now have underscores appended, and
//  scalar parameters are passed by reference.
//  The EXTERN statement used here stops C++ from "name mangling", which
//  would have destroyed the correspondence between the symbolic names
//  of the subroutines as recorded in the C++ and FORTRAN77 compiled codes.
//
extern "C" {
  void kronrod_ ( int *n, double *eps, double x[], double w1[],
    double w2[], int *ier );
  void timestamp_ ( void );
}
//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for KRONROD_PRB.
//
//  Discussion:
//
//    KRONROD_PRB tests the KRONROD library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  timestamp_ ( );
  cout << "\n";
  cout << "KRONROD_PRB:\n";
  cout << "  Combination C++ and FORTRAN77 version.\n";
  cout << "  Demonstrate how a C++ main program can call\n";
  cout << "  a FORTRAN77 library.\n";

  test01 ( );
  test02 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "KRONROD_PRB:\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp_ ( );

  return 0;
}
//****************************************************************************80

void test01 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests the code for the odd case N = 3.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  double eps;
  int i;
  int i2;
  int ier;
  int n = 3;
  double s;
  double *w1;
  double *w2;
  double wg[3] = {
    0.555555555555555555556,
    0.888888888888888888889,
    0.555555555555555555556 };
  double wk[7] = {
    0.104656226026467265194,
    0.268488089868333440729,
    0.401397414775962222905,
    0.450916538658474142345,
    0.401397414775962222905,
    0.268488089868333440729,
    0.104656226026467265194 };
  double *x;
  double xg[3] = {
   -0.77459666924148337704,
    0.0,
    0.77459666924148337704 };
  double xk[7]= {
   -0.96049126870802028342,
   -0.77459666924148337704,
   -0.43424374934680255800,
    0.0,
    0.43424374934680255800,
    0.77459666924148337704,
    0.96049126870802028342 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Request KRONROD to compute the Gauss rule\n";
  cout << "  of order 3, and the Kronrod extension of\n";
  cout << "  order 3+4=7.\n";
  cout << "\n";
  cout << "  Compare to exact data.\n";

  eps = 0.000001;
  w1 = new double[n+1];
  w2 = new double[n+1];
  x = new double[n+1];
//
//  To call the FORTRAN77 subroutine KRONROD, we must append
//  an underscore to the name, and pass scalars by reference.
//
  kronrod_ ( &n, &eps, x, w1, w2, &ier );

  cout << "\n";
  cout << "  KRONROD returns 3 vectors of length " << n + 1 << "\n";
  cout << "\n";
  cout << "     I      X               WK              WG\n";
  cout << "\n";
  for ( i = 1; i <= n + 1; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i-1]
         << "  " << setw(14) << w1[i-1]
         << "  " << setw(14) << w2[i-1] << "\n";
  }

  cout << "\n";
  cout << "               Gauss Abscissas\n";
  cout << "            Exact           Computed\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
      s = +1.0;
    }
    cout << "  " << setw(4) << i
         << "  " << setw(14) << xg[i-1]
         << "  " << setw(14) << s * x[i2-1] << "\n";
  }
  cout << "\n";
  cout << "               Gauss Weights\n";
  cout << "            Exact           Computed\n";
  cout << "\n";
  for ( i = 1; i <= n; i++ )
  {
    if ( 2 * i <= n + 1 )
    {
      i2 = 2 * i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - 2 * i;
    }
    cout << "  " << setw(4) << i
         << "  " << setw(14) << wg[i-1]
         << "  " << setw(14) << w2[i2-1] << "\n";
  }

  cout << "\n";
  cout << "             Gauss Kronrod Abscissas\n";
  cout << "            Exact           Computed\n";
  cout << "\n";
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
      s = -1.0;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
      s = +1.0;
    }
    cout << "  " << setw(4) << i
         << "  " << setw(14) << xk[i-1]
         << "  " << setw(14) << s * x[i2-1] << "\n";
  }
  cout << "\n";
  cout << "             Gauss Kronrod Weights\n";
  cout << "            Exact           Computed\n";
  cout << "\n";
  for ( i = 1; i <= 2 * n + 1; i++ )
  {
    if ( i <= n + 1 )
    {
      i2 = i;
    }
    else
    {
      i2 = 2 * ( n + 1 ) - i;
    }
    cout << "  " << setw(4) << i
         << "  " << setw(14) << wk[i-1]
         << "  " << setw(14) << w1[i2-1] << "\n";
  }

  delete [] w1;
  delete [] w2;
  delete [] x;

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests the code for the even case N = 4.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2010
//
//  Author:
//
//    John Burkardt
//
{
  double eps;
  int i;
  int i2;
  int ier;
  int n = 4;
  double s;
  double *w1;
  double *w2;
  double *x;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Request KRONROD to compute the Gauss rule\n";
  cout << "  of order 4, and the Kronrod extension of\n";
  cout << "  order 4+5=9.\n";

  eps = 0.000001;
  w1 = new double[n+1];
  w2 = new double[n+1];
  x = new double[n+1];
//
//  To call the FORTRAN77 subroutine KRONROD, we must append
//  an underscore to the name, and pass scalars by reference.
//
  kronrod_ ( &n, &eps, x, w1, w2, &ier );

  cout << "\n";
  cout << "  KRONROD returns 3 vectors of length " << n + 1 <<"\n";
  cout << "\n";
  cout << "     I      X               WK              WG\n";
  cout << "\n";
  for ( i = 1; i <= n + 1; i++ )
  {
    cout << "  " << setw(4) << i
         << "  " << setw(14) << x[i-1]
         << "  " << setw(14) << w1[i-1]
         << "  " << setw(14) << w2[i-1] << "\n";
  }

  delete [] w1;
  delete [] w2;
  delete [] x;

  return;
}
