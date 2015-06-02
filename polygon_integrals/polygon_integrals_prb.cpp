# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "polygon_integrals.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for POLYGON_INTEGRALS_PRB.
//
//  Discussion:
//
//    POLYGON_INTEGRALS_PRB tests POLYGON_INTEGRALS.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "POLYGON_INTEGRALS_PRB:\n";
  cout << "  C++ version\n";
  cout << "  Test the POLYGON_INTEGRALS library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "POLYGON_INTEGRALS_PRB:\n";
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
//    TEST01 carries out a test on a rectangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  double alpha_exact[6] = {
    1.0, 
    5.0, 4.0, 
    30.66666666666667, 22.0, 18.66666666666666 };
  double alpha_pq;
  int k;
  double mu_exact[6] = {
    1.0, 
    0.0, 0.0, 
    5.666666666666667, 2.0, 2.666666666666667 };
  double mu_pq;
  int n = 4;
  double nu_exact[6] = {
    40.0, 
    200.0, 160.0, 
    1226.66666666666667, 880.0, 746.66666666666666 };
  double nu_pq;
  int p;
  int q;
  int s;
  double x[4] = {
    2.0, 10.0, 8.0, 0.0 };
  double y[4] = {
    0.0,  4.0, 8.0, 4.0 };

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Check normalized moments of a rectangle.\n";
  cout << "\n";
  cout << "   P   Q             Nu(P,Q)\n";
  cout << "            Computed         Exact\n";
  cout << "\n";
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      nu_pq = moment ( n, x, y, p, q );
      cout << "  " << setw(2) << p
           << "  " << setw(2) << q
           << "  " << setw(14) << nu_pq
           << "  " << setw(14) << nu_exact[k] << "\n";
      k = k + 1;
    }
  }

  cout << "\n";
  cout << "   P   Q           Alpha(P,Q)\n";
  cout << "            Computed         Exact\n";
  cout << "\n";
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      alpha_pq = moment_normalized ( n, x,y, p, q );
      cout << "  " << setw(2) << p
           << "  " << setw(2) << q
           << "  " << setw(14) << alpha_pq
           << "  " << setw(14) << alpha_exact[k] << "\n";
      k = k + 1;
    }
  }

  cout << "\n";
  cout << "   P   Q             Mu(P,Q)\n";
  cout << "            Computed         Exact\n";
  cout << "\n";
  k = 0;
  for ( s = 0; s <= 2; s++ )
  {
    for ( p = s; 0 <= p; p-- )
    {
      q = s - p;
      mu_pq = moment_central ( n, x, y , p, q );
      cout << "  " << setw(2) << p
           << "  " << setw(2) << q
           << "  " << setw(14) << mu_pq
           << "  " << setw(14) << mu_exact[k] << "\n";
      k = k + 1;
    }
  }

  return;
}
