# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "subset_sum_serial.hpp"

int main ( );
void test01 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for SUBSET_SUM_SERIAL_PRB.
//
//  Discussion:
//
//    SUBSET_SUM_SERIAL_TEST tests the SUBSET_SUM_SERIAL library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "SUBSET_SUM_SERIAL_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the SUBSET_SUM_SERIAL library.\n";

  test01 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "SUBSET_SUM_SERIAL_PRB\n";
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
//    SUBSET_SUM_SERIAL_TEST01 tests the SUBSET_SUM_SERIAL program.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 July 2013
//
//  Author:
//
//    John Burkardt
//
{
  int *choice;
  int i;
  int n = 21;
  int target;
  int weight[21] = {
    518533, 1037066, 2074132, 1648264, 796528, 
   1593056,  686112, 1372224,  244448, 488896, 
    977792, 1955584, 1411168,  322336, 644672, 
   1289344,   78688,  157376,  314752, 629504, 
   1259008 };
  int w_sum;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Test the SUBSET_SUM_SERIAL function, which looks for a selection\n";
  cout << "  from a set of weights that adds up to a given target.\n";
//
//  Define the problem data.
//
  target = 2463098;
  cout << "\n";
  cout << "  Target value:\n";
  cout << "  " << target << "\n";
 
  cout << "\n";
  cout << "   I      W(I)\n";
  cout << "\n";
  for ( i = 0; i < n; i++ )
  {
    cout << "  " << setw(2) << i
         << "  " << setw(8) << weight[i] << "\n";
  }

  choice = subset_sum_serial ( n, weight, target );

  if ( choice[0] == -1 )
  {
    cout << "\n";
    cout << "  No solution was found.\n";
  }
  else
  {
    cout << "\n";
    cout << "   I*     W*\n";
    cout << "\n";
    w_sum = 0;
    for ( i = 0; i < n; i++ )
    {
      if ( choice[i] == 1 )
      {
        w_sum = w_sum + weight[i];
        cout << "  " << setw(2) << i
             << "  " << setw(8) << weight[i] << "\n";
      }
    }
    cout << "\n";
    cout << "  Sum:    " << w_sum << "\n";
    cout << "  Target: " << target << "\n";
  }

  delete [] choice;

  return;
}

