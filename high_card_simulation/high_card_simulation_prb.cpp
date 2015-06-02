# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>

using namespace std;

# include "high_card_simulation.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HIGH_CARD_SIMULATION_PRB.
//
//  Discussion:
//
//    HIGH_CARD_SIMULATION_PRB tests the HIGH_CARD_SIMULATION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "HIGH_CARD_SIMULATION_PRB\n";
  cout << "  C++ version.\n";
  cout << "  Test the HIGH_CARD_SIMULATION library.\n";

  test01 ( );
  test02 ( );
  test03 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "HIGH_CARD_SIMULATION_PRB\n";
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
//    TEST01 varies the skip number.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  int deck_size;
  int i;
  double p;
  int seed;
  int skip_num;
  int trial_num;

  deck_size = 100;
  trial_num = 100;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  Estimate the chances of picking the high\n";
  cout << "  card by skipping a given number of initial cards,\n";
  cout << "  using a deck of " << deck_size << " cards\n";
  cout << "  and simulating " << trial_num << " trials.\n";
  cout << "\n";
  cout << "  Skip   Deck Size    Chance of Winning\n";
  cout << "\n";

  seed = 123456789;

  for ( i = 0; i < 10; i++ )
  {
    skip_num = 1 + ( i * deck_size ) / 10;

    p = high_card_simulation ( deck_size, skip_num, trial_num, seed );

    cout << "  " << setw(3) << skip_num
         << "  " << setw(3) << deck_size
         << "  " << setw(14) << p << "\n";
  }
  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 plots the results for a deck of 100 cards.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  string command_filename = "test02_commands.txt";
  ofstream command_unit;
  string data_filename = "test02_data.txt";
  ofstream data_unit;
  int deck_size = 100;
  int i;
  double *p;
  int seed;
  int skip_num;
  int trial_num;

  trial_num = 1000;
  seed = 123456789;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  Compute the changes of picking the high card\n";
  cout << "  after skipping 0 through 99 cards,\n";
  cout << "  using a deck with " << deck_size << " cards\n";
  cout << "  and taking " << trial_num << " trials.\n";

  p = new double[deck_size];

  for ( skip_num = 0; skip_num < deck_size; skip_num++ )
  {
    p[skip_num] = high_card_simulation ( deck_size, skip_num, 
      trial_num, seed );
  }
/*
  Create the data file.
*/
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < deck_size; i++ )
  {
    data_unit << "  " << i
              << "  " << p[i] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
/*
  Write the command file.
*/
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'test02.png'\n";
  command_unit << "set xlabel '<--- Skip --->'\n";
  command_unit << "set ylabel '<--- P(Correct) --->'\n";
  command_unit << 
    "set title 'Estimated Prob of Correct Guess after Skipping K Cards'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'red'\n";

  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] p;

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 plots the results for a deck of 100 cards.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    26 February 2014
//
//  Author:
//
//    John Burkardt
//
{
  string command_filename = "test03_commands.txt";
  ofstream command_unit;
  string data_filename = "test03_data.txt";
  ofstream data_unit;
  int deck_size = 100;
  int i;
  double *p;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  HIGH_CARD_PROBABILITY computes the exact probability of \n";
  cout << "  winning a high card game with a deck of " << deck_size << " cards\n";
  cout << "  assuming we skip the first K cards and select the next card\n";
  cout << "  that is bigger.\n";

  p = high_card_probability ( deck_size );

  cout << "\n";
  cout << "    K   Prob(K)\n";
  cout << "\n";
  for ( i = 0; i < deck_size; i++ )
  {
    cout << "  " << setw(3) << i
         << "  " << setw(8) << p[i] << "\n";
  }
/*
  Create data file.
*/
  data_unit.open ( data_filename.c_str ( ) );
  for ( i = 0; i < deck_size; i++ )
  {
    data_unit << "  " << i
              << "  " << p[i] << "\n";
  }
  data_unit.close ( );
  cout << "\n";
  cout << "  Created graphics data file '" << data_filename << "'\n";
/*
  Create the command file.
*/
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "# " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "# Usage:\n";
  command_unit << "#  gnuplot < " << command_filename << "\n";
  command_unit << "#\n";
  command_unit << "set term png\n";
  command_unit << "set output 'test03.png'\n";
  command_unit << "set xlabel '<--- Skip --->'\n";
  command_unit << "set ylabel '<--- P(Correct) --->'\n";
  command_unit << 
    "set title 'Probability of Correct Guess after Skipping K Cards'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "plot '" << data_filename << "' using 1:2 lw 3 linecolor rgb 'red'\n";

  command_unit.close ( );
  cout << "  Created graphics command file '" << command_filename << "'\n";

  delete [] p;

  return;
}
