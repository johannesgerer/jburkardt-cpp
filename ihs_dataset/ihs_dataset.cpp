# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "ihs.H"

int main ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for IHS_DATASET.
//
//  Discussion:
//
//    IHS_DATASET generates an IHS dataset and writes it to a file.
//
//    This program is meant to be used interactively.  It's also
//    possible to prepare a simple input file beforehand and use it
//    in batch mode.
//
//    The program requests input values from the user:
//
//    * DIM_NUM, the spatial dimension,
//    * N, the number of points to generate,
//    * D, the duplication factor,
//    * SEED, a seed for the random number generator.
//
//    The program generates the data, writes it to the file
//
//      ihs_DIM_NUM_N.txt
//
//    where "DIM_NUM" and "N" are the numeric values specified by the user,
//    and then asks the user for more input.   To indicate that no further
//    computations are desired, it is enough to input a nonsensical
//    value, such as -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
{
  int d;
  int dim_num;
  char file_out_name[80];
  int n;
  int seed;
  int seed_init;
  int *r;
  char *string;
//
//  Print introduction and options.
//
  timestamp ( );

  cout << "\n";
  cout << "IHS_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate an IHS dataset.\n";
  cout << "\n";
  cout << "  This program is meant to be used interactively.\n";
  cout << "  It is also possible to prepare a simple input\n";
  cout << "  file beforehand and use it in batch mode.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * DIM_NUM, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * D, the duplication factor,\n";
  cout << "  * SEED, a seed for the random number generator.\n";
  cout << "\n";
  cout << "  The program generates the data, writes it to the file\n";
  cout << "\n";
  cout << "    ihs_DIM_NUM_N.txt\n";
  cout << "\n";
  cout << "  where \"DIM_NUM\" and \"N\" are the numeric values specified\n";
  cout << "  by the user, and then asks the user for more input.\n";
  cout << "\n";
  cout << "  To indicate that no further computations are\n";
  cout << "  desired, it is enough to input a nonsensical value,\n";
  cout << "  such as -1.\n";

  for ( ; ; )
  {
    cout << "  *\n";
    cout << " *\n";
    cout << "*  Ready to generate a new dataset:\n";
    cout << " *\n";
    cout << "  *\n";
    cout << "  Enter DIM_NUM, the spatial dimension:\n";
    cout << "  (2 is a small typical value).\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> dim_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "IHS_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read DIM_NUM.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input DIM_NUM = " << dim_num << "\n";

    if ( dim_num < 1 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of DIM_NUM = " << dim_num << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter N, the number of points to generate:\n";
    cout << "  (10 is a small typical value).\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> n;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read N.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input N = " << n << "\n";

    if ( n < 1 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of N = " << n << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }


    cout << "\n";
    cout << "  Enter D, the duplication factor:\n";
    cout << "  This must be at least 1, but not too large.\n";
    cout << "  (5 is a typical value).\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> d;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read D.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input D = " << d << "\n";

    if ( d < 1 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of D = " << d << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter SEED, a seed for the random number generator:\n";
    cout << "  (Try '123456789' if you have no preference.)\n";
    cout << "  (0 indicates you want a seed to be chosen for you.)\n";
    cout << "  (Any negative value terminates execution).\n";

    cin >> seed;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read SEED.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input SEED = " << seed << "\n";

    if ( seed < 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of SEED = " << seed << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    if ( seed == 0 )
    {
      seed = get_seed ( );
      cout << "\n";
      cout << "  The actual value of SEED will be = " << seed << "\n";
    }
//
//  Initialize the data.
//
    r = new int[dim_num*n];

    seed_init = seed;

    ihs ( dim_num, n, d, &seed, r );

    sprintf ( file_out_name, "ihs_%d_%d.txt", dim_num, n );

    ihs_write ( dim_num, n, d, seed_init, seed, r, file_out_name );

    delete [] r;

    cout << "\n";
    cout << "  The data was written to the file \"" << file_out_name << "\"\n";
  }

  cout << "\n";
  timestamp ( );

  return 0;
}
