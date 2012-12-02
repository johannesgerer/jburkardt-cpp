# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "hex_grid.H"

int main ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for HEX_GRID_DATASET.
//
//  Discussion:
//
//    HEX_GRID_DATASET generates a hexagonal grid dataset and writes it to a file.
//
//    This program is meant to be used interactively.  It's also
//    possible to prepare a simple input file beforehand and use it
//    in batch mode.
//
//    The program requests input values from the user:
//
//    * X1, Y1, the lower left corner of a bounding box.
//    * X2, Y2, the upper right corner of a bounding box.
//    * NODES_PER_LAYER, the number of nodes per layer.
//
//    The program will now define the dataset, and write it to a file.
//
//    The program will now request that you type "Y" if you want to 
//    set up another dataset.  Otherwise the program terminates.
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
# define NDIM 2

  double box[NDIM*2];
  char c;
  double hx;
  double hy;
  int layers;
  int n;
  int nodes_per_layer;
  char output_file_name[255];
  double *p;
  char *string;
  double x;
  double y;

  timestamp ( );

  cout << "\n";
  cout << "HEX_GRID_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a hexagonal grid dataset.\n";
  cout << "\n";
  cout << "  This program is meant to be used interactively.\n";
  cout << "  It is also possible to prepare a simple input\n";
  cout << "  file beforehand and use it in batch mode.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * X1, Y1, the lower left corner of the region,\n";
  cout << "  * X2, Y2, the upper right corner of the region,\n";
  cout << "  * NODES_PER_LAYER, the number of nodes per layer,\n";
  cout << "\n";
  cout << "  After the dataset of nodes is computed, it is\n";
  cout << "  written to a file, and another dataset may be made.\n";

  for ( ; ; )
  {
    cout << "  *\n";
    cout << " *\n";
    cout << "*  Ready to generate a new dataset:\n";
    cout << " *\n";
    cout << "  *\n";
    cout << "\n";
    cout << "  Enter X1, Y1, the lower left corner of the region:\n";
    cout << "  (Try '0 0' if you do not have a preference.)\n";
    cout << "  (Enter \"*\" or \"QUIT\" to terminate execution.)\n";

    cin >> x >> y;

    if ( !cin.good ( ) )
    {
      cin.clear ( );
      cout << "\n";
      cout << "HEX_GRID_DATASET - Warning!\n";
      cout << "  Terminating abnormally because of an I/O error\n";
      cout << "  while expecting input for X1, Y1.\n";
      break;
    }

    box[0+0*2] = x;
    box[1+0*2] = y;

    cout << "  User input X1 = " << x << "  "
         <<              "Y1 = " << y << "\n";

    cout << "\n";
    cout << "  Enter X2, Y2, the upper right corner of the region:\n";
    cout << "  (Try '10 10' if you do not have a preference.)\n";
    cout << "  (Enter \"*\" or \"QUIT\" to terminate execution.)\n";

    cin >> x >> y;

    if ( !cin.good ( ) )
    {
      cin.clear ( );
      cout << "\n";
      cout << "HEX_GRID_DATASET - Warning!\n";
      cout << "  Terminating abnormally because of an I/O error\n";
      cout << "  while expecting input for X2, Y2.\n";
      break;
    }

    cout << "  User input X2 = " << x << "  "
         <<              "Y2 = " << y << "\n";

    box[0+1*2] = x;
    box[1+1*2] = y;

    cout << "\n";
    cout << "  Enter NODES_PER_LAYER, the number of nodes in a layer.\n";
    cout << "  (Try '10' if you do not have a preference.)\n";
    cout << "  (1 or any smaller value terminates execution).\n";

    cin >> nodes_per_layer;

    if ( !cin.good ( ) )
    {
      cin.clear ( );
      cout << "\n";
      cout << "HEX_GRID_DATASET - Warning!\n";
      cout << "  Terminating abnormally because of an I/O error\n";
      cout << "  while expecting input for NODES_PER_LAYER.\n";
      break;
    }

    cout << "  User input NODES_PER_LAYER = " << nodes_per_layer << "\n";

    if ( nodes_per_layer <= 1 )
    {
      cout << "\n";
      cout << "HEX_GRID_DATASET\n";
      cout << "  The input value of NODES_PER_LAYER = " 
           << nodes_per_layer << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    layers = hex_grid_layers ( nodes_per_layer, box );

    cout << "\n";
    cout << "  The number of layers will be " << layers << "\n";

    hex_grid_h ( nodes_per_layer, box, &hx, &hy );

    cout << "  The X spacing will be " << hx << "\n";
    cout << "  The Y spacing will be " << hy << "\n";

    n = hex_grid_n ( nodes_per_layer, box );

    cout << "  The number of nodes " << n << "\n";

    p = hex_grid_points ( nodes_per_layer, layers, n, box );

    sprintf ( output_file_name, "hex_grid_%d_%d_%d.txt", 
      nodes_per_layer, layers, n );

    hex_grid_write ( n, nodes_per_layer, layers, hx, hy, box, 
      p, output_file_name );

    delete [] p;

    cout << "\n";
    cout << "  The data was written to the file \"" 
         << output_file_name << "\".\n";

    cout << "\n";
    cout << "  Enter \"Y\" if you want to define another dataset.\n";

    cin >> c;

    if ( !cin.good ( ) )
    {
      cin.clear ( );
      cout << "\n";
      cout << "HEX_GRID_DATASET - Warning!\n";
      cout << "  Terminating abnormally because of an I/O error\n";
      cout << "  while expecting user input.\n";
      break;
    }

    if ( c != 'y' && c != 'Y' )
    {
      cout << "\n";
      cout << "HEX_GRID_DATASET:\n";
      cout << "  Normal end of execution.\n";
      break;
    }

  }

  cout << "\n";
  timestamp ( );

  return 0;
# undef NDIM
}
