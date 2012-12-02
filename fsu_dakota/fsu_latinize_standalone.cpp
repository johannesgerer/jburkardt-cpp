# ifdef ANSI_HEADERS
#   include <cstdlib>
#   include <cmath>
#   include <ctime>
# else
#   include <stdlib.h>
#   include <math.h>
#   include <time.h>
# endif

# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "fsu.H"

int main ( int argc, char *argv[] );
void latinize_handle ( char *input_filename );

//****************************************************************************80

int main ( int argc, char *argv[] )

//****************************************************************************80
//
//  Purpose:
//
//    FSU_LATINIZE_STANDALONE is the main routine of a program to "latinize" a dataset.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Usage:
//
//    fsu_latinize_standalone input_filename
//
{
  int i;
  char input_filename[80];

  timestamp ( );

  cout << "\n";
  cout << "FSU_LATINIZE_STANDALONE (C++ version)\n";
  cout << "  Read a dataset of N points in M dimensions,\n";
  cout << "  modify it into a Latin hypercube,\n";
  cout << "  write the modified dataset to a file.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";

  if ( argc <= 1 ) 
  {
    cout << "\n";
    cout << "FSU_LATINIZE_STANDALONE:\n";
    cout << "  Please enter the name of a file to be analyzed.\n";

    cin.getline ( input_filename, sizeof ( input_filename ) );

    latinize_handle ( input_filename );

  }
  else 
  {
    for ( i = 1; i < argc; i++ ) 
    {
      latinize_handle ( argv[i] );
    }
  } 

  cout << "\n";
  cout << "FSU_LATINIZE_STANDALONE\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void latinize_handle ( char *input_filename )

//****************************************************************************80
//
//  Purpose:
//
//    LATINIZE_HANDLE handles a single file.
//
//  License:
//
//    Copyright (C) 2004  John Burkardt and Max Gunzburger
//
//    This library is free software; you can redistribute it and/or
//    modify it under the terms of the GNU Lesser General Public
//    License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//
//    This library is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    Lesser General Public License for more details.
//
//    You should have received a copy of the GNU Lesser General Public
//    License along with this library; if not, write to the Free Software
//    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//  Local parameters:
//
//    Local, int DIM_NUM, the spatial dimension of the point set.
//
//    Local, int N, the number of points.
//
//    Local, double Z[DIM_NUM*N], the point set.
//
//    Local, int NS, the number of sample points.
//
{
  int dim_num;
  int n;
  ofstream output;
  char *output_filename;
  double *table;
//
//  Need to create the output file name from the input filename.
//
  output_filename = file_name_ext_swap ( input_filename, "latin.txt" );

  dtable_header_read ( input_filename, &dim_num, &n );

  cout << "\n";
  cout << "  Read the header of \"" << input_filename << "\".\n";
  cout << "\n";
  cout << "  Spatial dimension DIM_NUM = " << dim_num << "\n";
  cout << "  Number of points N  = " << n << "\n";

  table = dtable_data_read ( input_filename, dim_num, n );

  cout << "\n";
  cout << "  Read the data in \"" << input_filename << "\".\n";

  r8mat_transpose_print_some ( dim_num, n, table, 1, 1, 5, 5, 
    "  Small portion of data read from file:" );

  fsu_latinize ( dim_num, n, table );

  cout << "\n";
  cout << "  Latinized the data.\n";

  r8mat_transpose_print_some ( dim_num, n, table, 1, 1, 5, 5, 
    "  Small portion of Latinized data:" );
//
//  Write the data to a file.
//
  output.open ( output_filename );

  if ( !output )
  {
    cout << "\n";
    cout << "LATINIZE_HANDLE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  dtable_header_write ( dim_num, n, output_filename, output );

  output << "#  This data was created\n";
  output << "#  the C++ FSU_LATINIZE_STANDALONE program.\n";
  output << "#  The data was read from " << input_filename 
         << " by DTABLE_DATA_READ.CC\n";
  output << "#  The data was latinized by FSU_LATINIZE.CC.\n";
  output << "#  The data was written to " << output_filename 
         << " by DTABLE_DATA_WRITE.CC.\n";
  output << "#\n";

  dtable_data_write ( dim_num, n, table, output );

  output.close ( );

  cout << "\n";
  cout << "  Wrote the latinized data to \"" << output_filename << "\".\n";

  delete [] table;

  return;
}
