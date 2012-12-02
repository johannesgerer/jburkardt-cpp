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

int main ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    FSU_HALTON_STANDALONE generates a Halton dataset and writes it to a file.
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
//  Discussion:
//
//    This program is meant to be used interactively.  It's also
//    possible to prepare a simple input file beforehand and use it
//    in batch mode.
//
//    The program requests input values from the user:
//
//    * DIM_NUM, the spatial dimension,
//    * N, the number of points to generate,
//    * STEP, the index of the first subsequence element to be computed.
//    * SEED(1:DIM_NUM), the Halton sequence index corresponding
//      to STEP = 0.
//    * LEAP(1:DIM_NUM), the successive jumps in the Halton sequence.
//    * BASE(1:DIM_NUM), the Halton bases (usually distinct primes).
//
//    The program generates the data, writes it to the file
//
//      halton_DIM_NUM_N.txt
//
//    where "DIM_NUM" and "N" are the numeric values specified by the user,
//    and then asks the user for more input.   To indicate that no further
//    computations are desired, it is enough to input a nonsensical
//    value, such as -1.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
{
  int *base;
  char command[255];
  char file_out_name[255];
  int i;
  int *leap;
  int n;
  int dim_num;
  double *r;
  int *seed;
  int step;
  char *string;

  timestamp ( );

  cout << "\n";
  cout << "FSU_HALTON_STANDALONE\n";
  cout << "  C++ version\n";
  cout << "  Generate a Halton dataset.\n";
  cout << "\n";
  cout << "  This program is meant to be used interactively.\n";
  cout << "  It is also possible to prepare a simple input\n";
  cout << "  file beforehand and use it in batch mode.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * DIM_NUM, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * STEP, the index of the first subsequence element.\n";
  cout << "  * SEED(1:DIM_NUM),the Halton sequence element\n";
  cout << "    corresponding to STEP = 0\n";
  cout << "  * LEAP(1:DIM_NUM), the succesive jumps in the Halton sequence.\n";
  cout << "  * BASE(1:DIM_NUM), the Halton bases,\n";
  cout << "    usually distinct primes.\n";
  cout << "\n";
  cout << "  The program generates the data, writes it to the file\n";
  cout << "\n";
  cout << "    halton_DIM_NUM_N.txt\n";
  cout << "\n";
  cout << "  where ""DIM_NUM"" and ""N"" are the numeric values specified\n";
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

    cout << "\n";
    cout << "  Enter DIM_NUM, the spatial dimension:\n";
    cout << "  (Try '2' if you don't have a preference.)\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> dim_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read DIM_NUM.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    if ( !halham_dim_num_check ( dim_num ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of DIM_NUM = " << dim_num << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter N, the number of points to generate:\n";
    cout << "  (Try '25' if you don't have a preference.)\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> n;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read N.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    if ( !halham_n_check ( n ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of N = " << n << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter STEP, the index of the first subsequence element:\n";
    cout << "  (Try '0' or '1' if you don't have a preference.\n";
    cout << "  (Any negative value terminates execution).\n";

    cin >> step;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read STEP.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    if ( !halham_step_check ( step ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of STEP = " << step << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    base = new int[dim_num];
    leap = new int[dim_num];
    r = new double[dim_num*n];
    seed = new int[dim_num];

    cout << "\n";
    cout << "  Enter SEED(1:DIM_NUM), the starting element index\n";
    cout << "  for each coordinate\n";
    cout << "  (Try '0 0 ... 0' if you don't have a preference).\n";
    cout << "  (a negative value terminates execution.)\n";

    for ( i = 0; i < dim_num; i++ )
    {
      cin >> seed[i];

      if ( cin.rdstate ( ) )
      {
        cin.clear ( );
  
        cout << "\n";
        cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
        cout << "  An I/O error occurred trying to read SEED[" << i << "].\n";
        cout << "  Abnormal end of execution.\n";
        break;
      }
    }

    i4vec_transpose_print ( dim_num, seed, "  User input:" );

    if ( !halham_seed_check ( dim_num, seed ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of SEED\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter LEAP(1:DIM_NUM), the leaping multiplier\n";
    cout << "  for each coordinate\n";
    cout << "  (Try '1 1 1 ... 1' if you don't have a preference).\n";
    cout << "  (another choice is any prime larger than all the bases).\n";
    cout << "  (any value of 0 or less terminates execution.)\n";

    for ( i = 0; i < dim_num; i++ )
    {
      cin >> leap[i];
      if ( cin.rdstate ( ) )
      {
        cin.clear ( );
  
        cout << "\n";
        cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
        cout << "  An I/O error occurred trying to read LEAP[" << i << "].\n";
        cout << "  Abnormal end of execution.\n";
        break;
      }
    }

    i4vec_transpose_print ( dim_num, leap, "  User input:" );

    if ( !halham_leap_check ( dim_num, leap ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of LEAP\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter BASE(1:DIM_NUM), the base for each coordinate,\n";
    cout << "  usually distinct primes.\n";
    cout << "  (Try '2 3 5 7 11 ...' if you don't have a preference).\n";
    cout << "  (any value of 1 or less terminates execution.)\n";

    for ( i = 0; i < dim_num; i++ )
    {
      cin >> base[i];
      if ( cin.rdstate ( ) )
      {
        cin.clear ( );
  
        cout << "\n";
        cout << "FSU_HALTON_STANDALONE - Fatal error!\n";
        cout << "  An I/O error occurred trying to read BASE[" << i << "].\n";
        cout << "  Abnormal end of execution.\n";
        break;
      }
    }

    i4vec_transpose_print ( dim_num, base, "  User input:" );

    if ( !halton_base_check ( dim_num, base ) )
    {
      cout << "\n";
      cout << "FSU_HALTON_STANDALONE\n";
      cout << "  The input value of BASE\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    fsu_halton ( dim_num, n, step, seed, leap, base, r );

    strcpy ( file_out_name, "halton_" );
    string = i4_to_s ( dim_num );
    strcat ( file_out_name, string );
    delete [] string;
    strcat ( file_out_name, "_" );
    string = i4_to_s ( n );
    strcat ( file_out_name, string );
    delete [] string;
    strcat ( file_out_name, ".txt" );

    halham_write ( dim_num, n, step, seed, leap, base, r, file_out_name );

    delete [] base;
    delete [] leap;
    delete [] r;
    delete [] seed;

    cout << "\n";
    cout << "  The data was written to the file \"" 
      << file_out_name << "\".\n";
  }

  cout << "\n";
  timestamp ( );

  return 0;
}
