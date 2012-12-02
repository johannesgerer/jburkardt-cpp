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
//    FSU_CVT_STANDALONE computes a centroidal Voronoi Tessellation dataset.
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
//    * SEED, a seed to use for random number generation;
//    * INIT, initialize the points:
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//    * IT_MAX, the maximum number of iterations;
//    * SAMPLE, how to conduct the sampling:
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//    * SAMPLE_NUM, the number of sampling points;
//    * BATCH, the number of sampling points to create at one time;
//    * OUTPUT, a file in which to store the data;
//
//    To indicate that no further computations are desired, it is 
//    enough to input a nonsensical value, such as -1.
//
//  Modified:
//
//    10 November 2006
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, and Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review, Volume 41, 1999, pages 637-676.
//
{
  int batch;
  char file_out_name[80];
  int i;
  int init;
  char init_string[80];
  int it_num;
  int it_max;
  int n;
  int dim_num;
  double *r;
  bool reset;
  int sample;
  char sample_string[80];
  int sample_num;
  int seed;
  int seed_init;
  char *string;
  bool success;
//
//  Print introduction and options.
//
  timestamp ( );

  cout << "\n";
  cout << "FSU_CVT_STANDALONE\n";
  cout << "  C++ version\n";
  cout << "  Generate a CVT dataset.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << "\n";
  cout << "\n";
  cout << "  This program is meant to be used interactively.\n";
  cout << "  It is also possible to prepare a simple input\n";
  cout << "  file beforehand and use it in batch mode.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * DIM_NUM, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * SEED, a seed to use for random number generation,\n";
  cout << "  * INIT, initialize the points:\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using FORTRAN RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "  * IT_MAX, the maximum number of iterations.\n";
  cout << "  * SAMPLE, how to conduct the sampling.\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using FORTRAN RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "  * SAMPLE_NUM, the number of sample points.\n";
  cout << "  * BATCH, the number of sample points to generate at one time.\n";
  cout << "  * OUTPUT, a file to store the data.\n";
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
    cout << "  (Try '2' if you don't have a preference.)\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> dim_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read DIM_NUM.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input DIM_NUM = " << dim_num << "\n";

    if ( dim_num < 1 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
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
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read N.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input N = " << n << "\n";

    if ( n < 1 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of N = " << n << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  Enter SEED, a seed for the random number generator:\n";
    cout << "  (Try '123456789' if you don't have a preference.)\n";
    cout << "  (Any negative value terminates execution).\n";

    cin >> seed;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read SEED.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input SEED = " << seed << "\n";

    if ( seed < 0 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of SEED = " << seed << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  INIT is the method of initializing the data:\n";
    cout << "\n";
    cout << "  GRID     by picking points from a grid;\n";
    cout << "  HALTON   from a Halton sequence;\n";
    cout << "  RANDOM   using C++ RANDOM function;\n";
    cout << "  UNIFORM  using a simple uniform RNG;\n";
    cout << "\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "  (Try 'RANDOM' if you don't have a preference.)\n";
    cout << "\n";
    cout << "  Enter INIT:\n";
    cin >> init_string;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read INIT.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input INIT = \"" << init_string << "\".\n";

    if ( s_eqi ( init_string, "RANDOM"  ) )
    {
      init = -1;
    }
    else if ( s_eqi ( init_string, "UNIFORM" ) )
    {
      init = 0;
    }
    else if ( s_eqi ( init_string, "HALTON"  ) )
    {
      init = 1;
    }
    else if ( s_eqi ( init_string, "GRID"    ) )
    {
      init = 2;
    }
    else
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of INIT\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  IT_MAX is the maximum number of iterations.\n";
    cout << "\n";
    cout << "  An iteration carries out the following steps:\n";
    cout << "  * the Voronoi region associated with each\n";
    cout << "    generator is estimated by sampling;\n";
    cout << "  * the centroid of each Voronoi region is estimated.\n";
    cout << "  * the generator is replaced by the centroid.\n";
    cout << "\n";
    cout << "  If \"enough\" sampling points are used,\n";
    cout << "  and \"enough\" iterations are taken, this process\n";
    cout << "  will converge.\n";
    cout << "\n";
    cout << "  (Try '50' if you don't have a preference.)\n";
    cout << "  (A negative value terminates execution).\n";
    cout << "\n";
    cout << "  Enter IT_MAX:\n";
    cin >> it_max;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read IT_MAX.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input IT_MAX = " << it_max << "\n";

    if ( it_max < 0 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of IT_MAX = " << it_max << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  SAMPLE is the method of sampling the region:\n";
    cout << "\n";
    cout << "  GRID     by picking points from a grid;\n";
    cout << "  HALTON   from a Halton sequence;\n";
    cout << "  RANDOM   using C++ RANDOM function;\n";
    cout << "  UNIFORM  using a simple uniform RNG;\n";
    cout << "\n";
    cout << "  (Try 'RANDOM' if you don't have a preference.)\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "\n";
    cout << "  Enter SAMPLE:\n";
    cin >> sample_string;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read SAMPLE.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }
    
    cout << "  User input SAMPLE = \"" << sample_string << "\".\n";

    if ( s_eqi ( sample_string, "RANDOM"  ) )
    {
      sample = -1;
    }
    else if ( s_eqi ( sample_string, "UNIFORM" ) )
    {
      sample = 0;
    }
    else if ( s_eqi ( sample_string, "HALTON"  ) )
    {
      sample = 1;
    }
    else if ( s_eqi ( sample_string, "GRID"    ) )
    {
      sample = 2;
    }
    else
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of SAMPLE\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  SAMPLE_NUM is the number of sample points.\n";
    cout << "\n";
    cout << "  The Voronoi regions will be explored by generating\n";
    cout << "  SAMPLE_NUM points.  For each sample point, the\n";
    cout << "  nearest generator is found.  Using more points\n";
    cout << "  gives a better estimate of these regions.\n";
    cout << "\n";
    cout << "  SAMPLE_NUM should be much larger than N, the\n";
    cout << "  number of generators.\n";
    cout << "\n";
    cout << "  (Try '10000' if you don't have a preference.)\n";
    cout << "  (A zero or negative value terminates execution.)\n";
    cout << "\n";
    cout << "  Enter SAMPLE_NUM:\n";

    cin >> sample_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read SAMPLE_NUM.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input SAMPLE_NUM = " << sample_num << "\n";

    if ( sample_num <= 0 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of SAMPLE_NUM = " << sample_num << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }
    cout << "\n";
    cout << "  BATCH is the number of sample points to create at one time.\n";
    cout << "\n";
    cout << "  BATCH should be between 1 and SAMPLE_NUM.\n";
    cout << "\n";
    cout << "  It is FASTER to set BATCH to SAMPLE_NUM;\n";
    cout << "  setting BATCH to 1 requires the least memory.\n";
    cout << "\n";
    cout << "  (Try '" << sample_num << "' if you don't have a preference.)\n";
    cout << "  (A zero or negative value terminates execution.)\n";
    cout << "\n";
    cout << "  Enter BATCH:\n";

    cin >> batch;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read BATCH.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input BATCH = " << batch << "\n";

    if ( batch <= 0 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of BATCH = " << batch << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  OUTPUT is a file in which to store the data.\n";
    cout << "\n";
    cout << "  (Try 'cvt.txt' if you don't have a preference.)\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "\n";
    cout << "  Enter OUTPUT:\n";
    cin >> file_out_name;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "FSU_CVT_STANDALONE - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read OUTPUT.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input OUTPUT = \"" << file_out_name << "\".\n";

    if ( s_len_trim ( file_out_name ) <= 0 )
    {
      cout << "\n";
      cout << "FSU_CVT_STANDALONE\n";
      cout << "  The input value of OUTPUT\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }
//
//  Initialize the data.
//
    r = new double[dim_num*n];

    seed_init = seed;

    fsu_cvt ( dim_num, n, batch, init, sample, sample_num, it_max, 
      &seed, r, &it_num );

    cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, 
      it_num, sample_string, sample_num, r, file_out_name );

    delete [] r;

    cout << "\n";
    cout << "  The data was written to the file \"" << file_out_name << "\"\n";
  }

  cout << "\n";
  timestamp ( );

  return 0;
}
