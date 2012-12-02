# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "cvt.H"

int main ( void );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for CVT_DATASET.
//
//  Discussion:
//
//    CVT_DATASET computes a centroidal Voronoi Tessellation dataset.
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
//      ** file, by reading data from file;
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//      ** USER, call the "user" routine;
//    * IT_MAX, the maximum number of iterations;
//    * IT_FIXED, the number of iterative steps to take
//      using a fixed set of sampling points.
//    * SAMPLE, how to conduct the sampling:
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//      ** USER, call the "user" routine;
//    * SAMPLE_NUM, the number of sampling points;
//    * BATCH, the number of sampling points to create at one time;
//    * OUTPUT, a file in which to store the data;
//
//    To indicate that no further computations are desired, it is 
//    enough to input a nonsensical value, such as -1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Qiang Du, Vance Faber, Max Gunzburger,
//    Centroidal Voronoi Tessellations: Applications and Algorithms,
//    SIAM Review,
//    Volume 41, 1999, pages 637-676.
//
{
# define DEBUG 1

  int batch;
  bool comment;
  double energy;
  char file_out_name[80];
  int i;
  int init;
  char init_string[80];
  bool initialize;
  double it_diff;
  int it_fixed;
  int it_num;
  int it_max;
  int n;
  int dim_num;
  double *r;
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
  cout << "CVT_DATASET\n";
  cout << "  C++ version\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  Generate a CVT dataset.\n";
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
  cout << "    ** file, read data from a file;\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using FORTRAN RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "    ** USER, call the \"USER\" routine;\n";
  cout << "  * IT_MAX, the maximum number of iterations.\n";
  cout << "  * IT_FIXED, the number of iterative steps to take\n";
  cout << "    using a fixed set of sampling points.\n";
  cout << "  * SAMPLE, how to conduct the sampling.\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using FORTRAN RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "    ** USER, call the \"USER\" routine;\n";
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
    cout << "  (Try '2' if you have no preference.)\n";
    cout << "  (0 or any negative value terminates execution).\n";

    cin >> dim_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
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
    cout << "  (Try '25' if you have no preference.)\n";
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
    cout << "  Enter SEED, a seed for the random number generator:\n";
    cout << "  (Try '123456789' if you have no preference.)\n";
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

    cout << "\n";
    cout << "  INIT is the method of initializing the data:\n";
    cout << "\n";
    cout << "  file     read data from a file;\n";
    cout << "  GRID     by picking points from a grid;\n";
    cout << "  HALTON   from a Halton sequence;\n";
    cout << "  RANDOM   using C++ RANDOM function;\n";
    cout << "  UNIFORM  using a simple uniform RNG;\n";
    cout << "  USER     call the \"USER\" routine;\n";
    cout << "\n";
    cout << "  (Try 'RANDOM' if you have no preference.)\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "\n";
    cout << "  Enter INIT:\n";
    cin >> init_string;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
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
    else if ( s_eqi ( init_string, "USER"    ) )
    {
      init = 3;
    }
    else if ( 0 < s_len_trim ( init_string ) )
    {
    }
    else
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
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
    cout << "  (Try '50' if you have no preference.)\n";
    cout << "  (A negative value terminates execution).\n";
    cout << "\n";
    cout << "  Enter IT_MAX:\n";
    cin >> it_max;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read IT_MAX.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input IT_MAX = " << it_max << "\n";

    if ( it_max < 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of IT_MAX = " << it_max << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  IT_FIXED is the number of consecutive iterations\n";
    cout << "  to take with a fixed set of sample points.\n";
    cout << "\n";
    cout << "  Setting IT_FIXED to 1 means a new set of sample\n";
    cout << "  points is generated on every iterative step;\n";
    cout << "  Setting IT_FIXED equal to IT_MAX means a single set\n";
    cout << "  of sample points is used for the entire iteration.\n";
    cout << "\n";
    cout << "  Any value between 1 and IT_MAX is reasonable.\n";
    cout << "\n";
    cout << "  (Try " << it_max << " if you have no preference.)\n";
    cout << "  (A 0 or negative value terminates execution).\n";
    cout << "\n";
    cout << "  Enter IT_FIXED:\n";
    cin >> it_fixed;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read IT_FIXED.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input IT_FIXED = " << it_fixed << "\n";

    if ( it_max <= 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of IT_FIXED = " << it_fixed << "\n";
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
    cout << "  USER     call the \"USER\" routine;\n";
    cout << "\n";
    cout << "  (Try 'RANDOM' if you have no preference.)\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "\n";
    cout << "  Enter SAMPLE:\n";
    cin >> sample_string;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
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
    else if ( s_eqi ( sample_string, "USER"    ) )
    {
      sample = 3;
    }
    else
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
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
    cout << "  (Try '10000' if you have no preference.)\n";
    cout << "  (A zero or negative value terminates execution.)\n";
    cout << "\n";
    cout << "  Enter SAMPLE_NUM:\n";

    cin >> sample_num;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read SAMPLE_NUM.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input SAMPLE_NUM = " << sample_num << "\n";

    if ( sample_num <= 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
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
    cout << "  (Try " << i_min ( sample_num, 1000 ) 
         << " if you have no preference.)\n";
    cout << "  (A zero or negative value terminates execution.)\n";
    cout << "\n";
    cout << "  Enter BATCH:\n";

    cin >> batch;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read BATCH.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input BATCH = " << batch << "\n";

    if ( batch <= 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of BATCH = " << batch << "\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }

    cout << "\n";
    cout << "  OUTPUT is a file in which to store the data.\n";
    cout << "\n";
    cout << "  (Try 'cvt.txt' if you have no preference.)\n";
    cout << "  (A blank value terminates execution).\n";
    cout << "\n";
    cout << "  Enter OUTPUT:\n";
    cin >> file_out_name;

    if ( cin.rdstate ( ) )
    {
      cin.clear ( );
  
      cout << "\n";
      cout << "CVT_DATASET - Fatal error!\n";
      cout << "  An I/O error occurred while trying to read OUTPUT.\n";
      cout << "  Abnormal end of execution.\n";
      break;
    }

    cout << "  User input OUTPUT = \"" << file_out_name << "\".\n";

    if ( s_len_trim ( file_out_name ) <= 0 )
    {
      cout << "\n";
      cout << "CVT_DATASET\n";
      cout << "  The input value of OUTPUT\n";
      cout << "  is interpreted as a request for termination.\n";
      cout << "  Normal end of execution.\n";
      break;
    }
//
//  Initialize the data.
//
    r = new double[dim_num*n];

    if ( init == 4 )
    {
      data_read ( sample_string, dim_num, n, r );
    }

    seed_init = seed;

    cvt ( dim_num, n, batch, init, sample, sample_num, it_max, it_fixed,
      &seed, r, &it_num, &it_diff, &energy );

    comment = false;

    cvt_write ( dim_num, n, batch, seed_init, seed, init_string, it_max, 
      it_fixed, it_num, it_diff, energy, sample_string, sample_num, r,
      file_out_name, comment );

    delete [] r;

    cout << "\n";
    cout << "  The data was written to the file \"" << file_out_name << "\"\n";
  }

  cout << "\n";
  timestamp ( );

  return 0;
# undef DEBUG
}
