# include <cstdlib>
# include <cmath>
# include <ctime>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <cstring>

using namespace std;

int main ( void );
char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );
double cluster_energy ( int dim_num, int n, double cell_generator[],
  int sample_num, int sample_function_cvt, int *seed );
void cvt_iteration ( int m, int n, double generator[], int sample_num,
 int sample_function_cvt, int *seed, double *change_l2 );
int find_closest ( int m, int n, double x[], double generator[] );
int i4_max ( int i1, int i2 );
int i4_min ( int i1, int i2 );
void i4_to_halton ( int seed, int base[], int ndim, double r[] );
void lcvt_write ( int dim_num, int n, int seed_start, int sample_function_init,
  char* file_in_name, int sample_function_cvt, int sample_num, int cvt_it,
  double cvt_energy, int latin_it, double latin_energy, double cell_generator[],
  char *file_out_name );
int prime ( int n );
double r8_epsilon ( );
double r8_uniform_01 ( int *seed );
void r8mat_latinize ( int m, int n, double table[] );
double *r8table_data_read ( char *input_filename, int m, int n );
int *r8vec_sort_heap_index_a ( int n, double a[] );
void region_sampler ( int m, int n, int n_total, double x[],
  int sample_function, bool reset, int *seed );
bool s_eqi ( char *s1, char *s2 );
int s_len_trim ( char* s );
double s_to_r8 ( char *s, int *lchar, bool *error );
bool s_to_r8vec ( char *s, int n, double rvec[] );
void timestamp ( );
void tuple_next_fast ( int m, int n, int rank, int x[] );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for LCVT_DATASET.
//
//  Discussion:
//
//    LCVT_DATASET computes a Latinized CVT dataset and writes it to a file.
//
//    This program is meant to be used interactively.  It's also
//    possible to prepare a simple input file beforehand and use it
//    in batch mode.
//
//    The program requests input values from the user:
//
//    * DIM_NUM, the spatial dimension;
//    * N, the number of points to generate;
//    * SEED_INIT, a seed to use for random number generation;
//    * INIT, initialize the points:
//      ** file, by reading data from file;
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//      ** USER, call the "user" routine;
//    * CVT_IT_NUM, the maximum number of iterations;
//    * SAMPLE, how to conduct the sampling:
//      ** GRID, picking points from a grid;
//      ** HALTON, from a Halton sequence;
//      ** RANDOM, using C++ RANDOM function;
//      ** UNIFORM, using a simple uniform RNG;
//      ** USER, call the "user" routine.
//    * SAMPLE_NUM, the number of sampling points;
//    * BATCH, the number of sampling points to create at one time.
//    * LAT_IT_NUM, the maximum number of iterations;
//    * OUTPUT, a file in which to store the data.
//
//    To indicate that no further computations are desired, it is
//    enough to input a nonsensical value, such as -1.
//
//  Modified:
//
//    12 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    double CELL_GENERATOR(M,N), the Voronoi cell generators
//    of the Voronoi tessellation, as approximated by the CVT algorithm.  This
//    is the output quantity of most interest.
//
//    int CVT_IT, the number of iterations used in the Centroidal
//    Voronoi Tesselation calculation.  The default value is 10.
//
//    int LATIN_IT, the number of Latin hypercube iterations to carry out.
//    This defaults to 5.
//
//    int M, the spatial dimension.
//
//    int N, the number of Voronoi cells to generate.
//
//    int SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
//    -1, the sampling function is C++ RANDOM_NUMBER,
//    0, the sampling function is UNIFORM,
//    1, the sampling function is HALTON,
//    2, the sampling function is GRID.
//
//    int SAMPLE_FUNCTION_INIT, specifies how the initial
//    generators are chosen:
//    -1, the initialization function is C++ RANDOM,
//    0, the initialization function is UNIFORM,
//    1, the initialization function is HALTON,
//    2, the initialization function is GRID,
//    3, the initial values are read in from a file.
//
//    int SAMPLE_NUM, the number of sampling points used on
//    each CVT iteration.  A typical value is 5000 * N.
//
//    int SEED, determines how to initialize the random number routine.
//    If SEED is zero, then RANDOM_INITIALIZE will make up a seed
//    from the current real time clock reading.
//    If SEED is nonzero, then a reproducible sequence of random numbers
//    defined by SEED will be chosen.
//    By default, SEED initially has a value chosen by RANDOM_INITIALIZE,
//    but the user can reset SEED at any time.
//
{
  int batch;
  double cvt_energy;
  double cvt_it_diff;
  int cvt_it;
  int cvt_it_num ;
  bool debug = true;
  int dim_num;
  char input_file_name[80];
  int init;
  char init_string[80];
  double lat_energy;
  int lat_it;
  int lat_it_num;
  int n;
  int n_total;
  char output_file_name[80];
  double *r;
  bool reset;
  int sample;
  int sample_num;
  char sample_string[80];
  int seed;
  int seed_init;
//
//  Print introduction and options.
//
  timestamp ( );

  cout << "\n";
  cout << "LCVT_DATASET\n";
  cout << "  C++ version\n";
  cout << "  Create a \"Latinized\" CVT datasets.\n";
  cout << "\n";
  cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
  cout << "\n";
  cout << "  This program is meant to be used interactively.\n";
  cout << "  It is also possible to prepare a simple input \n";
  cout << "  file beforehand and use it in batch mode.\n";
  cout << "\n";
  cout << "  The program requests input values from the user:\n";
  cout << "\n";
  cout << "  * DIM_NUM, the spatial dimension,\n";
  cout << "  * N, the number of points to generate,\n";
  cout << "  * SEED_INIT, a seed to use for random number generation,\n";
  cout << "  * INIT, initialize the points:\n";
  cout << "    ** file, read data from a file;\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using C++ RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "    ** USER, call the \"user\" routine;\n";
  cout << "  * CVT_IT_NUM, the maximum number of CVT iterations.\n";
  cout << "  * SAMPLE, how to conduct the sampling.\n";
  cout << "    ** GRID, by picking points from a grid;\n";
  cout << "    ** HALTON, from a Halton sequence;\n";
  cout << "    ** RANDOM, using C++ RANDOM function;\n";
  cout << "    ** UNIFORM, using a simple uniform RNG;\n";
  cout << "    ** USER, call the \"user\" routine;\n";
  cout << "  * SAMPLE_NUM, the number of sample points.\n";
  cout << "  * BATCH, number of sample points to create at one time.\n";
  cout << "  * LAT_IT_NUM, the number of Latinizing iterations.\n";
  cout << "  * OUTPUT, a file in which to store the data.\n";
  cout << "\n";
  cout << "  To indicate that no further computations are\n";
  cout << "  desired, it is enough to input a nonsensical value,\n";
  cout << "  such as -1.\n";

  cout << "  *\n";
  cout << " *\n";
  cout << "*  Ready to generate a new dataset:\n";
  cout << " *\n";
  cout << "  *\n";

  cout << "  Enter DIM_NUM, the spatial dimension:\n";
  cout << "  (Try \"2\" if you do not have a preference.)\n";
  cout << "  (0 or any negative value terminates execution).\n";

  cin >> dim_num;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for DIM_NUM.\n";
    exit ( 1 );
  }

  cout << "  User input DIM_NUM = " << dim_num << "\n";

  if ( dim_num < 1 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of DIM_NUM = " << dim_num << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Enter N, the number of points to generate:\n";
  cout << "  (Try \"25\" if you do not have a preference.)\n";
  cout << "  (0 or any negative value terminates execution).\n";

  cin >> n;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for N.\n";
    exit ( 1 );
  }

  cout << "  User input N = " << n << "\n";

  if ( n < 1 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of N = " << n << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  Enter SEED_INIT, a seed for the random number generator:\n";
  cout << "  (Try \"123456789\" if you do not have a preference.)\n";
  cout << "  (Any negative value terminates execution).\n";

  cin >> seed_init;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for SEED_INIT.\n";
    exit ( 1 );
  }

  cout << "  User input SEED_INIT = " << seed_init << "\n";

  if ( seed_init < 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of SEED_INIT = " << seed_init << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  seed = seed_init;

  cout << "\n";
  cout << "  INIT is the method of initializing the data:\n";
  cout << "\n";
  cout << "  file     read data from a file;\n";
  cout << "  GRID     by picking points from a grid;\n";
  cout << "  HALTON   from a Halton sequence;\n";
  cout << "  RANDOM   using C++ RANDOM function;\n";
  cout << "  UNIFORM  using a simple uniform RNG;\n";
  cout << "  USER     call the \"user\" routine;\n";
  cout << "\n";
  cout << "  (Try \"RANDOM\" if you do not have a preference.)\n";
  cout << "  (A blank value terminates execution).\n";
  cout << "\n";
  cout << "  Enter INIT:\n";

  cin >> init_string;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for INIT\n";
    exit ( 1 );
  }

  if ( s_eqi ( init_string, "RANDOM"  ) )
  {
    init = -1;
    cout << "  User input INIT = \"RANDOM\"\n";
  }
  else if ( s_eqi ( init_string, "UNIFORM" ) )
  {
    init = 0;
    cout << "  User input INIT = \"UNIFORM\"\n";
  }
  else if ( s_eqi ( init_string, "HALTON"  ) )
  {
    init = 1;
    cout << "  User input INIT = \"HALTON\".\n";
  }
  else if ( s_eqi ( init_string, "GRID"    ) )
  {
    init = 2;
    cout << "  User input INIT = \"GRID\".\n";
  }
  else if ( s_eqi ( init_string, "USER"    ) )
  {
    init = 3;
    cout << "  User input INIT = \"USER\".\n";
  }
  else if ( 0 < s_len_trim ( init_string ) )
  {
    init = 4;
    cout << "  User input INIT = FILE_NAME = \"" << init_string << "\".\n";
    strcpy ( input_file_name, init_string );
  }
  else
  {
    cout << "n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of INIT \n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  CVT_IT_NUM is the number of CVT iterations.\n";
  cout << "\n";
  cout << "  A CVT iteration carries out the following steps:\n";
  cout << "  * the Voronoi region associated with each\n";
  cout << "    generator is estimated by sampling;\n";
  cout << "  * the centroid of each Voronoi region is estimated.\n";
  cout << "  * the generator is replaced by the centroid.\n";
  cout << "\n";
  cout << "  If \"enough\" sampling points are used,\n";
  cout << "  and \"enough\" iterations are taken, this process\n";
  cout << "  will converge\n";
  cout << "\n";
  cout << "  (Try \"50\" if you do not have a preference.)\n";
  cout << "  (A negative value terminates execution).\n";
  cout << "\n";
  cout << "  Enter CVT_IT_NUM:\n";

  cin >> cvt_it_num;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for CVT_IT_NUM.\n";
    exit ( 1 );
  }

  cout << "  User input CVT_IT_NUM = " << cvt_it_num << "\n";

  if ( cvt_it_num < 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of CVT_IT_NUM = " << cvt_it_num << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  SAMPLE is the method of sampling the region:\n";
  cout << "\n";
  cout << "  GRID     by picking points from a grid;\n";
  cout << "  HALTON   from a Halton sequence;\n";
  cout << "  RANDOM   using C++ RANDOM function;\n";
  cout << "  UNIFORM  using a simple uniform RNG;\n";
  cout << "  USER     call the \"user\" routine;\n";
  cout << "\n";
  cout << "  (Try \"RANDOM\" if you do not have a preference.)\n";
  cout << "  (A blank value terminates execution).\n";
  cout << "\n";
  cout << "  Enter SAMPLE:\n";

  cin >> sample_string;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for SAMPLE.\n";
    exit ( 1 );
  }

  if ( s_eqi ( sample_string, "RANDOM"  ) )
  {
    cout << "  User input INIT = \"RANDOM\".\n";
    sample = -1;
  }
  else if ( s_eqi ( sample_string, "UNIFORM" ) )
  {
    cout << "  User input INIT = \"UNIFORM\".\n";
    sample = 0;
  }
  else if ( s_eqi ( sample_string, "HALTON"  ) )
  {
    cout << "  User input INIT = \"HALTON\".\n";
    sample = 1;
  }
  else if ( s_eqi ( sample_string, "GRID"    ) )
  {
    cout << "  User input INIT = \"GRID\".\n";
    sample = 2;
  }
  else if ( s_eqi ( sample_string, "USER" ) )
  {
    cout << "  User input INIT = \"USER\".\n";
    sample = 3;
  }
  else
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of SAMPLE \n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  SAMPLE_NUM is the number of sample points for CVT.\n";
  cout << "\n";
  cout << "  The Voronoi regions will be explored by generating\n";
  cout << "  SAMPLE_NUM points.  For each sample point, the\n";
  cout << "  nearest generator is found.  Using more points\n";
  cout << "  gives a better estimate of these regions.\n";
  cout << "\n";
  cout << "  SAMPLE_NUM should be much larger than N, the\n";
  cout << "  number of generators.\n";
  cout << "\n";
  cout << "  (Try \"10000\" if you do not have a preference.)\n";
  cout << "  (A zero or negative value terminates execution.)\n";
  cout << "\n";
  cout << "  Enter SAMPLE_NUM:\n";

  cin >> sample_num;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for SAMPLE_NUM.\n";
    exit ( 1 );
  }

  cout << "  User input SAMPLE_NUM = " << sample_num << "\n";

  if ( sample_num <= 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of SAMPLE_NUM = " << sample_num << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  BATCH is the number of sample points to create\n";
  cout << "  at one time.\n";
  cout << "\n";
  cout << "  BATCH should be between 1 and SAMPLE_NUM.\n";
  cout << "\n";
  cout << "  It is FASTER to set BATCH to SAMPLE_NUM;\n";
  cout << "  setting BATCH to 1 requires the least memory.\n";
  cout << "\n";
  cout << "  (Try \"" << i4_min ( sample_num, 1000 ) <<
    "\" if you do not have a preference.)\n";
  cout << "  (A zero or negative value terminates execution.)\n";
  cout << "\n";
  cout << "  Enter BATCH:\n";

  cin >> batch;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for SAMPLE_NUM.\n";
    exit ( 1 );
  }

  cout << "  User input BATCH = " << batch << "\n";

  if ( batch <= 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of BATCH = " << batch << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  LAT_IT_NUM is the number of Latinizing iterations.\n";
  cout << "\n";
  cout << "  Each step of the latinizing iteration begins\n";
  cout << "  by carrying out CVT_IT_NUM steps of CVT iteration,\n";
  cout << "  after which the data is \"latinized\".\n";
  cout << "\n";
  cout << "  Often, one latinizing step is enough.\n";
  cout << "\n";
  cout << "  In some cases, it may be worth while to carry\n";
  cout << "  out several latinizing steps; that is, the\n";
  cout << "  Latinized data is smoothed by another series\n";
  cout << "  of CVT steps, then latinized, and so on.\n";
  cout << "\n";
  cout << "  (Try \"1\" if you do not have a preference.)\n";
  cout << "  (A negative value terminates execution).\n";
  cout << "\n";
  cout << "  Enter LAT_IT_NUM:\n";

  cin >> lat_it_num;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for LAT_IT_NUM.\n";
    exit ( 1 );
  }

  cout << "  User input LAT_IT_NUM = " << lat_it_num << "\n";

  if ( cvt_it_num < 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of LAT_IT_NUM = " << lat_it_num << "\n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }

  cout << "\n";
  cout << "  OUTPUT is the name of a file into which\n";
  cout << "  the computed data may be stored.\n";
  cout << "\n";
  cout << "  (Try \"lcvt.txt\" if you do not have a preference.)\n";
  cout << "  (A blank value terminates execution).\n";
  cout << "\n";
  cout << "  Enter OUTPUT:\n";

  cin >> output_file_name;

  if ( cin.rdstate ( ) )
  {
    cin.clear ( );
    cout << "\n";
    cout << "LCVT_DATASET - Warning!\n";
    cout << "  Terminating abnormally because of an I/O error\n";
    cout << "  while expecting input for OUTPUT.\n";
    exit ( 1 );
  }

  cout << "  User input OUTPUT = \"" << output_file_name << "\".\n";

  if ( s_len_trim ( output_file_name ) <= 0 )
  {
    cout << "\n";
    cout << "LCVT_DATASET\n";
    cout << "  The input value of OUTPUT \n";
    cout << "  is interpreted as a request for termination.\n";
    cout << "  Normal end of execution.\n";
    exit ( 1 );
  }
//
//  Initialize the data.
//
  if ( init == 4 )
  {
    r = r8table_data_read ( input_file_name, dim_num, n );
  }
  else
  {
    n_total = n;
    r = new double[dim_num*n];
    reset = true;

    region_sampler ( dim_num, n, n_total, r, init, reset, &seed );
  }

  if ( debug )
  {
    cout << "\n";
    cout << "  Latin IT      CVT Energy    Latin Energy\n";
    cout << "\n";
  }

  for ( lat_it = 1; lat_it <= lat_it_num; lat_it++ )
  {
    if ( debug )
    {
      cout << "\n";
      cout << "    CVT IT  Change\n";
      cout << "\n";
    }

    for ( cvt_it = 1; cvt_it <= cvt_it_num; cvt_it++ )
    {
      cvt_iteration ( dim_num, n, r, sample_num, sample, &seed, &cvt_it_diff );

      if ( debug )
      {
        cout << "  " << setw(8) << cvt_it
             << "  " << setw(14) << cvt_it_diff << "\n";
      }

    }

    if ( debug )
    {
      cout << "\n";
    }

    cvt_energy = cluster_energy ( dim_num, n, r, sample_num, sample, &seed );

    r8mat_latinize ( dim_num, n, r );

    lat_energy = cluster_energy ( dim_num, n, r, sample_num, sample, &seed );

        cout << "  " << setw(8) << lat_it
             << "  " << setw(14) << cvt_energy
             << "  " << setw(14) << lat_energy << "\n";
  }
//
//  Write the data to a file.
//
  lcvt_write ( dim_num, n, seed_init, init, input_file_name, sample,
    sample_num, cvt_it_num, cvt_energy, lat_it_num, lat_energy,
    r, output_file_name );

  delete [] r;

  cout << "\n";
  cout << "  The data was written to the file \""
    << output_file_name << "\".\n";


  cout << "\n";
  cout << "LCVT_DATASET:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 )
  {
    c = c - 32;
  }

  return c;
}
//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, char C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 )
  {
    c1 = c1 - 32;
  }
  if ( 97 <= c2 && c2 <= 122 )
  {
    c2 = c2 - 32;
  }

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************80

double cluster_energy ( int dim_num, int n, double cell_generator[],
  int sample_num, int sample_function_cvt, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    CLUSTER_ENERGY returns the energy of a dataset.
//
//  Discussion:
//
//    The energy is the integral of the square of the distance from each point
//    in the region to its nearest generator.
//
//  Modified:
//
//    08 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of generators.
//
//    Input, double CELL_GENERATOR[DIM_NUM*N], the coordinates of the points.
//
//    Input, int SAMPLE_NUM, the number of sample points to use.
//
//    Input, int SAMPLE_FUNCTION_CVT, specifies how the sampling is done.
//    -1, 'RANDOM', using C++ RANDOM function;
//     0, 'UNIFORM', using a simple uniform RNG;
//     1, 'HALTON', from a Halton sequence;
//     2, 'GRID', points from a grid;
//     3, 'USER', call "user" routine.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double CLUSTER_ENERGY, the estimated energy.
//
{
  double energy;
  int i;
  int j;
  int nearest;
  bool reset;
  double *x;

  x = new double [dim_num];

  energy = 0.0;
  reset = true;

  for ( j = 0; j < sample_num; j++ )
  {
//
//  Generate a sampling point X.
//
    region_sampler ( dim_num, 1, sample_num, x, sample_function_cvt,
      reset, seed );

    reset = false;
//
//  Find the nearest cell generator.
//
    nearest = find_closest ( dim_num, n, x, cell_generator );

    for ( i = 0; i < dim_num; i++ )
    {
      energy = energy
        + pow ( x[i] - cell_generator[i+nearest*dim_num], 2 );
    }
  }
//
//  Add the contribution to the energy.
//
  energy = energy / ( double ) ( sample_num );

  delete [] x;

  return energy;
}
//****************************************************************************80

void cvt_iteration ( int m, int n, double generator[], int sample_num,
 int sample_function_cvt, int *seed, double *change_l2 )

//****************************************************************************80
//
//  Purpose:
//
//    CVT_ITERATION takes one step of the CVT iteration.
//
//  Discussion:
//
//    The routine is given a set of points, called "generators", which
//    define a tessellation of the region into Voronoi cells.  Each point
//    defines a cell.  Each cell, in turn, has a centroid, but it is
//    unlikely that the centroid and the generator coincide.
//
//    Each time this CVT iteration is carried out, an attempt is made
//    to modify the generators in such a way that they are closer and
//    closer to being the centroids of the Voronoi cells they generate.
//
//    A large number of sample points are generated, and the nearest generator
//    is determined.  A count is kept of how many points were nearest to each
//    generator.  Once the sampling is completed, the location of all the
//    generators is adjusted.  This step should decrease the discrepancy
//    between the generators and the centroids.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of Voronoi cells.
//
//    Input/output, double GENERATOR[M*N], the Voronoi
//    cell generators.  On output, these have been modified
//
//    Input, int SAMPLE_NUM, the number of sample points.
//
//    Input, int SAMPLE_FUNCTION_CVT, region sampling function:
//    -1, sampling function is RANDOM (C++ STDLIB library function);
//    0, sampling function is UNIFORM;
//    1, sampling function is HALTON;
//    2, sampling function is GRID;
//
//    Input/output, int *SEED, the random number seed.
//
//    Output, double *CHANGE_L2, the L2 norm of the difference between
//    the input and output data.
//
{
  double *generator2;
  int *counter;
  int i;
  int j;
  int k;
  int nearest;
  bool reset;
  double *x;

  generator2 = new double[m*n];
  for ( k = 0; k < m*n; k++ )
  {
    generator2[k] = 0.0;
  }

  counter = new int[n];
  for ( i = 0; i < n; i++ )
  {
    counter[i] = 0;
  }

  x = new double[m];

  reset = true;
//
//  If we are using the C++ random number generator, then initialize using the current seed.
//  (Currently, if we are using RANDOM for both the initializing and sampling, we make this
//  call twice, which is inefficient and possibly misleading.)
//
  if ( sample_function_cvt == -1 )
  {
    srandom ( *seed );
  }

  for ( j = 0; j < sample_num; j++ )
  {
//
//  Generate a sampling point X.
//
    region_sampler ( m, 1, sample_num, x, sample_function_cvt, reset, seed );

    reset = false;
//
//  Find the nearest cell generator.
//
    nearest = find_closest ( m, n, x, generator );
//
//  Add X to the averaging data for GENERATOR(*,NEAREST).
//
    for ( i = 0; i < m; i++ )
    {
      generator2[nearest*m+i] = generator2[nearest*m+i] + x[i];
    }

    counter[nearest] = counter[nearest] + 1;
  }
//
//  Compute the new generators.
//
  for ( j = 0; j < n; j++ )
  {
    if ( counter[j] != 0 )
    {
      for ( i = 0; i < m; i++ )
      {
        generator2[j*m+i] = generator2[j*m+i] / ( ( double ) counter[j] );
      }
    }
  }
//
//  Determine the change.
//
  *change_l2 = 0.0;
  for ( k = 0; k < m*n; k++ )
  {
    *change_l2 = *change_l2 + pow ( ( generator2[k] - generator[k] ), 2 );
  }
  *change_l2 = sqrt ( *change_l2 );
//
//  Update.
//
  for ( k = 0; k < m*n; k++ )
  {
    generator[k] = generator2[k];
  }

  delete [] counter;
  delete [] generator2;
  delete [] x;

  return;
}
//****************************************************************************80

int find_closest ( int m, int n, double x[], double generator[] )

//****************************************************************************80
//
//  Purpose:
//
//    FIND_CLOSEST finds the Voronoi cell generator closest to a point X.
//
//  Discussion:
//
//    This routine finds the closest Voronoi cell generator by checking every
//    one.  For problems with many cells, this process can take the bulk
//    of the CPU time.  Other approaches, which group the cell generators into
//    bins, can run faster by a large factor.
//
//  Modified:
//
//    24 September 2006
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of cell generators.
//
//    Input, double X[M], the point to be checked.
//
//    Input, double GENERATOR[M*N], the cell generators.
//
//    Output, int FIND_CLOSEST, the index of the nearest cell generators.
//
{
  double dist_min;
  double dist;
  int i;
  int j;
  int nearest;

  nearest = 0;
  dist_min = 0.0;

  for ( j = 0; j < n; j++ )
  {
    dist = 0.0;
    for ( i = 0; i < m; i++ )
    {
      dist = dist + pow ( x[i] - generator[i+j*m], 2 );
    }

    if ( j == 0 || dist < dist_min )
    {
      dist_min = dist;
      nearest = j;
    }

  }

  return nearest;
}
//****************************************************************************80

int i4_max ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
{
  int value;

  if ( i2 < i1 )
  {
    value = i1;
  }
  else
  {
    value = i2;
  }
  return value;
}
//****************************************************************************80

int i4_min ( int i1, int i2 )

//****************************************************************************80
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Modified:
//
//    11 September 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1 and I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of i1 and i2.
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//****************************************************************************80

void i4_to_halton ( int seed, int base[], int ndim, double r[] )

//****************************************************************************80
//
//  Purpose:
//
//    I4_TO_HALTON computes an element of a Halton sequence.
//
//  Reference:
//
//    John Halton,
//    On the efficiency of certain quasi-random sequences of points
//    in evaluating multi-dimensional integrals,
//    Numerische Mathematik,
//    Volume 2, pages 84-90, 1960.
//
//  Modified:
//
//    28 February 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int SEED, the index of the desired element.
//    SEED = 0 is allowed, and returns R = 0.
//
//    Input, int BASE[NDIM], the Halton bases, which are usually
//    distinct prime numbers.  Each base must be greater than 1.
//
//    Input, int NDIM, the dimension of the elements of the sequence.
//
//    Output, double R[NDIM], the SEED-th element of the Halton sequence
//    for the given bases.
//
{
  double base_inv;
  int digit;
  int i;
  int seed2;

  for ( i = 0; i < ndim; i++ )
  {
    if ( base[i] <= 1 )
    {
      cout << "\n";
      cout << "I4_TO_HALTON - Fatal error!\n";
      cout << "  An input base is less than or equal to 1.\n";
      cout << "  BASE[" << i << "] = " << base[i] << "\n";
      exit ( 1 );
    }
  }

  for ( i = 0; i < ndim; i++ )
  {
    seed2 = seed;
    base_inv = 1.0 / ( ( double ) base[i] );
    r[i] = 0.0;

    while ( seed2 != 0 )
    {
      digit = seed2 % base[i];
      r[i] = r[i] + ( ( double ) digit ) * base_inv;
      base_inv = base_inv / ( ( double ) base[i] );
      seed2 = seed2 / base[i];
    }
  }

  return;
}
//****************************************************************************80

void lcvt_write ( int dim_num, int n, int seed_start, int sample_function_init,
  char* file_in_name, int sample_function_cvt, int sample_num, int cvt_it,
  double cvt_energy, int latin_it, double latin_energy, double cell_generator[],
  char *file_out_name )

//****************************************************************************80
//
//  Purpose:
//
//    LCVT_WRITE writes a Latinized CVT dataset to a file.
//
//  Discussion:
//
//    The initial lines of the file are comments, which begin with a
//    "#" character.
//
//    Thereafter, each line of the file contains the M-dimensional
//    components of the next entry of the dataset.
//
//  Modified:
//
//    09 September 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, int N, the number of points.
//
//    Input, int SEED_START, the initial random number seed.
//
//    Input, int SAMPLE_FUNCTION_INIT, specifies how the initial
//    generators are chosen:
//    -1, the initialization function is RANDOM (C++ intrinsic),
//    0, the initialization function is UNIFORM,
//    1, the initialization function is HALTON,
//    2, the initialization function is GRID,
//    3, the initial values are read in from a file.
//
//    Input, char *FILE_IN_NAME, the name of the file
//    from which initialization values were read for the generators,
//    if SAMPLE_FUNCTION_INIT = 3.
//
//    Input, int SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
//    -1, the sampling function is RANDOM (C++ intrinsic),
//    0, the sampling function is UNIFORM,
//    1, the sampling function is HALTON,
//    2, the sampling function is GRID.
//
//    Input, int SAMPLE_NUM, the number of sampling points used on
//    each CVT iteration.
//
//    Input, int CVT_IT, the number of CVT iterations.
//
//    Input, double CVT_ENERGY, the energy of the final CVT dataset.
//
//    Input, int LATIN_IT, the number of Latin iterations.
//
//    Input, double LATIN_ENERGY, the energy of the Latinized
//    CVT dataset.
//
//    Input, double CELL_GENERATOR[DIM_NUM*N], the points.
//
//    Input, char *FILE_OUT_NAME, the name of
//    the output file.
//
{
  bool comment = true;
  ofstream file_out;
  int i;
  int j;
  char *s;

  file_out.open ( file_out_name );

  if ( !file_out )
  {
    cout << "\n";
    cout << "LCVT_WRITE - Fatal error!\n";
    cout << "  Could not open the output file.\n";
    exit ( 1 );
  }

  s = timestring ( );

  if ( comment )
  {
    file_out << "#  " << file_out_name << "\n";
    file_out << "#  created by routine LCVT_WRITE.C" << "\n";
    file_out << "#  at " << s << "\n";
    file_out << "#\n";

    file_out << "#  Dimension DIM_NUM =        "  << dim_num       << "\n";
    file_out << "#  Number of points N =       "  << n             << "\n";
    file_out << "#  EPSILON (unit roundoff) =  "  << r8_epsilon ( ) << "\n";

    if ( sample_function_init == 0 ||
         sample_function_init == 1 ||
         sample_function_cvt == 0 ||
         sample_function_cvt == 0 )
    {
      file_out << "#\n";
      file_out << "#  Initial SEED      =        "  << seed_start     << "\n";
    }

    file_out << "#\n";
    if ( sample_function_init == -1 )
    {
      file_out << "#  Initialization by RANDOM (C++ STDLIB intrinsic).\n";
    }
    else if ( sample_function_init == 0 )
    {
      file_out << "#  Initialization by UNIFORM.\n";
    }
    else if ( sample_function_init == 1 )
    {
      file_out << "#  Initialization by HALTON.\n";
    }
    else if ( sample_function_init == 2 )
    {
      file_out << "#  Initialization by GRID.\n";
    }
    else if ( sample_function_init == 3 )
    {
      file_out << "#  Initialization from file \"" << file_in_name << "\".\n";
    }

    if ( sample_function_cvt == -1 )
    {
      file_out << "#  Sampling by RANDOM (C++ STDLIB intrinsic).\n";
    }
    else if ( sample_function_cvt == 0 )
    {
      file_out << "#  Sampling by UNIFORM.\n";
    }
    else if ( sample_function_cvt == 1 )
    {
      file_out << "#  Sampling by HALTON.\n";
    }
    else if ( sample_function_cvt == 2 )
    {
      file_out << "#  Sampling by GRID.\n";
    }

    file_out << "#  Number of sample points    "  << sample_num << "\n";
    file_out << "#  Number of CVT iterations   "  << cvt_it << "\n";
    file_out << "#  Energy of CVT dataset      "  << cvt_energy << "\n";
    file_out << "#  Number of Latin iterations "  << latin_it << "\n";
    file_out << "#  Energy of Latin dataset    "  << latin_energy << "\n";

    file_out << "#\n";
  }

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < dim_num; i++ )
    {
      file_out << setw(10) << cell_generator[i+j*dim_num] << "  ";
    }
    file_out << "\n";
  }

  file_out.close ( );

  delete [] s;

  return;
}
//****************************************************************************80

int prime ( int n )

//****************************************************************************80
//
//  Purpose:
//
//    PRIME returns any of the first PRIME_MAX prime numbers.
//
//  Discussion:
//
//    PRIME_MAX is 1600, and the largest prime stored is 13499.
//
//    Thanks to Bart Vandewoestyne for pointing out a typo, 18 February 2005.
//
//  Modified:
//
//    18 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Milton Abramowitz and Irene Stegun,
//    Handbook of Mathematical Functions,
//    US Department of Commerce, 1964, pages 870-873.
//
//    Daniel Zwillinger,
//    CRC Standard Mathematical Tables and Formulae,
//    30th Edition,
//    CRC Press, 1996, pages 95-98.
//
//  Parameters:
//
//    Input, int N, the index of the desired prime number.
//    In general, is should be true that 0 <= N <= PRIME_MAX.
//    N = -1 returns PRIME_MAX, the index of the largest prime available.
//    N = 0 is legal, returning PRIME = 1.
//
//    Output, int PRIME, the N-th prime.  If N is out of range, PRIME
//    is returned as -1.
//
{
# define PRIME_MAX 1600

  int npvec[PRIME_MAX] = {
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409,
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463,
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541,
      547,  557,  563,  569,  571,  577,  587,  593,  599,  601,
      607,  613,  617,  619,  631,  641,  643,  647,  653,  659,
      661,  673,  677,  683,  691,  701,  709,  719,  727,  733,
      739,  743,  751,  757,  761,  769,  773,  787,  797,  809,
      811,  821,  823,  827,  829,  839,  853,  857,  859,  863,
      877,  881,  883,  887,  907,  911,  919,  929,  937,  941,
      947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013,
     1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
     1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
     1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223,
     1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291,
     1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373,
     1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451,
     1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511,
     1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583,
     1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657,
     1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733,
     1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811,
     1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889,
     1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987,
     1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053,
     2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129,
     2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213,
     2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287,
     2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357,
     2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423,
     2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531,
     2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617,
     2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687,
     2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741,
     2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819,
     2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903,
     2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999,
     3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079,
     3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181,
     3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257,
     3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331,
     3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413,
     3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511,
     3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571,
     3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643,
     3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727,
     3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821,
     3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907,
     3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989,
     4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057,
     4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139,
     4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231,
     4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297,
     4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409,
     4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493,
     4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583,
     4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657,
     4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751,
     4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831,
     4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937,
     4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003,
     5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087,
     5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179,
     5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279,
     5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387,
     5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443,
     5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521,
     5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639,
     5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693,
     5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791,
     5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857,
     5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939,
     5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053,
     6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133,
     6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221,
     6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301,
     6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367,
     6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473,
     6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571,
     6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673,
     6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761,
     6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833,
     6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917,
     6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997,
     7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103,
     7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207,
     7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297,
     7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411,
     7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499,
     7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561,
     7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643,
     7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723,
     7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829,
     7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919,
     7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017,
     8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111,
     8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219,
     8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291,
     8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387,
     8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501,
     8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597,
     8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677,
     8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741,
     8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831,
     8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929,
     8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011,
     9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109,
     9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199,
     9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283,
     9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377,
     9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439,
     9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533,
     9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631,
     9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733,
     9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811,
     9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887,
     9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007,
    10009,10037,10039,10061,10067,10069,10079,10091,10093,10099,
    10103,10111,10133,10139,10141,10151,10159,10163,10169,10177,
    10181,10193,10211,10223,10243,10247,10253,10259,10267,10271,
    10273,10289,10301,10303,10313,10321,10331,10333,10337,10343,
    10357,10369,10391,10399,10427,10429,10433,10453,10457,10459,
    10463,10477,10487,10499,10501,10513,10529,10531,10559,10567,
    10589,10597,10601,10607,10613,10627,10631,10639,10651,10657,
    10663,10667,10687,10691,10709,10711,10723,10729,10733,10739,
    10753,10771,10781,10789,10799,10831,10837,10847,10853,10859,
    10861,10867,10883,10889,10891,10903,10909,10937,10939,10949,
    10957,10973,10979,10987,10993,11003,11027,11047,11057,11059,
    11069,11071,11083,11087,11093,11113,11117,11119,11131,11149,
    11159,11161,11171,11173,11177,11197,11213,11239,11243,11251,
    11257,11261,11273,11279,11287,11299,11311,11317,11321,11329,
    11351,11353,11369,11383,11393,11399,11411,11423,11437,11443,
    11447,11467,11471,11483,11489,11491,11497,11503,11519,11527,
    11549,11551,11579,11587,11593,11597,11617,11621,11633,11657,
    11677,11681,11689,11699,11701,11717,11719,11731,11743,11777,
    11779,11783,11789,11801,11807,11813,11821,11827,11831,11833,
    11839,11863,11867,11887,11897,11903,11909,11923,11927,11933,
    11939,11941,11953,11959,11969,11971,11981,11987,12007,12011,
    12037,12041,12043,12049,12071,12073,12097,12101,12107,12109,
    12113,12119,12143,12149,12157,12161,12163,12197,12203,12211,
    12227,12239,12241,12251,12253,12263,12269,12277,12281,12289,
    12301,12323,12329,12343,12347,12373,12377,12379,12391,12401,
    12409,12413,12421,12433,12437,12451,12457,12473,12479,12487,
    12491,12497,12503,12511,12517,12527,12539,12541,12547,12553,
    12569,12577,12583,12589,12601,12611,12613,12619,12637,12641,
    12647,12653,12659,12671,12689,12697,12703,12713,12721,12739,
    12743,12757,12763,12781,12791,12799,12809,12821,12823,12829,
    12841,12853,12889,12893,12899,12907,12911,12917,12919,12923,
    12941,12953,12959,12967,12973,12979,12983,13001,13003,13007,
    13009,13033,13037,13043,13049,13063,13093,13099,13103,13109,
    13121,13127,13147,13151,13159,13163,13171,13177,13183,13187,
    13217,13219,13229,13241,13249,13259,13267,13291,13297,13309,
    13313,13327,13331,13337,13339,13367,13381,13397,13399,13411,
    13417,13421,13441,13451,13457,13463,13469,13477,13487,13499 };

  if ( n == -1 )
  {
    return PRIME_MAX;
  }
  else if ( n == 0 )
  {
    return 1;
  }
  else if ( n <= PRIME_MAX )
  {
    return npvec[n-1];
  }
  else
  {
    cout << "\n";
    cout << "PRIME - Fatal error!\n";
    cout << "  Unexpected input value of n = " << n << "\n";
    exit ( 1 );
  }

  return 0;
# undef PRIME_MAX
}
//****************************************************************************80

double r8_epsilon ( )

//****************************************************************************80
//
//  Purpose:
//
//    R8_EPSILON returns the R8 roundoff unit.
//
//  Discussion:
//
//    The roundoff unit is a number R which is a power of 2 with the
//    property that, to the precision of the computer's arithmetic,
//      1 < 1 + R
//    but
//      1 = ( 1 + R / 2 )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 September 2012
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Output, double R8_EPSILON, the R8 round-off unit.
//
{
  const double value = 2.220446049250313E-016;

  return value;
}
//****************************************************************************80

double r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Modified:
//
//    11 August 2004
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

void r8mat_latinize ( int m, int n, double table[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8MAT_LATINIZE "Latinizes" an R8MAT.
//
//  Discussion:
//
//    It is assumed, though not necessary, that the input dataset
//    has points that lie in the unit hypercube.
//
//    In any case, the output dataset will have this property.
//
//  Modified:
//
//    06 December 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of cells.
//
//    Input/output, double TABLE[M*N].  On input, the dataset to
//    be "Latinized".  On output, the Latinized dataset.
//
{
  double *column;
  int i;
  int *indx;
  int j;

  column = new double[n];

  for ( i = 0; i < m; i++ )
  {
    for ( j = 0; j < n; j++ )
    {
      column[j] = table[i+j*m];
    }
    indx = r8vec_sort_heap_index_a ( n, column );

    for ( j = 0; j < n; j++ )
    {
      table[i+indx[j]*m] = ( double ) ( 2 * j + 1 ) / ( double ) ( 2 * n );
    }

    delete [] indx;
  }

  delete [] column;

  return;
}
//****************************************************************************80

double *r8table_data_read ( char *input_filename, int m, int n )

//****************************************************************************80
//
//  Purpose:
//
//    R8TABLE_DATA_READ reads the data from an R8TABLE file.
//
//  Discussion:
//
//    The file is assumed to contain one record per line.
//
//    Records beginning with the '#' character are comments, and are ignored.
//    Blank lines are also ignored.
//
//    Each line that is not ignored is assumed to contain exactly (or at least)
//    M real numbers, representing the coordinates of a point.
//
//    There are assumed to be exactly (or at least) N such records.
//
//  Modified:
//
//    27 January 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *INPUT_FILENAME, the name of the input file.
//
//    Input, int M, the number of spatial dimensions.
//
//    Input, int N, the number of points.  The program
//    will stop reading data once N values have been read.
//
//    Output, double DTABLE_DATA_READ[M*N], the table data.
//
{
  bool error;
  ifstream input;
  int i;
  int j;
  char line[255];
  double *table;
  double *x;

  input.open ( input_filename );

  if ( !input )
  {
    cout << "\n";
    cout << "R8TABLE_DATA_READ - Fatal error!\n";
    cout << "  Could not open the input file: \"" << input_filename << "\"\n";
    return NULL;
  }

  table = new double[m*n];

  x = new double[m];

  j = 0;

  while ( j < n )
  {
    input.getline ( line, sizeof ( line ) );

    if ( input.eof ( ) )
    {
      break;
    }

    if ( line[0] == '#' || s_len_trim ( line ) == 0 )
    {
      continue;
    }

    error = s_to_r8vec ( line, m, x );

    if ( error )
    {
      continue;
    }

    for ( i = 0; i < m; i++ )
    {
      table[i+j*m] = x[i];
    }
    j = j + 1;

  }

  input.close ( );

  delete [] x;

  return table;
}
//****************************************************************************80

int *r8vec_sort_heap_index_a ( int n, double a[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_SORT_HEAP_INDEX_A does an indexed heap ascending sort of an R8VEC.
//
//  Discussion:
//
//    The sorting is not actually carried out.  Rather an index array is
//    created which defines the sorting.  This array may be used to sort
//    or index the array, or to sort or index related arrays keyed on the
//    original array.
//
//    Once the index array is computed, the sorting can be carried out
//    "implicitly:
//
//      A(INDX(I)), I = 1 to N is sorted,
//
//    after which A(I), I = 1 to N is sorted.
//
//  Modified:
//
//    30 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the array.
//
//    Input, double A[N], an array to be index-sorted.
//
//    Output, int R8VEC_SORT_HEAP_INDEX_A[N], contains the sort index.  The
//    I-th element of the sorted array is A(INDX(I)).
//
{
  double aval;
  int i;
  int *indx;
  int indxt;
  int ir;
  int j;
  int l;

  indx = new int[n];

  for ( i = 1; i <= n; i++ )
  {
    indx[i-1] = i;
  }

  l = n / 2 + 1;
  ir = n;

  for ( ; ; )
  {
    if ( 1 < l )
    {
      l = l - 1;
      indxt = indx[l-1];
      aval = a[indxt-1];
    }
    else
    {
      indxt = indx[ir-1];
      aval = a[indxt-1];
      indx[ir-1] = indx[0];
      ir = ir - 1;

      if ( ir == 1 )
      {
        indx[0] = indxt;
        for ( i = 0; i < n; i++ )
        {
          indx[i] = indx[i] - 1;
        }
        return indx;
      }

    }

    i = l;
    j = l + l;

    while ( j <= ir )
    {
      if ( j < ir )
      {
        if ( a[indx[j-1]-1] < a[indx[j]-1] )
        {
          j = j + 1;
        }
      }

      if ( aval < a[indx[j-1]-1] )
      {
        indx[i-1] = indx[j-1];
        i = j;
        j = j + j;
      }
      else
      {
        j = ir + 1;
      }
    }
    indx[i-1] = indxt;
  }
}
//****************************************************************************80

void region_sampler ( int m, int n, int n_total, double x[],
  int sample_function, bool reset, int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    REGION_SAMPLER returns a sample point in the physical region.
//
//  Discussion:
//
//    This routine original interfaced with a lower routine called
//    TEST_REGION, which tested whether the points generated in the
//    bounding box were actually inside a possibly smaller physical
//    region of interest.  It's been a long time since that option
//    was actually used, so it's been dropped.
//
//    A point is chosen in the bounding box, either by a uniform random
//    number generator, or from a vector Halton sequence.
//
//    The entries of the local vector HALTON_BASE should be distinct primes.
//    Right now, we're assuming M is no greater than 3.
//
//  Modified:
//
//    04 March 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of points to generate now.
//
//    Input, int N_TOTAL, the total number of points to generate.
//
//    Output, double X[M*N], the sample points.
//
//    Input, int SAMPLE_FUNCTION, region sampling function:
//    -1, sampling function is RANDOM (C++ STDLIB library function);
//    0, sampling function is UNIFORM;
//    1, sampling function is HALTON;
//    2, sampling function is GRID;
//    3, sample points are generated elsewhere, and this routine is skipped.
//
//    Input, bool RESET, if true, then this is the first call for a particular
//    calculation, and initialization should be taken care of.
//
//    Input/output, int *SEED, the random number seed.
//
{
  double exponent;
  static int *halton_base = NULL;
  static int halton_seed = 1;
  int i;
  int ival;
  int j;
  int k;
  static int ngrid;
  static int rank;
  static int *tuple = NULL;

  if ( sample_function == -1 )
  {
    for ( k = 0; k < m*n; k++ )
    {
      x[k] = ( double ) random ( ) / ( double ) RAND_MAX;
    }
  }
  else if ( sample_function == 0 )
  {
    for ( k = 0; k < m*n; k++ )
    {
      x[k] = r8_uniform_01 ( seed );
    }
  }
  else if ( sample_function == 1 )
  {
    if ( reset )
    {
      halton_seed = 1;
      reset = false;
      if ( halton_base )
      {
        delete [] halton_base;
      }
      halton_base = new int[m];
      for ( i = 0; i < m; i++ )
      {
        halton_base[i] = prime ( i+1 );
      }
    }
//
//  The unusual syntax X+J*M essentially means pass the address of the beginning
//  of the J-th vector of length M in X.
//
    for ( j = 0; j < n; j++ )
    {
      i4_to_halton ( halton_seed, halton_base, m, x+j*m );
      halton_seed = halton_seed + 1;
    }
  }
  else if ( sample_function == 2 )
  {
    if ( reset )
    {
      rank = 0;
      exponent = 1.0 / ( ( double ) ( m ) );

      ngrid = ( int ) pow ( ( double ) n_total, exponent );

      if ( pow ( ( double ) ngrid, m ) < n_total )
      {
        ngrid = ngrid + 1;
      }
      if ( tuple != NULL )
      {
        delete [] tuple;
      }
      tuple = new int[m];
      reset = false;
    }

    for ( j = 0; j < n; j++ )
    {
      tuple_next_fast ( ngrid, m, rank, tuple );
      rank = rank + 1;
      for ( i = 0; i < m; i++ )
      {
        x[j*m+i] = ( ( double ) ( 2 * tuple[i] - 1 ) )
                 / ( ( double ) ( 2 * ngrid ) );
      }
    }
  }
  else if ( sample_function == 3 )
  {
  }
  else
  {
    cout << "\n";
    cout << "REGION_SAMPLER - Fatal error!\n";
    cout << "  Illegal SAMPLE_FUNCTION = " << sample_function << "\n";
    exit ( 1 );
  }

  return;
}
//****************************************************************************80

bool s_eqi ( char *s1, char *s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Modified:
//
//    05 May 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S1, char *S2, pointers to two strings.
//
//    Output, bool S_EQI, is true if the strings are equal.
//
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  nchar = i4_min ( nchar1, nchar2 );
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < nchar1 )
  {
    for ( i = nchar; i < nchar1; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return false;
      }
    }
  }
  else if ( nchar < nchar2 )
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return false;
      }
    }
  }

  return true;

}
//****************************************************************************80

int s_len_trim ( char* s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Modified:
//
//    26 April 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, a pointer to a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;
  char* t;

  n = strlen ( s );
  t = s + strlen ( s ) - 1;

  while ( 0 < n )
  {
    if ( *t != ' ' )
    {
      return n;
    }
    t--;
    n--;
  }

  return n;
}
//****************************************************************************80

double s_to_r8 ( char *s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Examples:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Modified:
//
//    07 August 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_D, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80

bool s_to_r8vec ( char *s, int n, double rvec[] )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8VEC reads an R8VEC from a string.
//
//  Modified:
//
//    19 February 2001
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char *S, the string to be read.
//
//    Input, int N, the number of values expected.
//
//    Output, double RVEC[N], the values read from the string.
//
//    Output, bool S_TO_R8VEC, is true if an error occurred.
//
{
  bool error;
  int i;
  int lchar;
  double x;

  for ( i = 0; i < n; i++ )
  {
    rvec[i] = s_to_r8 ( s, &lchar, &error );

    if ( error )
    {
      return error;
    }

    s = s + lchar;

  }

  return error;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

void tuple_next_fast ( int m, int n, int rank, int x[] )

//****************************************************************************80
//
//  Purpose:
//
//    TUPLE_NEXT_FAST computes the next element of a tuple space, "fast".
//
//  Discussion:
//
//    The elements are N vectors.  Each entry is constrained to lie
//    between 1 and M.  The elements are produced one at a time.
//    The first element is
//      (1,1,...,1)
//    and the last element is
//      (M,M,...,M)
//    Intermediate elements are produced in lexicographic order.
//
//  Example:
//
//    N = 2,
//    M = 3
//
//    INPUT        OUTPUT
//    -------      -------
//    Rank  X      X
//    ----  ---    ---
//    0     * *    1 1
//    1     1 1    1 2
//    2     1 2    1 3
//    3     1 3    2 1
//    4     2 1    2 2
//    5     2 2    2 3
//    6     2 3    3 1
//    7     3 1    3 2
//    8     3 2    3 3
//    9     3 3    1 1
//
//  Modified:
//
//    28 April 2003
//
//  Parameters:
//
//    Input, int M, the maximum entry.
//
//    Input, int N, the number of components.
//
//    Input, int RANK, indicates the rank of the tuples.
//    On the very first call only, it is necessary that
//    the user set RANK = 0.
//
//    Input/output, int X[N], on input the previous tuple.
//    On output, the next tuple.
//
{
  static int *base = NULL;
  int i;

  if ( rank == 0 )
  {
    if ( base )
    {
      delete [] base;
    }
    base = new int[n];

    base[n-1] = 1;
    for ( i = n-2; 0 <= i; i-- )
    {
      base[i] = base[i+1] * m;
    }
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = ( ( rank / base[i] ) % m ) + 1;
  }

  return;
}
