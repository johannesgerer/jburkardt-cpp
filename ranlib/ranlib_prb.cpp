# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <cstring>

using namespace std;

# include "ranlib.hpp"
# include "rnglib.hpp"

int main ( );
void test_phrtsd ( char *phrase );
void test_bot ( );
void test_genbet ( char *phrase );
void test_ignbin ( char *phrase );
void test_genchi ( char *phrase );
void test_genexp ( char *phrase );
void test_genf ( char *phrase );
void test_gengam ( char *phrase );
void test_ignnbn ( char *phrase );
void test_gennch ( char *phrase );
void test_gennf ( char *phrase );
void test_gennor ( char *phrase );
void test_ignpoi ( char *phrase );
void test_genunf ( char *phrase );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RANLIB_PRB.
//
//  Discussion:
//
//    RANLIB_PRB tests the RANLIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  char phrase[] = "randomizer";

  timestamp ( );
  cout << "\n";
  cout << "RANLIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the RANLIB library.\n";

  test_phrtsd ( phrase );

  test_bot ( );

  test_genbet ( phrase );
  test_ignbin ( phrase );
  test_genchi ( phrase );
  test_genexp ( phrase );
  test_genf ( phrase );
  test_gengam ( phrase );
  test_ignnbn ( phrase );
  test_gennch ( phrase );
  test_gennf ( phrase );
  test_gennor ( phrase );
  test_ignpoi ( phrase );
  test_genunf ( phrase );
//
//  Terminate.
//
  cout << "\n";
  cout << "RANLIB_PRB\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test_phrtsd ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_PHRTSD tests PHRTSD, which generates two seeds from a phrase.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  int seed1;
  int seed2;

  cout << "\n";
  cout << "TEST_PHRTSD\n";
  cout << "  Test PHRTST,\n";
  cout << "  which generates two seeds from a phrase.\n";

  cout << "\n";
  cout << "  Randomizing phrase is \"" << phrase << "\"\n";

  phrtsd ( phrase, seed1, seed2 );

  cout << "\n";
  cout << "  Seed1 = " << seed1 << "\n";
  cout << "  Seed2 = " << seed2 << "\n";

  return;
}
//****************************************************************************80

void test_bot ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_BOT is a test program for the bottom level routines
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  int answer[10000];
  int genlst[5] = { 0, 4, 9, 19, 31 };
  int ians;
  int iblock;
  int igen;
  int itmp;
  int ix;
  int ixgen;
  int nbad;
  int seed1;
  int seed2;

  cout << "\n";
  cout << "TEST_BOT\n";
  cout << "  Test the lower level random number generators.\n";
  cout << "\n";
  cout << "  Five of the 32 generators will be tested.\n";
  cout << "  We generate 100000 numbers, reset the block\n";
  cout << "  and do it again.  No disagreements should occur.\n";
  cout << "\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set up all generators.
//
  seed1 = 12345;
  seed2 = 54321;
  set_initial_seed ( seed1, seed2 );
//
//  For a selected set of generators
//
  nbad = 0;

  for ( ixgen = 0; ixgen < 5; ixgen++ )
  {
    igen = genlst[ixgen];
    cgn_set ( igen );
    cout << "  Testing generator " << igen << "\n";
//
//  Use 10 blocks, and generate 1000 numbers per block
//
    init_generator ( 0 );

    for ( iblock = 0; iblock < 10; iblock++ )
    {
      for ( ians = 0; ians < 1000; ians++ )
      {
        ix = ians + iblock * 1000;
        answer[ix] = i4_uni ( );
      }
      init_generator ( 2 );
    }
//
//  Do it again and compare answers
//  Use 10 blocks, and generate 1000 numbers.
//
    init_generator ( 0 );

    for ( iblock = 0; iblock < 10; iblock++ )
    {
      for ( ians = 0; ians < 1000; ians++ )
      {
        ix = ians + iblock * 1000;
        itmp = i4_uni ( );

        if ( itmp != answer[ix] )
        {
          cout << "\n";
          cout << "TEST_BOT - Warning!\n";
          cout << "  Data disagreement:\n";
          cout << "  Block = " << iblock << "\n";
          cout << "  N within block = " << ians << "\n";
          cout << "  Index in ANSWER = " << ix << "\n";
          cout << "  First value =  " << answer[ix] << "\n";
          cout << "  Second value = " << itmp << "\n";

          nbad = nbad + 1;

          if ( 10 < nbad )
          {
            cout << "\n";
            cout << "TEST_BOT - Warning!\n";
            cout << "  More than 10 mismatches!\n";
            cout << "  Tests terminated early.\n";
            return;
          }
        }
      }
      init_generator ( 2 );
    }
  }
  return;
}
//****************************************************************************80

void test_genbet ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENBET tests GENBET, which generates Beta deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float *array;
  float av;
  float avtr;
  float b;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  string pdf = "bet";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENBET\n";
  cout << "  Test GENBET,\n";
  cout << "  which generates Beta deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = 1.0;
  high = 10.0;
  b = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = genbet ( a, b );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = a;
  param[1] = b;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         "
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_ignbin ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_IGNBIN tests IGNBIN, which generates Binomial deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 10000;
  int nn;
  float param[2];
  char pdf[] = "bin";
  float pp;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_IGNBIN\n";
  cout << "  Test IGNBIN,\n";
  cout << "  which generates binomial deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 0.5;
  high = 20.0;
  nn = ( int ) genunf ( low, high );

  low = 0.0;
  high = 1.0;
  pp = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  NN = " << nn << "\n";
  cout << "  PP = " << pp << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignbin ( nn, pp );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = ( float ) ( nn );
  param[1] = pp;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_genchi ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENCHI tests GENCHI, which generates Chi-Square deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float df;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[1];
  char pdf[] = "chi";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENCHI\n";
  cout << "  Test GENCHI,\n";
  cout << "  which generates Chi-square deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 1.0;
  high = 10.0;
  df = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  DF = " << df << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = genchi ( df );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = df;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_genexp ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENEXP tests GENEXP, which generates exponential deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[2];
  char pdf[] = "exp";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENEXP\n";
  cout << "  Test GENEXP,\n";
  cout << "  which generates exponential deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low =  0.5;
  high = 10.0;
  mu = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  MU = " << mu << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = genexp ( mu );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = mu;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_genf ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENF tests GENF, which generates F deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float dfd;
  float dfn;
  float high;
  int i;
  float low;
  int n = 10000;
  float param[2];
  char pdf[] = "f";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENF\n";
  cout << "  Test GENF,\n";
  cout << "  which generates F deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 3.0;
  high = 10.0;
  dfn = genunf ( low, high );

  low = 5.0;
  high = 10.0;
  dfd = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  DFN =   " << dfn << "\n";
  cout << "  DFD =   " << dfd << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = genf ( dfn, dfd );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = dfn;
  param[1] = dfd;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_gengam ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENGAM tests GENGAM, which generates Gamma deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "gam";
  float r;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENGAM\n";
  cout << "  Test GENGAM,\n";
  cout << "  which generates Gamma deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = 1.0;
  high = 10.0;
  r = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  R = " << r << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = gengam ( a, r );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = a;
  param[1] = r;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_ignnbn ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_IGNNBN tests IGNNBN, which generates Negative Binomial deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  int n = 10000;
  int nn;
  float param[2];
  char pdf[] = "nbn";
  float pp;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_IGNNBN\n";
  cout << "  Test IGNNBN,\n";
  cout << "  which generates negative binomial deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 3.0;
  high = 20.0;
  nn = ( int ) genunf ( low, high );

  low = 0.0;
  high = 1.0;
  pp = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  NN = " << nn << "\n";
  cout << "  PP = " << pp << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignnbn ( nn, pp );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = ( float ) ( nn );
  param[1] = pp;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_gennch ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENNCH tests GENNCH, which generates noncentral Chi-Square deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float df;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "nch";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;
  float xnonc;

  cout << "\n";
  cout << "TEST_GENNCH\n";
  cout << "  Test GENNCH,\n";
  cout << "  which generates noncentral Chi-square deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 2.0;
  high = 10.0;
  df = genunf ( low, high );

  low = 0.0;
  high = 2.0;
  xnonc = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  DF =    " << df << "\n";
  cout << "  XNONC = " << xnonc << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennch ( df, xnonc );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = df;
  param[1] = xnonc;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_gennf ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENNF tests GENNF, which generates noncentral F deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float dfd;
  float dfn;
  float high;
  int i;
  float low;
  int n = 10000;
  float param[3];
  char pdf[] = "nf";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;
  float xnonc;

  cout << "\n";
  cout << "TEST_GENNF\n";
  cout << "  Test GENNF,\n";
  cout << "  which generates noncentral F deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 3.0;
  high = 10.0;
  dfn = genunf ( low, high );

  low = 5.0;
  high = 10.0;
  dfd = genunf ( low, high );

  low = 0.0;
  high = 2.0;
  xnonc = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  DFN =   " << dfn << "\n";
  cout << "  DFD =   " << dfd << "\n";
  cout << "  XNONC = " << xnonc << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennf ( dfn, dfd, xnonc );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = dfn;
  param[1] = dfd;
  param[2] = xnonc;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_gennor ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENNOR tests GENNOR, which generates normal deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[2];
  char pdf[] = "nor";
  float sd;
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENNOR\n";
  cout << "  Test GENNOR,\n";
  cout << "  which generates normal deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = -10.0;
  high = 10.0;
  mu = genunf ( low, high );

  low = 0.25;
  high = 4.0;
  sd = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << " '  MU =   " << mu << "\n";
  cout << " '  SD =   " << sd << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = gennor ( mu, sd );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = mu;
  param[1] = sd;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_ignpoi ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_IGNPOI tests IGNPOI, which generates Poisson deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float *array;
  float av;
  float avtr;
  float high;
  int i;
  float low;
  float mu;
  int n = 1000;
  float param[1];
  char pdf[] = "poi";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_IGNPOI\n";
  cout << "  Test IGNPOI,\n";
  cout << "  which generates Poisson deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 0.5;
  high = 20.0;
  mu = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  MU = " << mu << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = ( float ) ignpoi ( mu );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = mu;

  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         " 
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
//****************************************************************************80

void test_genunf ( char *phrase )

//****************************************************************************80
//
//  Purpose:
//
//    TEST_GENUNF tests GENUNF, which generates uniform deviates.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 April 2013
//
//  Author:
//
//    John Burkardt
//
{
  float a;
  float *array;
  float av;
  float avtr;
  float b;
  float high;
  int i;
  float low;
  int n = 1000;
  float param[2];
  char pdf[] = "unf";
  int seed1;
  int seed2;
  float var;
  float vartr;
  float xmax;
  float xmin;

  cout << "\n";
  cout << "TEST_GENUNF\n";
  cout << "  Test GENUNF,\n";
  cout << "  which generates uniform deviates.\n";
//
//  Initialize the generators.
//
  initialize ( );
//
//  Set the seeds based on the phrase.
//
  phrtsd ( phrase, seed1, seed2 );
//
//  Initialize all generators.
//
  set_initial_seed ( seed1, seed2 );
//
//  Select the parameters at random within a given range.
//
  low = 1.0;
  high = 10.0;
  a = genunf ( low, high );

  low = a + 1.0;
  high = a + 10.0;
  b = genunf ( low, high );

  cout << "\n";
  cout << "  N = " << n << "\n";
  cout << "\n";
  cout << "  Parameters:\n";
  cout << "\n";
  cout << "  A = " << a << "\n";
  cout << "  B = " << b << "\n";
//
//  Generate N samples.
//
  array = new float[n];
  for ( i = 0; i < n; i++ )
  {
    array[i] = genunf ( a, b );
  }
//
//  Compute statistics on the samples.
//
  stats ( array, n, av, var, xmin, xmax );
//
//  Request expected value of statistics for this distribution.
//
  param[0] = a;
  param[1] = b;
  trstat ( pdf, param, avtr, vartr );

  cout << "\n";
  cout << "  Sample data range:         "
       << "  " << setw(14) << xmin
       << "  " << setw(14) << xmax << "\n";
  cout << "  Sample mean, variance:     "
       << "  " << setw(14) << av
       << "  " << setw(14) << var << "\n";
  cout << "  Distribution mean, variance"
       << "  " << setw(14) << avtr
       << "  " << setw(14) << vartr << "\n";

  delete [] array;

  return;
}
