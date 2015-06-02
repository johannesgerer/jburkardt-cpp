# include <cstdlib>
# include <iostream>
# include <iomanip>

using namespace std;

# include "rnglib.hpp"

int main ( );
void test01 ( );
void test02 ( );
void test03 ( );
void test04 ( );

//****************************************************************************80

int main ( )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for RNGLIB_PRB.
//
//  Discussion:
//
//    RNGLIB_PRB tests the RNGLIB library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    28 March 2013
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );
  cout << "\n";
  cout << "RNGLIB_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the RNGLIB library.\n";
//
//  Call tests.
//
  test01 ( );
  test02 ( );
  test03 ( );
  test04 ( );
//
//  Terminate.
//
  cout << "\n";
  cout << "RNGLIB_PRB\n";
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
//    TEST01 calls I4_UNI 10 times, just to show how it is done.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int g;
  int i;
  int j;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4_UNI ( ) returns a random positive integer\n";
  cout << "  using the current generator.\n";
//
//  Initialize the package.
//
  cout << "\n";
  cout << "  INITIALIZE initializes the random number generator.\n";
  cout << "  It only needs to be called once before using the package.\n";

  initialize ( );
//
//  Set the generator index to 0.
//
  g = 0;
  cgn_set ( g );

  cout << "\n";
  cout << "  Current generator index = " << g << "\n";
//
//  Now call I4_UNI() several times.
//
  cout << "\n";
  cout << "   I     I4_UNI ( )\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++)
  {
    j = i4_uni ( );
    cout << "  " << setw(2) << i
         << "  " << setw(12) << j << "\n";
  }

  return;
}
//****************************************************************************80

void test02 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 calls R4_UNI_01 10 times, just to show how it is done.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int g;
  int i;
  float u;

  cout << "\n";
  cout << "TEST02\n";
  cout << "  R4_UNI_01 ( ) returns a random real number\n";
  cout << "  in [0,1] using the current generator.\n";
//
//  Initialize the package.
//
  cout << "\n";
  cout << "  INITIALIZE initializes the random number generator.\n";
  cout << "  It only needs to be called once before using the package.\n";

  initialize ( );
//
//  Set the current generator index.
//
  g = 0;
  cgn_set ( g );
  cout << "\n";
  cout << "  Current generator index = " << g << "\n";

  cout << "\n";
  cout << "   I     R4_UNI_01 ( )\n";
  cout << "\n";

  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << u << "\n";
  }

  return;
}
//****************************************************************************80

void test03 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST03 demonstrates how the seed can be reset to its initial or last value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int g;
  int i;
  float u;

  cout << "\n";
  cout << "TEST03\n";
  cout << "  R4_UNI_01 ( ) returns a random real number\n";
  cout << "  in [0,1] using the current generator.\n";
//
//  Initialize the package.
//
  cout << "\n";
  cout << "  INITIALIZE initializes the random number generator.\n";
  cout << "  It only needs to be called once before using the package.\n";

  initialize ( );
//
//  Set the current generator index.
//
  g = 16;
  cgn_set ( g );
  cout << "\n";
  cout << "  Current generator index = " << g << "\n";
//
//  Force the generator to start at the initial seed.
//
  cout << "\n";
  cout << "  INIT_GENERATOR can reset the seed to the initial value,\n";
  cout << "  the last (previous) value, or a new seed.\n";

  cout << "\n";
  cout << "  INIT_GENERATOR ( 0 ) starts at the initial seed.\n";

  init_generator ( 0 );

  cout << "\n";
  cout << "   I    R4_UNI_01 ( )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << u << "\n";
  }

  cout << "\n";
  cout << "  Calling INIT_GENERATOR ( 0 ) again restarts\n";
  cout << "  at the initial seed.\n";

  init_generator ( 0 );

  cout << "\n";
  cout << "   I    R4_UNI_01 ( )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << u << "\n";
  }

  cout << "\n";
  cout << "  Calling INIT_GENERATOR ( 2 ) restarts\n";
  cout << "  at a new \"far ahead\" seed.\n";

  init_generator ( 2 );

  cout << "\n";
  cout << "   I    R4_UNI_01 ( )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << u << "\n";
  }

  cout << "\n";
  cout << "  Calling INIT_GENERATOR ( 1 ) restarts\n";
  cout << "  at the last seed (in this case, the \"far ahead\"\n";
  cout << "  seed specified on the previous call.)\n";

  cout << "\n";
  cout << "   I    R4_UNI_01 ( )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    u = r4_uni_01 ( );
    cout << "  " << setw(2) << i
         << "  " << setw(14) << u << "\n";
    if ( ( i % 3 ) == 0 )
    {
      init_generator ( 1 );
      cout << "  (Reset to last seed)\n";
    }
  }

  return;
}
//****************************************************************************80

void test04 ( )

//****************************************************************************80
//
//  Purpose:
//
//    TEST04 demonstrates the use of multiple streams.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//     05 August 2013
//
//  Author:
//
//    John Burkardt
//
{
  int g[3];
  int i;
  int j;
  float u[3];

  cout << "\n";
  cout << "TEST04\n";
  cout << "  R4_UNI_01 ( ) returns a random real number\n";
  cout << "  in [0,1] using the current generator.\n";
//
//  Initialize the package.
//
  cout << "\n";
  cout << "  INITIALIZE initializes the random number generator.\n";
  cout << "  It only needs to be called once before using the package.\n";

  initialize ( );

  cout << "\n";
  cout << "  Let us call generators #2, #5 and #8.\n";
//
//  Use three separate generators.
//  Force them to start at their initial seeds.
//
  g[0] = 2;
  g[1] = 5;
  g[2] = 8;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  Initialize generator " << g[i] << "\n";
    cgn_set ( g[i] );
    init_generator ( 0 );
  }
//
//  Call the generators in the order 2, 5, 8.
//
  cout << "\n";
  cout << "   I    R4_UNI_01 ( 2 )  R4_UNI_01 ( 5 )  R4_UNI_01 ( 8 )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(2) << i;
    for ( j = 0; j < 3; j++ )
    {
      cgn_set ( g[j] );
      u[j] = r4_uni_01 ( );
      cout << "  " << setw(14) << u[j];
    }
    cout << "\n";
  }
//
//  Restart the generators at their initial seeds.
//
  g[0] = 5;
  g[1] = 8;
  g[2] = 2;
  cout << "\n";
  for ( i = 0; i < 3; i++ )
  {
    cout << "  Reinitialize generator " << g[i] << "\n";
    cgn_set ( g[i] );
    init_generator ( 0 );
  }
//
//  Call the generators in the order 5, 8, 2.
//
  cout << "\n";
  cout << "   I    R4_UNI_01 ( 5 )  R4_UNI_01 ( 8 )  R4_UNI_01 ( 2 )\n";
  cout << "\n";
  for ( i = 1; i <= 10; i++ )
  {
    cout << "  " << setw(2) << i;
    for ( j = 0; j < 3; j++ )
    {
      cgn_set ( g[j] );
      u[j] = r4_uni_01 ( );
      cout << "  " << setw(14) << u[j];
    }
    cout << "\n";
  }
  return;
}
