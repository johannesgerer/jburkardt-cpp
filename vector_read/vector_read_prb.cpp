# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>

using namespace std;

# include "vector_read.hpp"

int main ( void );
void test01 ( char *input_file_name );

//****************************************************************************80

int main ( void )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for VECTOR_READ_PRB.
//
//  Discussion:
//
//    VECTOR_READ_PRB tests routines from VECTOR_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  timestamp ( );

  cout << "\n";
  cout << "VECTOR_READ_PRB\n";
  cout << "  C++ version\n";
  cout << "  Test the VECTOR_READ library.\n";

  test01 ( "vector_read_prb.inp" );

  cout << "\n";
  cout << "VECTOR_READ_PRB:\n";
  cout << "  Normal end of execution.\n";

  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void test01 ( char *input_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests I4VEC_READ.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 March 2006
//
//  Author:
//
//    John Burkardt
//
{
  int i;
  ifstream input_file;
  int *i4vec;
  int i4vec_length;
  int i4vec_num;

  cout << "\n";
  cout << "TEST01\n";
  cout << "  I4VEC_READ reads integer vectors.\n";
  cout << "\n";

  i4vec_num = 0;

  input_file.open ( input_file_name );

  if ( !input_file )
  {
    cout << "\n";
    cout << "TEST01 - Fatal error!\n";
    cout << "  Could not open \"" << input_file_name << "\".\n";
    return;
  }

  for ( ; ; )
  {
    i4vec_length = i4vec_read ( input_file, &i4vec );

    if ( i4vec_length < 0 )
    {
      break;
    }

    i4vec_num = i4vec_num + 1;
    cout << "  " << setw(8) << i4vec_num 
         << "  " << setw(8) << i4vec_length;
    for ( i = 0; i < i4vec_length; i++ )
    {
      cout << "  " << i4vec[i];
    }
    cout << "\n";

    delete [] i4vec;
  }

  cout << "\n";
  cout << "  Number of vectors read was " << i4vec_num << "\n";

  input_file.close ( );

  return;
}
